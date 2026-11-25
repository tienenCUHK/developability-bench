#!/usr/bin/env python3
"""lam_analyze: lambda sweep + global-vs-positional ablation for the surface arm.

Deliverable 1: per-lambda table (designs, folded backbones, scTM, surface apolar,
aggregation 5-mer, charge, pI, diversity) and the per-backbone OLS slope k of
achieved property on lambda (surface and scTM), with a backbone bootstrap CI.

Deliverable 2: at the best lambda, compare the exposure-masked arm with a
global (all-position) apolar penalty arm, including paired delta scTM / delta
surface apolar with backbone bootstrap CIs.

Folds: results/p16/scTM_*.tsv (vanilla + g_surf10 reused as surf10) and
       results/lam_sweep/scTM_*.tsv (new arms).
Writes results/lam_sweep/designs_analyzed.tsv and results/lam_sweep/AGGREGATE_LAM.md.
"""
import argparse
import collections
import csv
import glob
import os
import statistics as st
import sys

ROOT = os.path.expanduser("~/invfold_survey")
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from lam_gen import (chain_a_positions, burial_of, read_fasta,  # noqa: E402
                     APOLAR)
from dev_bench import net_charge, isoelectric  # noqa: E402

P16 = f"{ROOT}/results/p16"
LAM = f"{ROOT}/results/lam_sweep"
RUNS_LAM = f"{ROOT}/runs/lam_sweep"
RUNS_P16 = f"{ROOT}/runs/bench30_p16"
DATA = f"{ROOT}/data"

# arm name -> (mode, lambda)
ARM_CFG = {
    "surf2": ("surf", 2.0),
    "surf5": ("surf", 5.0),
    "surf10": ("surf", 10.0),
    "surf20": ("surf", 20.0),
}
LAMS = [2.0, 5.0, 10.0, 20.0]
BASE_FOLD_SOURCE = {"g_surf10": "surf10"}
BASELINE = "vanilla"


# ---------------------------------------------------------------- properties
def cell_props(pdb, seqs):
    """Per-design descriptors with burial aligned to ProteinMPNN's chain order."""
    pos = chain_a_positions(pdb)
    native = "".join(aa for aa, _ in pos)
    known = [(i, ca) for i, (_, ca) in enumerate(pos) if ca is not None]
    coords = [ca for _, ca in known]
    bur = burial_of(coords)
    bur_pos = [None] * len(pos)
    for (i, _), b in zip(known, bur):
        bur_pos[i] = b
    known_b = sorted(bur)
    med = known_b[len(known_b) // 2] if known_b else None
    exposed = [i for i, b in enumerate(bur_pos) if b is not None and b <= med]
    rows = []
    L = min(len(s) for s in seqs)
    for s in seqs:
        Ls = len(s)
        rec = (sum(1 for a, b in zip(s, native) if a == b) / max(1, min(Ls, len(native))))
        win = [sum(1 for a in s[i:i + 5] if a in APOLAR) / 5
               for i in range(max(0, Ls - 4))]
        agg5 = sum(1 for w in win if w >= 0.8) / max(1, len(win))
        surf = "NA"
        if exposed:
            surf = round(sum(1 for i in exposed if i < Ls and s[i] in APOLAR) / len(exposed), 3)
        rows.append(dict(length=Ls, charge74=round(net_charge(s, 7.4), 2),
                         pI=round(isoelectric(s), 2), surf_apolar=surf,
                         agg5frac=round(agg5, 3), recovery=round(rec, 4)))
    return pos, L, rows


def analyze_fastas(sources):
    """sources: list of (root, arm_glob, out_arm, lam). Returns designs rows."""
    out = []
    for root, arm_dir, out_arm, lam in sources:
        for fa in sorted(glob.glob(f"{root}/{arm_dir}/*/seqs/*.fa")):
            tgt = fa.split("/")[-3]
            recs = read_fasta(fa)[1:]
            if not recs:
                continue
            seqs = [s for _, s in recs]
            pos, L, rows = cell_props(f"{DATA}/{tgt}.pdb", seqs)
            if L < 2:
                dv = float("nan")
            else:
                A = [list(s[:L]) for s in seqs]
                dv = st.mean(sum(a != b for a, b in zip(A[i], A[j])) / L
                             for i, j in itertools_comb(len(A)))
            for s, r in zip(seqs, rows):
                out.append(dict(arm=out_arm, lam=lam, target=tgt, seq=s,
                                diversity=round(dv, 4), **r))
    return out


def itertools_comb(n):
    for i in range(n):
        for j in range(i + 1, n):
            yield i, j


# ---------------------------------------------------------------- fold loading
def load_folds_dedup():
    seen = set()
    sc = collections.defaultdict(dict)
    def add(root, mapping):
        for f in sorted(glob.glob(f"{root}/scTM_*.tsv")):
            for r in csv.DictReader(open(f), delimiter="\t"):
                arm = mapping.get(r["model"], r["model"])
                k = (arm, r["target"], int(r["idx"]))
                if k in seen:
                    continue
                seen.add(k)
                sc[arm].setdefault(r["target"], []).append(float(r["scTM"]))
    add(P16, {**BASE_FOLD_SOURCE, **{a: a for a in ARM_CFG}})
    add(LAM, {a: a for a in ARM_CFG})
    return sc


# ---------------------------------------------------------------- bootstrap
def boot_mean(vals, rng, B):
    vals = np.asarray(vals, dtype=float)
    n = len(vals)
    idx = rng.integers(0, n, (B, n))
    bs = vals[idx].mean(1)
    return float(vals.mean()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)), n


def boot_paired(a, b, rng, B):
    """bootstrap of mean(b - a) over paired backbones."""
    d = np.asarray(b, dtype=float) - np.asarray(a, dtype=float)
    n = len(d)
    idx = rng.integers(0, n, (B, n))
    bs = d[idx].mean(1)
    n_le = int((bs <= 0).sum())
    n_ge = int((bs >= 0).sum())
    p = min(1.0, (1.0 + 2.0 * min(n_le, n_ge)) / (B + 1.0))
    return float(d.mean()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)), n, p


def ols_slope(xs, ys):
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    if len(xs) < 2:
        return None
    b, a = np.polyfit(xs, ys, 1)
    pred = a + b * xs
    ss_res = float(((ys - pred) ** 2).sum())
    ss_tot = float(((ys - ys.mean()) ** 2).sum())
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(b), r2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--surfglob-lam", type=float, default=None,
                    help="lambda of the global-vs-positional ablation arm "
                         "(default: auto-detect surfglobN arm directories)")
    args = ap.parse_args()
    os.makedirs(LAM, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    # ---- which surfglob lambda is present? ----
    glob_arms = {}
    for d in glob.glob(f"{RUNS_LAM}/surfglob*"):
        name = os.path.basename(d)
        try:
            lam = float(name.replace("surfglob", ""))
        except ValueError:
            continue
        glob_arms[name] = lam
    if args.surfglob_lam is not None:
        glob_arms = {f"surfglob{args.surfglob_lam:g}": args.surfglob_lam}

    # ---- designs (properties) ----
    cfg = dict(ARM_CFG)
    for name, lam in glob_arms.items():
        cfg[name] = ("surfglob", lam)
    sources = []
    for arm, (mode, lam) in ARM_CFG.items():
        if arm == "surf10":
            # reuse p16 g_surf10 (33 backbones) + newly generated 1TEN/1QYS
            sources.append((RUNS_P16, "g_surf10", "surf10", 10.0))
        else:
            sources.append((RUNS_LAM, arm, arm, lam))
    for name in glob_arms:
        sources.append((RUNS_LAM, name, name, glob_arms[name]))
    designs = analyze_fastas(sources)
    # add the 2 new surf10 cells (lam root) only once
    existing = {(r["arm"], r["target"]) for r in designs}
    extra = analyze_fastas([(RUNS_LAM, "surf10", "surf10", 10.0)])
    designs += [r for r in extra if (r["arm"], r["target"]) not in existing]
    # baseline properties from p16 (vanilla)
    base_designs = list(csv.DictReader(open(f"{P16}/designs_analyzed.tsv"), delimiter="\t"))
    base_designs = [r for r in base_designs if r["arm"] == BASELINE]
    fieldnames = list(designs[0].keys())
    with open(f"{LAM}/designs_analyzed.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(designs)
    by = collections.defaultdict(list)
    for r in designs:
        by[r["arm"]].append(r)

    # ---- folds ----
    sc = load_folds_dedup()
    arm_t = {a: {t: st.mean(v) for t, v in d.items()} for a, d in sc.items()}
    base_t = arm_t.get(BASELINE, {})

    # ---- per-lambda table ----
    rows = []
    for lam in LAMS:
        arm = f"surf{lam:g}"
        if arm not in by:
            continue
        v = by[arm]
        per = [arm_t[arm][t] for t in sorted(arm_t.get(arm, {}))]
        if not per:
            continue
        sa = [float(x["surf_apolar"]) for x in v if x["surf_apolar"] != "NA"]
        common = [t for t in sorted(base_t) if arm_t[arm].get(t)]
        if common:
            d_ci = boot_paired([base_t[t] for t in common],
                               [arm_t[arm][t] for t in common], rng, args.B)
        else:
            d_ci = (float("nan"),) * 5
        rows.append(dict(lam=lam, arm=arm, n_des=len(v), nT=len(per),
                         scTM=st.mean(per), surf=st.mean(sa),
                         agg5=st.mean(float(x["agg5frac"]) for x in v),
                         charge=st.mean(float(x["charge74"]) for x in v),
                         pI=st.mean(float(x["pI"]) for x in v),
                         diversity=st.mean(float(x["diversity"]) for x in v),
                         recovery=st.mean(float(x["recovery"]) for x in v),
                         d_scTM=d_ci[0], d_lo=d_ci[1], d_hi=d_ci[2], d_n=d_ci[3]))

    # ---- k: per-backbone OLS slope of property on lambda ----
    # achieved surface apolar per (lambda, backbone) from designs
    surf_bt = collections.defaultdict(dict)
    for r in designs:
        if r["surf_apolar"] == "NA":
            continue
        surf_bt[r["arm"]].setdefault(r["target"], []).append(float(r["surf_apolar"]))
    surf_bt = {a: {t: st.mean(v) for t, v in d.items()} for a, d in surf_bt.items()}
    ks = {"surf": [], "scTM": []}
    r2s = {"surf": [], "scTM": []}
    used_backbones = []
    for t in sorted(base_t):
        lams_ok, surf_vals, sc_vals = [], [], []
        ok = True
        for lam in LAMS:
            arm = f"surf{lam:g}"
            if t not in surf_bt.get(arm, {}) or t not in arm_t.get(arm, {}):
                ok = False
                break
            lams_ok.append(lam)
            surf_vals.append(surf_bt[arm][t])
            sc_vals.append(arm_t[arm][t])
        if not ok or len(lams_ok) < 2:
            continue
        used_backbones.append(t)
        for key, vals in (("surf", surf_vals), ("scTM", sc_vals)):
            s, r2 = ols_slope(lams_ok, vals)
            ks[key].append(s)
            r2s[key].append(r2)
    kstats = {}
    for key in ("surf", "scTM"):
        if ks[key]:
            m, lo, hi, n = boot_mean(ks[key], rng, args.B)
            d = np.asarray(ks[key])
            idx = rng.integers(0, len(d), (args.B, len(d)))
            bs = d[idx].mean(1)
            p = min(1.0, (1 + 2 * min(int((bs <= 0).sum()), int((bs >= 0).sum()))) / (args.B + 1))
            kstats[key] = dict(mean=m, lo=lo, hi=hi, n=n, p=p,
                               median=float(np.median(d)), r2=float(np.nanmean(r2s[key])))

    # ---- best lambda (rule: strongest control whose paired d scTM CI lower
    #      bound is >= -0.02, i.e. no meaningful foldability loss) ----
    best = None
    for r in sorted(rows, key=lambda x: -x["lam"]):
        if r["d_n"] >= 5 and r["d_lo"] >= -0.02:
            best = r
            break
    if best is None and rows:
        best = max(rows, key=lambda x: x["lam"])
    best_lam = best["lam"] if best else None

    # ---- deliverable 2: masked vs global at best lambda ----
    glob_row = None
    if best_lam is not None:
        gname = f"surfglob{best_lam:g}"
        mask_arm = f"surf{best_lam:g}"
        if gname in by and mask_arm in by:
            gv = by[gname]
            mper = [arm_t[mask_arm][t] for t in sorted(arm_t.get(mask_arm, {}))]
            gper = [arm_t[gname][t] for t in sorted(arm_t.get(gname, {}))]
            common = [t for t in sorted(arm_t.get(gname, {}))
                      if arm_t.get(mask_arm, {}).get(t)]
            if mper and gper and common:
                dsc = boot_paired([arm_t[mask_arm][t] for t in common],
                                  [arm_t[gname][t] for t in common], rng, args.B)
                ms_bt = collections.defaultdict(list)
                for r in by[mask_arm]:
                    if r["surf_apolar"] != "NA":
                        ms_bt[r["target"]].append(float(r["surf_apolar"]))
                gs_bt = collections.defaultdict(list)
                for r in gv:
                    if r["surf_apolar"] != "NA":
                        gs_bt[r["target"]].append(float(r["surf_apolar"]))
                csurf = [t for t in sorted(gs_bt) if t in ms_bt]
                if csurf:
                    dsurf = boot_paired([st.mean(ms_bt[t]) for t in csurf],
                                        [st.mean(gs_bt[t]) for t in csurf],
                                        rng, args.B)
                    msurf = [float(x["surf_apolar"]) for x in by[mask_arm]
                             if x["surf_apolar"] != "NA"]
                    gsurf = [float(x["surf_apolar"]) for x in gv
                             if x["surf_apolar"] != "NA"]
                    glob_row = dict(lam=best_lam, arm=gname, n_des=len(gv),
                                    nT=len(gper), scTM=st.mean(gper),
                                    surf=st.mean(gsurf),
                                    agg5=st.mean(float(x["agg5frac"]) for x in gv),
                                    charge=st.mean(float(x["charge74"]) for x in gv),
                                    pI=st.mean(float(x["pI"]) for x in gv),
                                    diversity=st.mean(float(x["diversity"]) for x in gv),
                                    recovery=st.mean(float(x["recovery"]) for x in gv),
                                    masked_surf=st.mean(msurf),
                                    masked_scTM=st.mean(mper),
                                    d_scTM=dsc, d_surf=dsurf)

    # ---- markdown ----
    L = []
    L.append("# Lambda sweep + global-vs-positional ablation — surface objective\n")
    L.append(f"- protocol: 16 designs per (lambda, backbone), ProteinMPNN v_48_020, "
             f"T=0.1, seed 0, 35 backbones; T=0.1; scTM = ESMFold refold vs chain-A "
             f"backbone (same folds/protocol as DevMPNN-Bench)")
    L.append(f"- lambda=10 reuses the 16-design g_surf10 cells from the main p16 run "
             f"(33 backbones) plus newly sampled 1TEN/1QYS, so all lambda values cover "
             f"35 backbones")
    L.append(f"- bootstrap: B={args.B}, unit = backbone, seed={args.seed}; "
             f"95% percentile CI; two-sided bootstrap p with add-one correction")
    L.append("")
    L.append("## Deliverable 1 — surface arm, lambda sweep (35 backbones)\n")
    L.append("| lambda | n des | nT folded | scTM | surf apolar | agg5 | charge | pI | diversity | recovery | d scTM vs vanilla [95% CI] |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        ci = (f"{r['d_scTM']:+.4f} [{r['d_lo']:+.4f}, {r['d_hi']:+.4f}]"
              if r["d_n"] >= 5 else "n/a")
        L.append(f"| {r['lam']:g} | {r['n_des']} | {r['nT']} | {r['scTM']:.4f} | "
                 f"{r['surf']:.3f} | {r['agg5']:.3f} | {r['charge']:+.2f} | "
                 f"{r['pI']:.2f} | {r['diversity']:.3f} | {r['recovery']:.3f} | {ci} |")
    L.append("")
    L.append("## Per-backbone controllability slope k (achieved property vs lambda)\n")
    if kstats:
        L.append(f"Per backbone: OLS slope over lambda in {LAMS}; estimate = mean over "
                 f"n={kstats['surf']['n']} backbones with bootstrap CI.\n")
        L.append("| quantity | k (per lambda unit) | 95% CI | median per-backbone k | mean R^2 | n backbones | p (k=0) |")
        L.append("|---|---|---|---|---|---|---|")
        for key, label in (("surf", "surface apolar"), ("scTM", "scTM")):
            s = kstats[key]
            L.append(f"| {label} | {s['mean']:+.5f} | [{s['lo']:+.5f}, {s['hi']:+.5f}] | "
                     f"{s['median']:+.5f} | {s['r2']:.3f} | {s['n']} | {s['p']:.4f} |")
        L.append("")
        ksurf = kstats["surf"]["mean"]
        ksc = kstats["scTM"]["mean"]
        L.append(f"Direction: surface apolar changes {ksurf:+.5f} per unit lambda "
                 f"(negative = the guidance moves the achieved property in the intended "
                 f"direction); scTM changes {ksc:+.5f} per unit lambda.")
    else:
        L.append("Not enough data to fit k.")
    L.append("")
    L.append("## Deliverable 2 — exposure mask vs global apolar penalty\n")
    if glob_row is not None:
        L.append(f"Both arms at lambda = {glob_row['lam']:g}; paired over backbones.\n")
        L.append("| arm | n des | nT folded | scTM | surf apolar | agg5 | charge | pI | diversity |")
        L.append("|---|---|---|---|---|---|---|---|---|")
        br = best
        L.append(f"| mask (surf:{best_lam:g}) | {br['n_des']} | {br['nT']} | "
                 f"{br['scTM']:.4f} | {br['surf']:.3f} | {br['agg5']:.3f} | "
                 f"{br['charge']:+.2f} | {br['pI']:.2f} | {br['diversity']:.3f} |")
        g = glob_row
        L.append(f"| global (surfglob:{g['lam']:g}) | {g['n_des']} | {g['nT']} | "
                 f"{g['scTM']:.4f} | {g['surf']:.3f} | {g['agg5']:.3f} | "
                 f"{g['charge']:+.2f} | {g['pI']:.2f} | {g['diversity']:.3f} |")
        L.append("")
        d, lo, hi, n, p = g["d_scTM"]
        L.append(f"- paired **d scTM (global - mask)** = {d:+.4f} [{lo:+.4f}, {hi:+.4f}], "
                 f"nT={n}, raw p={p:.4f}")
        d, lo, hi, n, p = g["d_surf"]
        L.append(f"- paired **d surface apolar (global - mask)** = {d:+.4f} "
                 f"[{lo:+.4f}, {hi:+.4f}], nT={n}, raw p={p:.4f} "
                 f"(negative = global controls better)")
    else:
        L.append("No global-penalty arm found yet (surfglob<lambda> directory missing).")
    L.append("")
    # one-line answer
    if kstats and glob_row is not None and best is not None:
        d_dsc = glob_row["d_scTM"]
        d_dsurf = glob_row["d_surf"]
        fold_txt = ("Foldability does not degrade with lambda" if abs(ksc) < 0.002
                    else "Foldability degrades with lambda")
        d_agg = glob_row["agg5"] - best["agg5"]
        if abs(d_dsurf[0]) < 0.01 and d_dsc[0] > -0.01:
            mask_txt = (f"the exposure mask does not matter at this lambda "
                        f"(global - mask: d scTM = {d_dsc[0]:+.4f}, "
                        f"d surface = {d_dsurf[0]:+.4f})")
        else:
            mask_txt = (f"the exposure mask does real work: dropping it (global penalty) "
                        f"additionally lowers surface apolar by {-d_dsurf[0]:.4f} and "
                        f"aggregation 5-mer by {-d_agg:.4f} but costs {d_dsc[0]:+.4f} scTM "
                        f"(95% CI [{d_dsc[1]:+.4f}, {d_dsc[2]:+.4f}], p = {d_dsc[4]:.4f})")
        L.append("## One-line answer\n")
        L.append(f"{fold_txt} (k_scTM = {ksc:+.5f} per lambda unit, 95% CI "
                 f"[{kstats['scTM']['lo']:+.5f}, {kstats['scTM']['hi']:+.5f}]); at the "
                 f"best lambda = {best_lam:g}, {mask_txt}.")
        L.append("")
    if best_lam is not None:
        L.append(f"## Best lambda rule\n\nLargest lambda whose paired d scTM vs vanilla "
                 f"has 95% CI lower bound >= -0.02 (no meaningful foldability cost); "
                 f"selected lambda = **{best_lam:g}**.\n")
    else:
        L.append("## Best lambda rule\n\nNot determinable yet (folds incomplete).\n")
    txt = "\n".join(L)
    open(f"{LAM}/AGGREGATE_LAM.md", "w").write(txt)
    print(txt)
    print("WROTE", f"{LAM}/AGGREGATE_LAM.md")


if __name__ == "__main__":
    main()
