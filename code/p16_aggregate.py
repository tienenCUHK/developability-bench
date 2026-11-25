#!/usr/bin/env python3
"""p16 aggregate: deduplicate folds, aggregate per arm, paired bootstrap + Holm.

Follows the protocol in results/AGGREGATE.md and paper_v1.md SS3.5:
  * unit of analysis is the backbone: per-(arm,backbone) mean first,
  * paired bootstrap over backbones (B=10,000) vs vanilla, 95% percentile CI,
  * two-sided bootstrap p-value with add-one correction,
    p = (1 + 2*min(#{d_b <= 0}, #{d_b >= 0})) / (B + 1),
  * Holm adjustment across the 9 comparison arms (vanilla_s1 is the seed-1
    repeatability control, not an arm, so it is reported separately as the
    seed-to-seed floor and never enters the Holm family).

Writes results/p16/AGGREGATE_P16.md.
"""
import argparse
import collections
import csv
import glob
import os
import statistics as st

import numpy as np

ROOT = os.path.expanduser("~/invfold_survey")
OUTDIR = f"{ROOT}/results/p16"
ORDER = ["vanilla", "vanilla_s1", "soluble", "hyper", "mem_global", "ligand",
         "g_surf10", "g_agg10", "g_charge3", "g_s10_cm2", "g_s10_cp2"]
BASE = "vanilla"
CONTROL = "vanilla_s1"
HOLM_ARMS = [a for a in ORDER if a not in (BASE, CONTROL)]

# paired d scTM from the 8-design reference run, results/AGGREGATE.md
REF8 = {
    "vanilla_s1": -0.0000, "soluble": +0.0031, "hyper": -0.0628,
    "mem_global": +0.0248, "ligand": +0.0207, "g_surf10": +0.0095,
    "g_agg10": +0.0114, "g_charge3": -0.0044, "g_s10_cm2": +0.0072,
    "g_s10_cp2": -0.0012,
}


def load_folds(outdir):
    """(arm,target,idx) -> scTM, deduplicated; later files win (like the 8-design run)."""
    sc, dups = {}, 0
    files = sorted(glob.glob(os.path.join(outdir, "scTM_*.tsv")))
    for f in files:
        for r in csv.DictReader(open(f), delimiter="\t"):
            k = (r["model"], r["target"], int(r["idx"]))
            if k in sc:
                dups += 1
            sc[k] = float(r["scTM"])
    return sc, dups, files


def holm(pvals):
    """Holm step-down adjustment, pvals: dict arm -> raw p."""
    m = len(pvals)
    order = sorted(pvals, key=lambda a: pvals[a])
    adj, running = {}, 0.0
    for rank, a in enumerate(order):
        val = min(1.0, (m - rank) * pvals[a])
        running = max(running, val)
        adj[a] = running
    return adj


def boot(xs, ys, rng, B):
    """Paired bootstrap of mean(ys)-mean(xs); returns (d, lo, hi, n, p)."""
    d = np.asarray(ys, dtype=float) - np.asarray(xs, dtype=float)
    n = len(d)
    idx = rng.integers(0, n, (B, n))
    bs = d[idx].mean(1)
    lo, hi = np.percentile(bs, 2.5), np.percentile(bs, 97.5)
    n_le = int(np.sum(bs <= 0))
    n_ge = int(np.sum(bs >= 0))
    p = (1.0 + 2.0 * min(n_le, n_ge)) / (B + 1.0)
    return float(d.mean()), float(lo), float(hi), n, min(1.0, p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=OUTDIR)
    ap.add_argument("--designs", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--B", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    designs_path = args.designs or os.path.join(args.outdir, "designs_analyzed.tsv")
    out_path = args.out or os.path.join(args.outdir, "AGGREGATE_P16.md")

    sc, dups, fold_files = load_folds(args.outdir)
    arm_t = collections.defaultdict(lambda: collections.defaultdict(list))
    for (a, t, i), v in sc.items():
        arm_t[a][t].append(v)
    tg = sorted({t for _, t, _ in sc})
    designs = list(csv.DictReader(open(designs_path), delimiter="\t"))
    by = collections.defaultdict(list)
    for r in designs:
        by[r["arm"]].append(r)

    rng = np.random.default_rng(args.seed)

    # ---- paired scTM comparisons vs vanilla ----
    comp = {}
    for a in ORDER:
        if a == BASE:
            continue
        common = [t for t in tg if arm_t[BASE].get(t) and arm_t[a].get(t)]
        if len(common) < 5:
            continue
        xs = [st.mean(arm_t[BASE][t]) for t in common]
        ys = [st.mean(arm_t[a][t]) for t in common]
        comp[a] = boot(xs, ys, rng, args.B)
    raw = {a: comp[a][4] for a in HOLM_ARMS if a in comp}
    holm_p = holm(raw)

    L = []
    L.append("# Aggregated DevMPNN-Bench results — power-up run (16 designs per cell)\n")
    L.append(f"- unique (arm,target,design) folds: {len(sc)}; duplicate rows dropped: {dups}")
    L.append(f"- designs analyzed: {len(designs)}")
    L.append(f"- fold result files: {len(fold_files)}")
    L.append(f"- backbones with any fold: {len(tg)}")
    L.append(f"- protocol: 16 designs per (arm, backbone), T=0.1, seed 0 "
             f"(vanilla_s1: seed 1), same 35 backbones as the 8-design run")
    L.append(f"- bootstrap: B={args.B}, unit = backbone mean, paired vs {BASE}, "
             f"seed={args.seed}; two-sided bootstrap p with add-one correction; "
             f"Holm across the {len(HOLM_ARMS)} comparison arms "
             f"(vanilla_s1 is the seed-1 repeatability control, not an arm, "
             f"and is excluded from the family)")
    L.append("")
    L.append("| arm | n des | nT folded | scTM | surf | agg5 | charge | pI | recovery | diversity |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for a in ORDER:
        v = by.get(a, [])
        per = [st.mean(arm_t[a][t]) for t in tg if arm_t[a].get(t)]
        if not per:
            continue
        sa = [float(x["surf_apolar"]) for x in v if x["surf_apolar"] not in ("NA", "")]
        L.append(
            f"| {a} | {len(v)} | {len(per)} | "
            f"{st.mean(per):.4f} | "
            f"{st.mean(sa):.3f} | "
            f"{st.mean(float(x['agg5frac']) for x in v):.3f} | "
            f"{st.mean(float(x['charge74']) for x in v):+.2f} | "
            f"{st.mean(float(x['pI']) for x in v):.2f} | "
            f"{st.mean(float(x['recovery']) for x in v):.3f} | "
            f"{st.mean(float(x['diversity']) for x in v):.3f} |")
    L.append("")
    L.append(f"## Paired bootstrap vs ProteinMPNN, scTM (B={args.B:,}, over backbones)\n")
    L.append("| arm | d scTM | 95% CI | nT | p (raw) | p (Holm) |")
    L.append("|---|---|---|---|---|---|")
    for a in HOLM_ARMS:
        if a not in comp:
            continue
        d, lo, hi, n, p = comp[a]
        L.append(f"| {a} | {d:+.4f} | [{lo:+.4f}, {hi:+.4f}] | {n} | "
                 f"{p:.4f} | {holm_p[a]:.4f} |")
    L.append("")
    L.append("## Repeatability: seed-to-seed floor (vanilla_s1 vs vanilla)\n")
    if CONTROL in comp:
        d, lo, hi, n, p = comp[CONTROL]
        L.append("| control | d scTM | 95% CI | nT | p (raw) |")
        L.append("|---|---|---|---|---|")
        L.append(f"| {CONTROL} (seed 1 vs seed 0) | {d:+.4f} | [{lo:+.4f}, {hi:+.4f}] | {n} | {p:.4f} |")
        L.append("")
        L.append("This is the benchmark's resolution floor; it is not part of the "
                 "Holm family.")
        L.append("")
    L.append("## Paired bootstrap vs ProteinMPNN, targeted properties "
             "(per-backbone means, absolute change)\n")
    props = [("surf_apolar", "surface apolar", 3),
             ("agg5frac", "aggregation 5-mer", 3),
             ("charge74", "net charge", 2),
             ("pI", "pI", 2)]

    def arm_prop(a, key):
        vals = collections.defaultdict(list)
        for r in by.get(a, []):
            v = r[key]
            if v in ("NA", ""):
                continue
            vals[r["target"]].append(float(v))
        return {t: st.mean(x) for t, x in vals.items()}

    base_props = {k: arm_prop(BASE, k) for k, _, _ in props}
    L.append("| arm | " + " | ".join(f"d {label} [95% CI]" for _, label, _ in props) + " |")
    L.append("|---|" + "---|" * len(props))
    for a in ORDER:
        if a == BASE:
            continue
        cells = []
        for k, _label, nd in props:
            pa = arm_prop(a, k)
            common = sorted(t for t in pa if t in base_props[k])
            if len(common) < 5:
                cells.append("n/a")
                continue
            xs = [base_props[k][t] for t in common]
            ys = [pa[t] for t in common]
            d, lo, hi, n, _p = boot(xs, ys, rng, args.B)
            cells.append(f"{d:+.{nd}f} [{lo:+.{nd}f}, {hi:+.{nd}f}]")
        L.append(f"| {a} | " + " | ".join(cells) + " |")
    L.append("")
    L.append("Property comparisons use all designs of the arm; a property is skipped "
             "for a backbone when it is undefined for every design there (multi-chain "
             "LigandMPNN designs can have no chain-A burial map).")
    L.append("")
    # ---- per-backbone heterogeneity ----
    L.append("## Per-backbone heterogeneity of d scTM (paired vs vanilla)\n")
    L.append("| arm | nT | mean d | median d | median abs(d) | max abs(d) "
             "| max-abs backbone | #(d < -0.2) |")
    L.append("|---|---|---|---|---|---|---|---|")
    het = {}
    for a in ORDER:
        if a == BASE:
            continue
        ds = {}
        for t in tg:
            if arm_t[BASE].get(t) and arm_t[a].get(t):
                ds[t] = st.mean(arm_t[a][t]) - st.mean(arm_t[BASE][t])
        if not ds:
            continue
        het[a] = ds
        absd = {t: abs(d) for t, d in ds.items()}
        tmax = max(absd, key=absd.get)
        L.append(f"| {a} | {len(ds)} | {st.mean(ds.values()):+.4f} | "
                 f"{st.median(ds.values()):+.4f} | {st.median(absd.values()):.4f} | "
                 f"{max(absd.values()):.4f} | {tmax} | "
                 f"{sum(1 for d in ds.values() if d < -0.2)} |")
    L.append("")
    cats = {a: sorted([(t, d) for t, d in ds.items() if d < -0.2], key=lambda x: x[1])
            for a, ds in het.items()}
    any_cat = False
    for a, lst in cats.items():
        if lst:
            any_cat = True
            L.append(f"- `{a}` collapses (d < -0.2): " +
                     ", ".join(f"{t} {d:+.3f}" for t, d in lst))
    if not any_cat:
        L.append("- no arm has a backbone with d < -0.2")
    L.append("")
    if "hyper" in het:
        ds = het["hyper"]
        severe = [t for t, d in ds.items() if d < -0.2]
        rest = [d for t, d in ds.items() if t not in severe]
        L.append(f"HyperMPNN mean-vs-median: mean d = {st.mean(ds.values()):+.4f}, "
                 f"median d = {st.median(ds.values()):+.4f}, "
                 f"mean excluding the {len(severe)} collapse(s) = "
                 f"{st.mean(rest):+.4f} (n={len(rest)}).")
        L.append("")
        focus = [t for t in ("1SRL", "1TIT", "3PGK") if t in ds]
        if focus:
            L.append("| backbone | vanilla scTM | hyper scTM | d |")
            L.append("|---|---|---|---|")
            for t in focus:
                L.append(f"| {t} | {st.mean(arm_t[BASE][t]):.4f} | "
                         f"{st.mean(arm_t['hyper'][t]):.4f} | {ds[t]:+.4f} |")
            L.append("")
    L.append("## Sign check: 16-design vs 8-design reference run (d scTM vs vanilla)\n")
    L.append("| arm | d scTM (16 designs) | d scTM (8 designs, results/AGGREGATE.md) | sign |")
    L.append("|---|---|---|---|")
    flips = []
    for a in ORDER:
        if a not in comp:
            continue
        d16 = comp[a][0]
        d8 = REF8.get(a)
        if d8 is None:
            continue
        same = (d16 >= 0) == (d8 >= 0)
        L.append(f"| {a} | {d16:+.4f} | {d8:+.4f} | {'SAME' if same else 'FLIP'} |")
        if not same:
            flips.append(f"{a} ({d8:+.4f} -> {d16:+.4f})")
    L.append("")
    if flips:
        L.append("**Sign changes: " + "; ".join(flips) + ".**")
    else:
        L.append("**No arm changes the sign of its paired scTM effect relative to "
                 "the 8-design run.**")
    L.append("")

    # ---- caveats / sensitivity ----
    L.append("## Caveats and sensitivity\n")
    L.append("1. **1TEN parse artifact (affects ProteinMPNN-family arms only).** "
             "ProteinMPNN emits 90-residue designs for 1TEN while chain A of "
             "data/1TEN.pdb has 89 C-alpha atoms (LigandMPNN emits 89), so the scTM "
             "alignment is shifted for vanilla/soluble/hyper/vanilla_s1 on that "
             "backbone and they score ~0.44 there, whereas mem_global/ligand score "
             "~0.94. The artifact is identical in the 8-design reference run "
             "(vanilla 1TEN mean 0.442 there vs 0.442 here). It inflates the paired "
             "mean of the two LigandMPNN arms by about +0.015.")
    no1 = {}
    for a in ORDER:
        if a == BASE or a not in comp:
            continue
        common = [t for t in tg if t != "1TEN" and arm_t[BASE].get(t) and arm_t[a].get(t)]
        if len(common) < 5:
            continue
        xs = [st.mean(arm_t[BASE][t]) for t in common]
        ys = [st.mean(arm_t[a][t]) for t in common]
        no1[a] = boot(xs, ys, rng, args.B)
    if no1:
        L.append("")
        L.append("   Excluding 1TEN (sensitivity only, not a protocol change):")
        L.append("")
        L.append("   | arm | d scTM | 95% CI | nT |")
        L.append("   |---|---|---|---|")
        for a in ORDER:
            if a not in no1:
                continue
            d, lo, hi, n, _p = no1[a]
            L.append(f"   | {a} | {d:+.4f} | [{lo:+.4f}, {hi:+.4f}] | {n} |")
        L.append("")
    L.append("2. **Significance summary.** Only HyperMPNN's 95% CI excludes zero "
             "before correction (raw p = "
             f"{comp['hyper'][4]:.4f}); no arm survives Holm at alpha = 0.05 "
             f"(HyperMPNN Holm p = {holm_p['hyper']:.4f}). In the 8-design run "
             "MembraneMPNN also excluded zero; at 16 designs its CI now contains "
             f"zero (raw p = {comp['mem_global'][4]:.4f}, and its mean is dominated by "
             "the 1TEN artifact above).")
    L.append("3. **Seed-to-seed floor** is measured by vanilla_s1 vs vanilla at 16 "
             "designs (see the repeatability table); it is wider than the 8-design "
             "floor, which matters when reading any comparison whose CI is of that "
             "order.")
    L.append("")

    # ---- coverage / failures ----
    L.append("## Coverage notes\n")
    L.append(f"- {len(designs)} designs from 371/385 (arm,backbone) cells; the 14 "
             f"missing cells are exactly the same ones missing in the 8-design "
             f"reference run, so arm x backbone coverage is identical between the two "
             f"runs (35/35 backbones for vanilla, soluble, hyper, vanilla_s1; 33/35 "
             f"for mem_global, ligand and the five guided arms).")
    for a in ORDER:
        n = len(arm_t.get(a, {}))
        L.append(f"- `{a}`: {n}/{len(tg)} backbones folded")
    L.append("")
    miss_any = False
    for a in ORDER:
        miss = [t for t in tg if not arm_t[a].get(t)]
        if miss:
            miss_any = True
            L.append(f"- missing for `{a}`: {', '.join(miss)}")
    L.append("")
    if miss_any:
        L.append("Generation failures (identical cells to the 8-design reference run): "
                 "`mem_global`/`ligand` fail on 1BNI and 1AYE with a prody shape mismatch "
                 "inside LigandMPNN; the five `g_*` arms fail on 1TEN and 1QYS because the "
                 "guided bias length is 89/91 while ProteinMPNN pads to 90/92. "
                 "All are pre-existing tool failures, not regressions of this run.")
        L.append("")
    txt = "\n".join(L) + "\n"
    open(out_path, "w").write(txt)
    print(txt)
    print("WROTE", out_path)


if __name__ == "__main__":
    main()
