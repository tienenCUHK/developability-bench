#!/usr/bin/env python3
"""Plug-and-play guidance for ProteinMPNN via --bias_by_res_jsonl.

ProteinMPNN's sampler adds the per-residue bias *after* dividing by temperature:

    probs = softmax(logits/T + bias_by_res/T)

so a bias value b produces an effective logit shift of b/T.  We therefore
parameterise guidance by the effective logit shift `lam` and store b = lam*T.

Two kinds of objective:
  * compositional (linear in amino-acid counts) -> the exact guidance is a
    position-independent per-AA tilt: b_a = -lam * w_a.  This is the closed-form
    "ProteinGuide for a linear predictor" case.
  * position-dependent (surface hydrophobicity) -> per-residue tilt, using the
    burial of the *input* backbone (CA neighbour count), which is available
    before any sequence is designed.

Usage:
  python guided_run.py --pdb data/1UBQ.pdb --property gravy --lam 2 --temp 0.1 --n 16 --out runs/g/1UBQ_gravy_2
"""
import argparse
import json
import os
import subprocess
import sys

AA21 = "ACDEFGHIKLMNPQRSTVWYX"
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2,
          I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9,
          Y=-1.3, V=4.2)
APOLAR = set("GPAVILMFYW")
CHARGE = {a: (1.0 if a in "KR" else (-1.0 if a in "DE" else 0.0)) for a in AA21}
MPNN = os.path.expanduser("~/invfold_survey/repos/ProteinMPNN/protein_mpnn_run.py")


def burial(pdb, chain="A", cutoff=10.0):
    xs = []
    for line in open(pdb, errors="ignore"):
        if line.startswith("ATOM") and line[21] == chain and line[12:16].strip() == "CA":
            xs.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    n = len(xs)
    out = []
    for i in range(n):
        c = 0
        for j in range(n):
            if i == j:
                continue
            dx, dy, dz = (xs[i][k] - xs[j][k] for k in range(3))
            if dx * dx + dy * dy + dz * dz <= cutoff * cutoff:
                c += 1
        out.append(c)
    return out


def build_bias(pdb, chain, prop, lam, temp, weights=None):
    L = len(burial(pdb, chain))
    bur = burial(pdb, chain)
    order = sorted(range(L), key=lambda i: bur[i])
    exposed = set(order[: int(0.5 * L)])  # the least-buried half
    wmap = {prop: lam}
    if weights:
        wmap.update(weights)
    bias = []
    for i in range(L):
        row = []
        for a in AA21:
            tot = 0.0
            for _p, _l in wmap.items():
                if _p == "agg":
                    # windowed exposure: penalise apolar residues in exposed 5-residue stretches
                    lo, hi = max(0, i - 2), min(L, i + 3)
                    w_i = sum(1 for j in range(lo, hi) if j in exposed) / (hi - lo)
                    tot += -_l * (1.0 if a in APOLAR else 0.0) * w_i
                elif _p == "gravy":
                    tot += -_l * KD.get(a, 0.0)
                elif _p == "charge":
                    tot += -_l * CHARGE.get(a, 0.0)
                elif _p == "surf":
                    tot += -_l * (1.0 if a in APOLAR else 0.0) * (1.0 if i in exposed else 0.0)
            row.append(tot * temp)
        bias.append(row)
    return bias


def _legacy_build_bias(pdb, chain, prop, lam, temp):
    L = len(burial(pdb, chain))
    bur = burial(pdb, chain)
    order = sorted(range(L), key=lambda i: bur[i])
    exposed = set(order[: int(0.5 * L)])
    bias = []
    for i in range(L):
        row = []
        for a in AA21:
            if prop == "gravy":
                v = -lam * KD.get(a, 0.0)                       # >0 lam => more hydrophilic
            elif prop == "charge":
                v = -lam * CHARGE.get(a, 0.0)                   # >0 lam => more negative
            elif prop == "surf":
                v = -lam * (1.0 if a in APOLAR else 0.0) * (1.0 if i in exposed else 0.0)
            else:
                raise SystemExit(f"unknown property {prop}")
            row.append(v * temp)      # sampler divides by T, so pre-multiply
        bias.append(row)
    return bias


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--chain", default="A")
    ap.add_argument("--property", default="surf")
    ap.add_argument("--lam", type=float, default=0.0)
    ap.add_argument("--weights", default="", help='e.g. "surf:10,charge:2" (overrides --property/--lam)')
    ap.add_argument("--temp", type=float, default=0.1)
    ap.add_argument("--n", type=int, default=16)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    name = os.path.basename(args.pdb).replace(".pdb", "")
    wm = {}
    for kv in [x for x in args.weights.split(",") if x]:
        k, v = kv.split(":")
        wm[k] = float(v)
    bias = build_bias(args.pdb, args.chain, args.property, args.lam, args.temp, wm or None)
    os.makedirs(args.out, exist_ok=True)
    jf = os.path.join(args.out, "bias.json")
    json.dump({name: {args.chain: bias}}, open(jf, "w"))

    cmd = [sys.executable, MPNN, "--pdb_path", args.pdb, "--pdb_path_chains", args.chain,
           "--num_seq_per_target", str(args.n), "--batch_size", str(args.n),
           "--sampling_temp", str(args.temp), "--out_folder", args.out,
           "--model_name", "v_48_020", "--seed", str(args.seed),
           "--bias_by_res_jsonl", jf]
    r = subprocess.run(cmd, capture_output=True, text=True)
    fa = os.path.join(args.out, "seqs", f"{name}.fa")
    if r.returncode != 0 or not os.path.exists(fa):
        print("FAILED", r.stderr[-400:], flush=True)
        sys.exit(1)
    print("OK", args.out, flush=True)


if __name__ == "__main__":
    main()
