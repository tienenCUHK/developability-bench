#!/usr/bin/env python3
"""p16 analyze: descriptors + recovery + diversity for every p16 design.

Adapted from analyze_bench30.py (same descriptor code, same diversity/recovery
definitions, same chain-A native reference); ROOT -> runs/bench30_p16 and output
-> results/p16/designs_analyzed.tsv.

Run with a Python that has torch and fair-esm installed (the descriptor pass imports the ESM stack).
"""
import collections
import csv
import glob
import itertools
import os
import statistics as st
import sys

ROOT = os.path.expanduser("~/invfold_survey")
sys.path.insert(0, ROOT)

from dev_bench import descriptors, ca_burial, read_fasta  # noqa: E402
from sd_eval import read_ca  # noqa: E402

RUNS = f"{ROOT}/runs/bench30_p16"
DATA = f"{ROOT}/data"
OUTDIR = f"{ROOT}/results/p16"
ARMS = ["vanilla", "soluble", "hyper", "vanilla_s1", "mem_global", "ligand",
        "g_surf10", "g_agg10", "g_charge3", "g_s10_cm2", "g_s10_cp2"]

native = {}
burial = {}


def nat(t):
    if t not in native:
        _, seq = read_ca(os.path.join(DATA, f"{t}.pdb"))
        native[t] = seq
        burial[t] = ca_burial(os.path.join(DATA, f"{t}.pdb"))
    return native[t], burial[t]


def main():
    rows = []
    for fa in sorted(glob.glob(os.path.join(RUNS, "*", "*", "seqs", "*.fa"))):
        parts = fa.split("/")
        arm, tgt = parts[-4], parts[-3]
        if arm not in ARMS:
            continue
        natseq, bur = nat(tgt)
        seqs = [(n, s) for n, s in read_fasta(fa)[1:]]
        if not seqs:
            print(f"skip empty {arm}/{tgt}", flush=True)
            continue
        L = min(len(s) for _, s in seqs)
        A = [list(s[:L]) for _, s in seqs]
        dv = (st.mean(sum(a != b for a, b in zip(A[i], A[j])) / L
                      for i, j in itertools.combinations(range(len(A)), 2))
              if len(A) > 1 else float("nan"))
        for name, s in seqs:
            rec = sum(1 for a, b in zip(s, natseq) if a == b) / max(1, min(len(s), len(natseq)))
            rows.append(dict(arm=arm, target=tgt, seq=s, recovery=round(rec, 4),
                             diversity=round(dv, 4), **descriptors(s, bur)))
        print(f"{arm:11s} {tgt:5s} n={len(seqs)}", flush=True)
    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, "designs_analyzed.tsv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("WROTE", out, len(rows), flush=True)


if __name__ == "__main__":
    main()
