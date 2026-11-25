#!/usr/bin/env python3
"""Fold every design of the DevMPNN benchmark and report scTM per (model, target)."""
import csv
import glob
import os
import sys
import time

import numpy as np

for _p in ("~/invfold_survey/shims", "~/invfold_survey/repos/einops",
           "~/invfold_survey/repos/ml_collections", "~/invfold_survey/repos/openfold",
           "~/invfold_survey/repos/esm"):
    _p = os.path.expanduser(_p)
    if _p not in sys.path:
        sys.path.insert(0, _p)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import torch  # noqa: E402
from tm import tmscore  # noqa: E402
from sd_eval import read_ca, read_fasta, diversity  # noqa: E402

ROOT = os.path.expanduser(os.environ.get("FOLD_ROOT", "~/invfold_survey/runs/dev"))
DATA = os.path.expanduser("~/invfold_survey/data")


def main():
    models = sys.argv[1].split(",") if len(sys.argv) > 1 else None
    import esm.pretrained as P
    tag = sys.argv[2] if len(sys.argv) > 2 else "all"
    outfile = os.path.join(ROOT, f"scTM_{tag}.tsv")
    rows, done = [], set()
    if os.path.exists(outfile):
        for r in csv.DictReader(open(outfile), delimiter="\t"):
            rows.append((r["model"], r["target"], int(r["idx"]), float(r["scTM"])))
            done.add((r["model"], r["target"]))
        print(f"resume: {len(done)} (model,target) cells already folded", flush=True)
    else:
        with open(outfile, "w") as f:
            f.write("model\ttarget\tidx\tscTM\n")
    model = P.esmfold_v1().eval().cuda()
    for fa in sorted(glob.glob(os.path.join(ROOT, "*", "*", "seqs", "*.fa"))):
        parts = fa.split("/")
        m, tgt = parts[-4], parts[-3]
        if models and m not in models:
            continue
        if (m, tgt) in done:
            continue
        ref = os.path.join(DATA, f"{tgt}.pdb")
        if not os.path.exists(ref):
            continue
        ref_ca, _ = read_ca(ref)
        if len(ref_ca) == 0:
            print(f"skip {m}/{tgt}: no chain-A CA atoms", flush=True)
            continue
        recs = read_fasta(fa)[1:]
        if not recs:
            continue
        D = diversity([s for _, s in recs])
        t0 = time.time()
        for i, (name, s) in enumerate(recs):
            with torch.no_grad():
                out = model.infer([s], num_recycles=3)
            pos = out["positions"].float().cpu().numpy()
            while pos.ndim > 4:
                pos = pos[0]
            if pos.ndim == 4:
                pos = pos[0]
            ca = pos[:, 1, :]
            n = min(len(ca), len(ref_ca))
            sc = tmscore(ca[:n].astype(np.float64), ref_ca[:n], L=len(ref_ca))
            rows.append((m, tgt, i, sc))
        print(f"{m:12s} {tgt:6s} n={len(recs)} D={D:.3f} scTM={np.mean([r[3] for r in rows[-len(recs):]]):.4f} "
              f"sec={time.time()-t0:.0f}", flush=True)
        with open(outfile, "a") as f:
            for r in rows[-len(recs):]:
                f.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.4f}\n")
    print("DONE", outfile, len(rows), flush=True)


if __name__ == "__main__":
    main()
