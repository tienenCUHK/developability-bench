#!/usr/bin/env python3
"""p16 fold: resumable, work-stealing ESMFold fold of runs/bench30_p16.

Adapted from dev_fold.py (same ESMFold protocol: esmfold_v1, num_recycles=3,
scTM against chain-A CA of the input backbone, first fasta record = native).
Differences: output goes to results/p16/scTM_{tag}.tsv, and cells are claimed
with lock files so several processes (one per GPU 0-3) can share the whole
11-arm x 35-backbone work list dynamically.  Safe to re-launch after ssh drops:
a cell counts as done when its tsv rows cover every design in its fasta.

Usage:
  CUDA_VISIBLE_DEVICES=0 nohup ~/anaconda3/bin/python p16_fold.py --tag g0 ...
"""
import argparse
import glob
import os
import sys
import time
import traceback

ROOT = os.path.expanduser("~/invfold_survey")
RUNS = f"{ROOT}/runs/bench30_p16"
OUT = f"{ROOT}/results/p16"
LOCKS = f"{OUT}/locks"
DATA = f"{ROOT}/data"

for _p in ("~/invfold_survey/shims", "~/invfold_survey/repos/einops",
           "~/invfold_survey/repos/ml_collections", "~/invfold_survey/repos/openfold",
           "~/invfold_survey/repos/esm"):
    _p = os.path.expanduser(_p)
    if _p not in sys.path:
        sys.path.insert(0, _p)
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402
import torch  # noqa: E402
from tm import tmscore  # noqa: E402
from sd_eval import read_ca, read_fasta  # noqa: E402

ARMS = ["vanilla", "soluble", "hyper", "vanilla_s1", "mem_global", "ligand",
        "g_surf10", "g_agg10", "g_charge3", "g_s10_cm2", "g_s10_cp2"]


def cells():
    """[(arm, target, n_designs, length)] for every non-empty fasta."""
    out = []
    for arm in ARMS:
        for fa in sorted(glob.glob(f"{RUNS}/{arm}/*/seqs/*.fa")):
            parts = fa.split("/")
            tgt = parts[-3]
            recs = read_fasta(fa)[1:]
            if not recs:
                continue
            out.append((arm, tgt, len(recs), max(len(s) for _, s in recs)))
    return out


def load_done():
    """(arm,target) -> set(idx) already in any results/p16/scTM_*.tsv."""
    done = {}
    for f in glob.glob(f"{OUT}/scTM_*.tsv"):
        try:
            for line in open(f):
                p = line.rstrip("\n").split("\t")
                if len(p) != 4 or p[0] == "model":
                    continue
                done.setdefault((p[0], p[1]), set()).add(int(p[2]))
        except Exception:
            pass
    return done


def claim(key):
    path = os.path.join(LOCKS, key.replace("/", "__") + ".lock")
    for _ in range(2):
        try:
            fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.write(fd, f"{os.getpid()} {time.time():.0f}".encode())
            os.close(fd)
            return path
        except FileExistsError:
            try:
                age = time.time() - os.path.getmtime(path)
            except FileNotFoundError:
                continue
            if age > 1800:  # dead worker: stale lock
                try:
                    os.unlink(path)
                except OSError:
                    pass
                continue
            return None
        except OSError:
            return None
    return None


def release(path):
    try:
        os.unlink(path)
    except OSError:
        pass


def fold_cell(model, arm, tgt, recs):
    ref = f"{DATA}/{tgt}.pdb"
    ref_ca, _ = read_ca(ref)
    rows = []
    for i, (_, s) in enumerate(recs):
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
        rows.append((arm, tgt, i, sc))
    return rows


def append_rows(outfile, rows):
    new = not os.path.exists(outfile) or os.path.getsize(outfile) == 0
    with open(outfile, "a") as f:
        if new:
            f.write("model\ttarget\tidx\tscTM\n")
        for r in rows:
            f.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.4f}\n")
        f.flush()
        os.fsync(f.fileno())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="g0")
    ap.add_argument("--idle-sleep", type=float, default=10.0)
    ap.add_argument("--idle-rounds", type=int, default=60)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(LOCKS, exist_ok=True)
    outfile = f"{OUT}/scTM_{args.tag}.tsv"

    todo = cells()
    done = load_done()
    pending = [c for c in todo
               if len(done.get((c[0], c[1]), set())) < c[2]]
    bad = [c for c in todo
           if len(done.get((c[0], c[1]), set())) > c[2]]
    # biggest cells first so the long tail is spread over workers
    pending.sort(key=lambda c: -c[3] * c[2])
    print(f"P16_FOLD tag={args.tag} CUDA={os.environ.get('CUDA_VISIBLE_DEVICES')} "
          f"cells_total={len(todo)} pending={len(pending)} "
          f"already_done={len(todo)-len(pending)} odd_counts={bad[:5]}", flush=True)
    if not pending:
        print("P16_FOLD_DONE no pending cells", flush=True)
        return

    import esm.pretrained as P
    tload = time.time()
    model = P.esmfold_v1().eval().cuda()
    print(f"model loaded in {time.time()-tload:.0f}s", flush=True)

    n_ok = n_fail = 0
    failed = []
    idle = 0
    t00 = time.time()
    while True:
        if idle == 0 and n_ok + n_fail > 0:
            done = load_done()
        pending = [c for c in pending
                   if len(done.get((c[0], c[1]), set())) < c[2]]
        if not pending:
            print("P16_FOLD_DONE all cells folded", flush=True)
            break
        worked = False
        for cell in pending:
            arm, tgt, ndes, _L = cell
            lk = claim(f"{arm}/{tgt}")
            if lk is None:
                continue
            try:
                # re-check under the lock (another worker may have finished it)
                done = load_done()
                if len(done.get((arm, tgt), set())) >= ndes:
                    worked = True
                    break
                fa = f"{RUNS}/{arm}/{tgt}/seqs/{tgt}.fa"
                recs = read_fasta(fa)[1:]
                if not recs:
                    worked = True
                    break
                t0 = time.time()
                rows = fold_cell(model, arm, tgt, recs)
                append_rows(outfile, rows)
                d = np.mean([r[3] for r in rows])
                print(f"{arm:11s} {tgt:5s} n={len(recs)} scTM={d:.4f} "
                      f"sec={time.time()-t0:.0f}", flush=True)
                n_ok += 1
                worked = True
                idle = 0
                break
            except Exception:
                n_fail += 1
                failed.append(f"{arm}/{tgt}")
                print(f"FOLD_FAIL {arm}/{tgt}\n{traceback.format_exc()[-800:]}", flush=True)
                worked = True
                idle = 0
                break
            finally:
                release(lk)
        if worked:
            continue
        # nothing claimable now (everything left is locked by a live worker)
        idle += 1
        if idle > args.idle_rounds:
            print(f"P16_FOLD_DONE idle exit; remaining locked by others. "
                  f"ok={n_ok} fail={n_fail}", flush=True)
            break
        time.sleep(args.idle_sleep)
    print(f"P16_FOLD_END tag={args.tag} ok={n_ok} fail={n_fail} "
          f"failed={failed} elapsed={time.time()-t00:.0f}s", flush=True)


if __name__ == "__main__":
    main()
