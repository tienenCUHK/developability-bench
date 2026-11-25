#!/usr/bin/env python3
"""lam sweep fold: resumable, work-stealing ESMFold fold of runs/lam_sweep.

Same fold protocol as dev_fold.py (esmfold_v1, num_recycles=3, scTM against the
chain-A CA of the input backbone), but cells are claimed through lock files so
several workers (one per GPU 0-3) share the work list dynamically.  Writes
results/lam_sweep/scTM_{tag}.tsv and is safe to re-launch after ssh drops: a
cell counts as done when its rows cover every design in its fasta.

Usage:
  CUDA_VISIBLE_DEVICES=0 nohup ~/anaconda3/bin/python lam_fold.py --tag l0 \
      --arms surf2,surf5,surf20,surf10 > logs/lam_fold_l0.log 2>&1 &
"""
import argparse
import glob
import os
import sys
import time
import traceback

ROOT = os.path.expanduser("~/invfold_survey")
RUNS = f"{ROOT}/runs/lam_sweep"
OUT = f"{ROOT}/results/lam_sweep"
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


def cells(arms):
    out = []
    for fa in sorted(glob.glob(f"{RUNS}/*/*/seqs/*.fa")):
        parts = fa.split("/")
        arm, tgt = parts[-4], parts[-3]
        if arms and arm not in arms:
            continue
        recs = read_fasta(fa)[1:]
        if not recs:
            continue
        out.append((arm, tgt, len(recs), max(len(s) for _, s in recs)))
    return out


def load_done():
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
            if age > 1800:
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
    ref_ca, _ = read_ca(f"{DATA}/{tgt}.pdb")
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
    ap.add_argument("--tag", default="l0")
    ap.add_argument("--arms", default="")
    ap.add_argument("--idle-sleep", type=float, default=10.0)
    ap.add_argument("--idle-rounds", type=int, default=45)
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(LOCKS, exist_ok=True)
    outfile = f"{OUT}/scTM_{args.tag}.tsv"
    arms = args.arms.split(",") if args.arms else None

    todo = cells(arms)
    done = load_done()
    pending = [c for c in todo if len(done.get((c[0], c[1]), set())) < c[2]]
    pending.sort(key=lambda c: -c[3] * c[2])
    print(f"LAM_FOLD tag={args.tag} CUDA={os.environ.get('CUDA_VISIBLE_DEVICES')} "
          f"cells_total={len(todo)} pending={len(pending)}", flush=True)
    if not pending:
        print("LAM_FOLD_DONE no pending cells", flush=True)
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
        done = load_done()
        pending = [c for c in pending if len(done.get((c[0], c[1]), set())) < c[2]]
        if not pending:
            print("LAM_FOLD_DONE all cells folded", flush=True)
            break
        worked = False
        for arm, tgt, ndes, _L in pending:
            lk = claim(f"{arm}/{tgt}")
            if lk is None:
                continue
            try:
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
                print(f"{arm:11s} {tgt:5s} n={len(recs)} "
                      f"scTM={np.mean([r[3] for r in rows]):.4f} "
                      f"sec={time.time()-t0:.0f}", flush=True)
                n_ok += 1
                worked = idle = 0
                break
            except Exception:
                n_fail += 1
                failed.append(f"{arm}/{tgt}")
                print(f"FOLD_FAIL {arm}/{tgt}\n{traceback.format_exc()[-800:]}", flush=True)
                worked = idle = 0
                break
            finally:
                release(lk)
        if worked:
            continue
        idle += 1
        if idle > args.idle_rounds:
            print(f"LAM_FOLD_DONE idle exit; ok={n_ok} fail={n_fail}", flush=True)
            break
        time.sleep(args.idle_sleep)
    print(f"LAM_FOLD_END tag={args.tag} ok={n_ok} fail={n_fail} failed={failed} "
          f"elapsed={time.time()-t00:.0f}s", flush=True)


if __name__ == "__main__":
    main()
