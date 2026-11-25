#!/usr/bin/env python3
"""p16 generation: 16 designs per (arm, backbone) for DevMPNN-Bench.

Scales the 8-design protocol (bench30_grid.py / agg_grid.py / hyper_grid.py /
seed2_grid.py) to 16 designs per cell, changing ONLY n=8 -> n=16.  Same weights,
same temperature 0.1, same seeds (0; 1 for vanilla_s1), same 35 backbones.

Runs are sharded round-robin over cells so four processes (one per GPU 0-3) can
run concurrently.  Idempotent: a cell whose fasta already holds the native
record + 16 designs is skipped, so the launcher can be re-run after ssh drops.

Usage (per GPU shard, always CUDA_VISIBLE_DEVICES in 0..3 only):
  CUDA_VISIBLE_DEVICES=0 ~/anaconda3/bin/python p16_gen.py --shard 0
  ...
  CUDA_VISIBLE_DEVICES=3 ~/anaconda3/bin/python p16_gen.py --shard 3
"""
import argparse
import os
import subprocess
import sys
import time

ROOT = os.path.expanduser("~/invfold_survey")
PY_ = sys.executable                         # must be a torch env python
LIG_PY = os.path.expanduser("~/anaconda3/envs/antifold/bin/python")  # has prody
MPNN = f"{ROOT}/repos/ProteinMPNN/protein_mpnn_run.py"
LIG = f"{ROOT}/repos/LigandMPNN/run.py"
LP = f"{ROOT}/repos/LigandMPNN/model_params"
GUIDED = f"{ROOT}/guided_run.py"
RUNS = f"{ROOT}/runs/bench30_p16"

N = 16
TARGETS = [ln.strip() for ln in open(f"{ROOT}/data/targets30.txt") if ln.strip()]

# arm -> (kind, spec); mirrors the reference 8-design grid scripts exactly
ARMS = {
    "vanilla":    ("pmpnn", "vanilla_model_weights", 0),
    "soluble":    ("pmpnn", "soluble_model_weights", 0),
    "hyper":      ("hyper", "", 0),
    "vanilla_s1": ("pmpnn", "vanilla_model_weights", 1),
    "mem_global": ("lig", "global_label_membrane_mpnn", 0),
    "ligand":     ("lig", "ligand_mpnn", 0),
    "g_surf10":   ("guide", "surf:10", 0),
    "g_agg10":    ("guide", "agg:10", 0),
    "g_charge3":  ("guide", "charge:3", 0),
    "g_s10_cm2":  ("guide", "surf:10,charge:-2", 0),
    "g_s10_cp2":  ("guide", "surf:10,charge:2", 0),
}
ORDER = ["vanilla", "soluble", "hyper", "vanilla_s1", "mem_global", "ligand",
         "g_surf10", "g_agg10", "g_charge3", "g_s10_cm2", "g_s10_cp2"]


def read_fasta(path):
    recs, name, buf = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name:
                recs.append((name, "".join(buf)))
            name, buf = line[1:], []
        elif line:
            buf.append(line)
    if name:
        recs.append((name, "".join(buf)))
    return recs


def cell_done(fa):
    """True iff the fasta has the native record + N designs."""
    if not os.path.exists(fa):
        return False
    try:
        return len(read_fasta(fa)) >= N + 1
    except Exception:
        return False


def run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True,
                          env={**os.environ})


def generate(arm, target):
    kind, spec, seed = ARMS[arm]
    out = f"{RUNS}/{arm}/{target}"
    fa = f"{out}/seqs/{target}.fa"
    if cell_done(fa):
        return "skip", None
    os.makedirs(out, exist_ok=True)
    pdb = f"{ROOT}/data/{target}.pdb"
    if kind in ("pmpnn", "hyper"):
        wdir = (f"{ROOT}/repos/HyperMPNN_weights" if kind == "hyper"
                else f"{ROOT}/repos/ProteinMPNN/{spec}")
        mname = "v48_020_epoch300_hyper" if kind == "hyper" else "v_48_020"
        cmd = [PY_, MPNN, "--pdb_path", pdb, "--pdb_path_chains", "A",
               "--num_seq_per_target", str(N), "--batch_size", str(N),
               "--sampling_temp", "0.1", "--out_folder", out,
               "--model_name", mname, "--seed", str(seed),
               "--path_to_model_weights", wdir]
        r = run(cmd)
    elif kind == "lig":
        ck = {"global_label_membrane_mpnn": "global_label_membrane_mpnn_v_48_020.pt",
              "ligand_mpnn": "ligandmpnn_v_32_020_25.pt"}[spec]
        cmd = [LIG_PY, LIG, "--model_type", spec,
               "--checkpoint_" + spec, os.path.join(LP, ck),
               "--pdb_path", pdb, "--out_folder", out,
               "--number_of_batches", "1", "--batch_size", str(N),
               "--temperature", "0.1", "--seed", str(seed)]
        if spec == "global_label_membrane_mpnn":
            cmd += ["--global_transmembrane_label", "0"]
        r = run(cmd)
    else:  # guide
        cmd = [PY_, GUIDED, "--pdb", pdb, "--weights", spec,
               "--n", str(N), "--out", out]
        r = run(cmd)
    if cell_done(fa):
        return "ok", r
    return "FAIL", r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=4)
    ap.add_argument("--arms", default="")  # optional comma list override
    ap.add_argument("--retries", type=int, default=2)
    args = ap.parse_args()
    arms = args.arms.split(",") if args.arms else ORDER
    cells = [(a, t) for a in arms for t in TARGETS]
    mine = [c for i, c in enumerate(cells) if i % args.nshards == args.shard]
    print(f"P16_GEN shard={args.shard}/{args.nshards} cells={len(mine)} "
          f"CUDA={os.environ.get('CUDA_VISIBLE_DEVICES')}", flush=True)
    n_ok = n_skip = 0
    fails = []
    t00 = time.time()
    for arm, target in mine:
        last = None
        for attempt in range(1, args.retries + 1):
            t0 = time.time()
            status, r = generate(arm, target)
            dt = time.time() - t0
            if status in ("ok", "skip"):
                if status == "skip":
                    n_skip += 1
                else:
                    n_ok += 1
                print(f"{arm:11s} {target:5s} {status:4s} {dt:6.1f}s", flush=True)
                last = status
                break
            last = status
            err = (r.stderr or r.stdout or "").strip().splitlines()[-3:] if r is not None else []
            print(f"{arm:11s} {target:5s} FAIL attempt{attempt} {dt:6.1f}s :: {' | '.join(err)}",
                  flush=True)
            time.sleep(1)
        if last not in ("ok", "skip"):
            fails.append((arm, target))
    print(f"P16_GEN_DONE shard={args.shard} ok={n_ok} skip={n_skip} "
          f"fail={len(fails)} {fails} elapsed={time.time()-t00:.0f}s", flush=True)


if __name__ == "__main__":
    main()
