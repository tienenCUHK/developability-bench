#!/usr/bin/env python3
"""lam sweep generation: surface-objective guidance at several lambdas.

Same sampler protocol as guided_run.py (ProteinMPNN v_48_020, T=0.1, seed 0,
16 sequences, bias_by_res_jsonl), with two differences:
  * `--mode surfglob` applies the apolar penalty at ALL positions (no exposure
    mask) -- the global-versus-positional ablation;
  * the bias is aligned to ProteinMPNN's own parsed chain length, which fixes
    the two backbones (1TEN, 1QYS) where guided_run.py's CA-count bias had the
    wrong length: 1QYS has an MSE (HETATM) residue and 1TEN has one residue
    without a CA atom.  For every other backbone this reproduces the p16
    g_surf10 bias exactly (checked against the saved bias.json).

Output: runs/lam_sweep/{arm_name}/{target}/seqs/{target}.fa
Usage:  CUDA_VISIBLE_DEVICES=0 ~/anaconda3/bin/python lam_gen.py \
            --mode surf --lam 2 --arm-name surf2
"""
import argparse
import json
import os
import subprocess
import sys
import time

ROOT = os.path.expanduser("~/invfold_survey")
PY_ = sys.executable
MPNN = f"{ROOT}/repos/ProteinMPNN/protein_mpnn_run.py"
RUNS = f"{ROOT}/runs/lam_sweep"
N = 16
TEMP = 0.1
SEED = 0
TARGETS = [ln.strip() for ln in open(f"{ROOT}/data/targets30.txt") if ln.strip()]

AA21 = "ACDEFGHIKLMNPQRSTVWYX"
APOLAR = set("GPAVILMFYW")
AA3_TO_1 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
            "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
            "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
            "MSE": "M"}


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


def chain_a_positions(pdb):
    """Replicate ProteinMPNN's parse_PDB_biounits for chain A.

    Returns a list aligned to ProteinMPNN's own residue order of
    (aa_1letter, CA_xyz_or_None).  MSE HETATM records are treated as MET and
    residues absent from the PDB are returned as ('X', None) placeholders,
    exactly like parse_PDB_biounits.
    """
    xyz, seqd = {}, {}
    min_resn, max_resn = 10 ** 6, -10 ** 6
    for raw in open(pdb, "rb"):
        line = raw.decode("utf-8", "ignore").rstrip()
        if line[:6] == "HETATM" and line[17:20] == "MSE":
            line = line.replace("HETATM", "ATOM  ").replace("MSE", "MET")
        if line[:4] != "ATOM":
            continue
        if line[21:22] != "A":
            continue
        atom = line[12:16].strip()
        resi = line[17:20]
        resn_s = line[22:27].strip()
        try:
            x, y, z = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        if resn_s[-1].isalpha():
            resa, resn = resn_s[-1], int(resn_s[:-1]) - 1
        else:
            resa, resn = "", int(resn_s) - 1
        min_resn = min(min_resn, resn)
        max_resn = max(max_resn, resn)
        xyz.setdefault(resn, {}).setdefault(resa, {}).setdefault(atom, (x, y, z))
        seqd.setdefault(resn, {}).setdefault(resa, resi)
    if max_resn < min_resn:
        return []
    pos = []
    for resn in range(min_resn, max_resn + 1):
        if resn in seqd:
            for k in sorted(seqd[resn]):
                aa = AA3_TO_1.get(seqd[resn][k], "X")
                ca = xyz.get(resn, {}).get(k, {}).get("CA")
                pos.append((aa, ca))
        else:
            pos.append(("X", None))
    return pos


def burial_of(coords, cutoff=10.0):
    n = len(coords)
    out = []
    for i in range(n):
        c = 0
        for j in range(n):
            if i == j:
                continue
            dx = coords[i][0] - coords[j][0]
            dy = coords[i][1] - coords[j][1]
            dz = coords[i][2] - coords[j][2]
            if dx * dx + dy * dy + dz * dz <= cutoff * cutoff:
                c += 1
        out.append(c)
    return out


def build_bias(pdb, lam, temp, mode):
    """Per-position 21-AA bias, aligned to ProteinMPNN's parsed chain A.

    surf     : -lam * [apolar] on the exposed half only (rank order of burial),
    surfglob : -lam * [apolar] on every position.
    Values are pre-multiplied by temp because the sampler divides by temp.
    Positions whose CA is unknown get a zero row (no guidance).
    """
    pos = chain_a_positions(pdb)
    idxs = [i for i, (_, ca) in enumerate(pos) if ca is not None]
    coords = [pos[i][1] for i in idxs]
    bur = burial_of(coords)
    n_exp = int(0.5 * len(bur))
    order = sorted(range(len(bur)), key=lambda j: bur[j])
    exposed_pos = {idxs[j] for j in order[:n_exp]}
    bias = []
    for i in range(len(pos)):
        apply_here = (mode == "surfglob") or (i in exposed_pos)
        row = [(-lam * (1.0 if a in APOLAR else 0.0) * temp) if apply_here else 0.0
               for a in AA21]
        bias.append(row)
    return bias, pos


def cell_done(fa):
    if not os.path.exists(fa):
        return False
    try:
        return len(read_fasta(fa)) >= N + 1
    except Exception:
        return False


def run_cell(arm, mode, lam, target):
    out = f"{RUNS}/{arm}/{target}"
    fa = f"{out}/seqs/{target}.fa"
    if cell_done(fa):
        return "skip", None
    os.makedirs(out, exist_ok=True)
    pdb = f"{ROOT}/data/{target}.pdb"
    bias, _pos = build_bias(pdb, lam, TEMP, mode)
    jf = os.path.join(out, "bias.json")
    json.dump({target: {"A": bias}}, open(jf, "w"))
    cmd = [PY_, MPNN, "--pdb_path", pdb, "--pdb_path_chains", "A",
           "--num_seq_per_target", str(N), "--batch_size", str(N),
           "--sampling_temp", str(TEMP), "--out_folder", out,
           "--model_name", "v_48_020", "--seed", str(SEED),
           "--bias_by_res_jsonl", jf]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if cell_done(fa):
        return "ok", r
    return "FAIL", r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["surf", "surfglob"], default="surf")
    ap.add_argument("--lam", type=float, required=True)
    ap.add_argument("--arm-name", required=True)
    ap.add_argument("--targets", default="")
    ap.add_argument("--retries", type=int, default=2)
    args = ap.parse_args()
    targets = args.targets.split(",") if args.targets else TARGETS
    print(f"LAM_GEN arm={args.arm_name} mode={args.mode} lam={args.lam} "
          f"cells={len(targets)} CUDA={os.environ.get('CUDA_VISIBLE_DEVICES')}",
          flush=True)
    ok = skip = 0
    fails = []
    t00 = time.time()
    for t in targets:
        last = None
        for attempt in range(1, args.retries + 1):
            t0 = time.time()
            status, r = run_cell(args.arm_name, args.mode, args.lam, t)
            dt = time.time() - t0
            if status in ("ok", "skip"):
                ok += status == "ok"
                skip += status == "skip"
                print(f"{args.arm_name:11s} {t:5s} {status:4s} {dt:6.1f}s", flush=True)
                last = status
                break
            last = status
            err = (r.stderr or r.stdout or "").strip().splitlines()[-3:] if r is not None else []
            print(f"{args.arm_name:11s} {t:5s} FAIL attempt{attempt} {dt:6.1f}s :: "
                  f"{' | '.join(err)}", flush=True)
        if last not in ("ok", "skip"):
            fails.append(t)
    print(f"LAM_GEN_DONE arm={args.arm_name} ok={ok} skip={skip} fail={len(fails)} "
          f"{fails} elapsed={time.time()-t00:.0f}s", flush=True)


if __name__ == "__main__":
    main()
