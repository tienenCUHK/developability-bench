import glob, os, sys, csv
import numpy as np
sys.path.insert(0, os.path.expanduser("~/invfold_survey"))
from tm import tmscore
from sd_eval import read_ca

ROOT = os.path.expanduser("~/invfold_survey")
DATA = f"{ROOT}/data"

def parse_cif_ca(path, want_chain=None):
    cols, rows, in_loop = [], [], False
    for line in open(path, errors="ignore"):
        if line.startswith("_atom_site."):
            cols.append(line.strip().split(".")[1]); in_loop = True; continue
        if in_loop and line.startswith("#"):
            break
        if in_loop and cols:
            parts = line.split()
            if len(parts) == len(cols):
                rows.append(parts)
    if not rows:
        return np.zeros((0,3))
    idx = {c:i for i,c in enumerate(cols)}
    ca = []
    chains = sorted({r[idx["label_asym_id"]] for r in rows}) if "label_asym_id" in idx else []
    chain = want_chain or (chains[0] if chains else None)
    for r in rows:
        if r[idx.get("label_atom_id", idx.get("auth_atom_id"))] != "CA":
            continue
        if chain and "label_asym_id" in idx and r[idx["label_asym_id"]] != chain:
            continue
        ca.append([float(r[idx["Cartn_x"]]), float(idx and r[idx["Cartn_y"]]), float(r[idx["Cartn_z"]])])
    return np.array(ca, dtype=np.float64)

def main():
    out = f"{ROOT}/runs/px/out"
    res = []
    for d in sorted(glob.glob(f"{out}/*/seed_1/predictions")):
        name = d.split("/")[-3]
        arm, tgt, idx = name.split("__")
        cif = sorted(glob.glob(f"{d}/*sample_1.cif")) or sorted(glob.glob(f"{d}/*.cif"))
        if not cif:
            continue
        ref, _ = read_ca(f"{DATA}/{tgt}.pdb")
        if len(ref) == 0:
            continue
        ca = parse_cif_ca(cif[0])
        if len(ca) == 0:
            print("empty", name, flush=True); continue
        n = min(len(ca), len(ref))
        sc = tmscore(ca[:n], ref[:n], L=len(ref))
        res.append((arm, tgt, int(idx), sc))
    with open(f"{ROOT}/runs/px/px_scTM.tsv","w") as f:
        f.write("arm\ttarget\tidx\tscTM\n")
        for r in res:
            f.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.4f}\n")
    print("WROTE", len(res), flush=True)

if __name__ == "__main__":
    main()
