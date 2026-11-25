"""Shrake-Rupley SASA on the backbone, to validate the CA-burial exposure proxy.

For each backbone: compute per-residue SASA from N,CA,C,O atoms, define exposed residues
as the least-buried half by SASA, and compare with the CA-neighbour-count definition used
throughout the benchmark.
"""
import os, sys, glob, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.expanduser("~/invfold_survey/data")
RADII = dict(N=1.55, CA=1.70, C=1.70, O=1.52)
PROBE = 1.4

def read_atoms(pdb, chain="A"):
    out = {}
    for line in open(pdb, errors="ignore"):
        if line.startswith("ATOM") and line[21] == chain:
            name = line[12:16].strip()
            if name in RADII:
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                out.setdefault(line[22:27], {})[name] = xyz
    return [(k, v) for k, v in out.items() if "CA" in v]

def sasa_per_residue(residues, n_sphere=92):
    # Fibonacci sphere
    i = np.arange(n_sphere) + 0.5
    phi = np.arccos(1 - 2*i/n_sphere); theta = np.pi*(1+5**0.5)*i
    sphere = np.stack([np.cos(theta)*np.sin(phi), np.sin(theta)*np.sin(phi), np.cos(phi)], 1)
    all_atoms, owner = [], []
    for ri, (_, atoms) in enumerate(residues):
        for nm, xyz in atoms.items():
            all_atoms.append(xyz); owner.append((ri, RADII[nm]+PROBE))
    all_atoms = np.array(all_atoms); owner = np.array(owner)
    sasa = np.zeros(len(residues))
    for ri, (_, atoms) in enumerate(residues):
        for nm, xyz in atoms.items():
            r = RADII[nm]+PROBE
            pts = xyz + r*sphere
            d = np.linalg.norm(pts[:,None,:]-all_atoms[None,:,:], axis=2)
            buried = (d < owner[None,:,1]).any(1)
            sasa[ri] += (r*r*4*np.pi/n_sphere)*(~buried).sum()
    return sasa

def ca_burial(pdb, chain="A", cutoff=10.0):
    xs=[]
    for line in open(pdb, errors="ignore"):
        if line.startswith("ATOM") and line[21]==chain and line[12:16].strip()=="CA":
            xs.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    xs=np.array(xs); n=len(xs)
    d=np.linalg.norm(xs[:,None,:]-xs[None,:,:],axis=2)
    return (d<=cutoff).sum(1)-1

if __name__ == "__main__":
    ids = ["1UBQ","1PGB","1FKB","3CHY","1AKE","2LZM"]
    print(f"{'pdb':6s} {'L':>4s} {'spearman':>9s} {'overlap_exposed':>15s}")
    from scipy.stats import spearmanr
    for t in ids:
        p = f"{DATA}/{t}.pdb"
        if not os.path.exists(p): continue
        res = read_atoms(p)
        s = sasa_per_residue(res)
        b = ca_burial(p)
        n = min(len(s), len(b))
        rho = spearmanr(s[:n], b[:n])[0]
        exp_sasa = set(np.argsort(s[:n])[: n//2])
        exp_bur = set(np.argsort(b[:n])[: n//2])
        ov = len(exp_sasa & exp_bur)/ (n//2)
        print(f"{t:6s} {n:4d} {rho:9.3f} {ov:15.2f}", flush=True)
