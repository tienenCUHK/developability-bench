import os, sys, numpy as np
sys.path.insert(0, os.path.expanduser("~/invfold_survey"))
from sasa_validate import read_atoms, sasa_per_residue, ca_burial
from scipy.stats import spearmanr
DATA = os.path.expanduser("~/invfold_survey/data")
print("pdb     L   rho(SASA,burial)  overlap_exposed")
for t in ["1UBQ","1PGB","1FKB","3CHY","1AKE","2LZM"]:
    p = f"{DATA}/{t}.pdb"
    s = sasa_per_residue(read_atoms(p)); b = ca_burial(p)
    n = min(len(s), len(b)); rho = spearmanr(s[:n], b[:n])[0]
    exp_s = set(np.argsort(s[:n])[n//2:])   # highest SASA = exposed
    exp_b = set(np.argsort(b[:n])[:n//2])   # lowest CA-neighbour count = exposed
    print(f"{t:6s} {n:4d} {rho:16.3f} {len(exp_s & exp_b)/(n//2):15.2f}")
