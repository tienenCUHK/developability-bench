import itertools, os, subprocess, sys
PY_ = sys.executable
ROOT = os.path.expanduser("~/invfold_survey")
TARGETS = ["1UBQ", "1PGB", "1FKB", "3CHY", "2LZM", "1AKE"]
GRID = {"surf": [0, 2, 5, 10], "gravy": [0, 2, 4, 8], "charge": [-3, 0, 3, 6]}
n = int(os.environ.get("NSEQ", "8"))
for prop, lams in GRID.items():
    for lam in lams:
        for t in TARGETS:
            out = f"{ROOT}/runs/grid/{prop}_{lam}/{t}"
            if os.path.exists(f"{out}/seqs/{t}.fa"):
                continue
            cmd = [PY_, f"{ROOT}/guided_run.py", "--pdb", f"{ROOT}/data/{t}.pdb",
                   "--property", prop, "--lam", str(lam), "--n", str(n), "--out", out]
            r = subprocess.run(cmd, capture_output=True, text=True, env={**os.environ, "CUDA_VISIBLE_DEVICES": "0"})
            print(f"{prop:7s} lam={lam:<5} {t:5s} {'ok' if r.returncode==0 else 'FAIL'}", flush=True)
print("GRID_DONE", flush=True)
