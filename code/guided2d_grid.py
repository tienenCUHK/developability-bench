import os, subprocess, sys
PY_ = sys.executable
ROOT = os.path.expanduser("~/invfold_survey")
TARGETS = ["1UBQ", "1PGB", "1FKB", "3CHY", "2LZM", "1AKE"]
SURF = [0, 5, 10]
CHARGE = [-2, 0, 2]
N = 8
for s in SURF:
    for c in CHARGE:
        cfg = f"2d_s{s}_c{c}"
        for t in TARGETS:
            out = f"{ROOT}/runs/grid2d/{cfg}/{t}"
            if os.path.exists(f"{out}/seqs/{t}.fa"):
                continue
            cmd = [PY_, f"{ROOT}/guided_run.py", "--pdb", f"{ROOT}/data/{t}.pdb",
                   "--weights", f"surf:{s},charge:{c}", "--n", str(N), "--out", out]
            r = subprocess.run(cmd, capture_output=True, text=True,
                               env={**os.environ, "CUDA_VISIBLE_DEVICES": "0"})
            print(f"{cfg:14s} {t:5s} {'ok' if r.returncode==0 else 'FAIL'}", flush=True)
print("GRID2D_DONE", flush=True)
