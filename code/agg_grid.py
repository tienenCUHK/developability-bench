import os, subprocess, sys
ROOT = os.path.expanduser("~/invfold_survey")
PY_ = sys.executable
TARGETS = [l.strip() for l in open(f"{ROOT}/data/targets30.txt") if l.strip()]
for t in TARGETS:
    out = f"{ROOT}/runs/bench30/g_agg10/{t}"
    if os.path.exists(f"{out}/seqs/{t}.fa"):
        continue
    cmd = [PY_, f"{ROOT}/guided_run.py", "--pdb", f"{ROOT}/data/{t}.pdb",
           "--weights", "agg:10", "--n", "8", "--out", out]
    r = subprocess.run(cmd, capture_output=True, text=True, env={**os.environ, "CUDA_VISIBLE_DEVICES": "0"})
    print(t, "ok" if os.path.exists(f"{out}/seqs/{t}.fa") else "FAIL", flush=True)
print("AGG_DONE", flush=True)
