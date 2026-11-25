import os, subprocess, sys
ROOT=os.path.expanduser("~/invfold_survey"); PY_=sys.executable
TARGETS=[l.strip() for l in open(f"{ROOT}/data/targets30.txt") if l.strip()]
for t in TARGETS:
    for T in ("0.3",):
        out=f"{ROOT}/runs/bench30/g_surf10_T{T}/{t}"
        if os.path.exists(f"{out}/seqs/{t}.fa"): continue
        cmd=[PY_, f"{ROOT}/guided_run.py","--pdb",f"{ROOT}/data/{t}.pdb","--weights","surf:10",
             "--n","8","--temp",T,"--out",out]
        r=subprocess.run(cmd,capture_output=True,text=True,env={**os.environ,"CUDA_VISIBLE_DEVICES":"0"})
        print(t,T,"ok" if os.path.exists(f"{out}/seqs/{t}.fa") else "FAIL",flush=True)
print("TEMP_DONE",flush=True)
