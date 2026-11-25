import os, subprocess, sys
ROOT = os.path.expanduser("~/invfold_survey"); PY_ = sys.executable
TARGETS = [l.strip() for l in open(f"{ROOT}/data/targets30.txt") if l.strip()]
for t in TARGETS:
    out = f"{ROOT}/runs/bench30/hyper/{t}"
    if os.path.exists(f"{out}/seqs/{t}.fa"):
        continue
    os.makedirs(out, exist_ok=True)
    cmd = [PY_, f"{ROOT}/repos/ProteinMPNN/protein_mpnn_run.py", "--pdb_path", f"{ROOT}/data/{t}.pdb",
           "--pdb_path_chains", "A", "--num_seq_per_target", "8", "--batch_size", "8",
           "--sampling_temp", "0.1", "--out_folder", out,
           "--model_name", "v48_020_epoch300_hyper",
           "--path_to_model_weights", f"{ROOT}/repos/HyperMPNN_weights", "--seed", "0"]
    r = subprocess.run(cmd, capture_output=True, text=True, env={**os.environ, "CUDA_VISIBLE_DEVICES": "0"})
    print(t, "ok" if os.path.exists(f"{out}/seqs/{t}.fa") else "FAIL", flush=True)
print("HYPER_DONE", flush=True)
