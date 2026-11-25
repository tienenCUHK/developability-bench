import os, subprocess, sys
ROOT = os.path.expanduser("~/invfold_survey")
PY_ = sys.executable
MPNN = f"{ROOT}/repos/ProteinMPNN/protein_mpnn_run.py"
TARGETS = [l.strip() for l in open(f"{ROOT}/data/targets30.txt") if l.strip()]
for t in TARGETS:
    out = f"{ROOT}/runs/bench30/vanilla_s1/{t}"
    if os.path.exists(f"{out}/seqs/{t}.fa"):
        continue
    os.makedirs(out, exist_ok=True)
    cmd = [PY_, MPNN, "--pdb_path", f"{ROOT}/data/{t}.pdb", "--pdb_path_chains", "A",
           "--num_seq_per_target", "8", "--batch_size", "8", "--sampling_temp", "0.1",
           "--out_folder", out, "--model_name", "v_48_020", "--seed", "1",
           "--path_to_model_weights", f"{ROOT}/repos/ProteinMPNN/vanilla_model_weights"]
    r = subprocess.run(cmd, capture_output=True, text=True, env={**os.environ, "CUDA_VISIBLE_DEVICES": "0"})
    print(t, "ok" if os.path.exists(f"{out}/seqs/{t}.fa") else "FAIL", flush=True)
print("SEED2_DONE", flush=True)
