import os, subprocess, sys
PY_ = sys.executable
LIG_PY = os.path.expanduser("~/anaconda3/envs/antifold/bin/python")
ROOT = os.path.expanduser("~/invfold_survey")
MPNN = f"{ROOT}/repos/ProteinMPNN/protein_mpnn_run.py"
LIG = f"{ROOT}/repos/LigandMPNN/run.py"
LP = f"{ROOT}/repos/LigandMPNN/model_params"
TARGETS = [l.strip() for l in open(f"{ROOT}/data/targets30.txt") if l.strip()]
N = 8
ARMS = {
    "vanilla":   ("pmpnn", "vanilla_model_weights"),
    "soluble":   ("pmpnn", "soluble_model_weights"),
    "mem_global":("lig", "global_label_membrane_mpnn"),
    "ligand":    ("lig", "ligand_mpnn"),
    "g_surf10":  ("guide", "surf:10"),
    "g_charge3": ("guide", "charge:3"),
    "g_s10_cm2": ("guide", "surf:10,charge:-2"),
    "g_s10_cp2": ("guide", "surf:10,charge:2"),
}
def run(cmd, env=None):
    return subprocess.run(cmd, capture_output=True, text=True,
                          env={**os.environ, "CUDA_VISIBLE_DEVICES": "0", **(env or {})})
for arm, (kind, spec) in ARMS.items():
    for t in TARGETS:
        out = f"{ROOT}/runs/bench30/{arm}/{t}"
        if os.path.exists(f"{out}/seqs/{t}.fa"):
            continue
        os.makedirs(out, exist_ok=True)
        pdb = f"{ROOT}/data/{t}.pdb"
        if kind == "pmpnn":
            cmd = [PY_, MPNN, "--pdb_path", pdb, "--pdb_path_chains", "A",
                   "--num_seq_per_target", str(N), "--batch_size", str(N),
                   "--sampling_temp", "0.1", "--out_folder", out, "--model_name", "v_48_020",
                   "--seed", "0", "--path_to_model_weights", f"{ROOT}/repos/ProteinMPNN/{spec}"]
            r = run(cmd)
            ok = os.path.exists(f"{out}/seqs/{t}.fa")
        elif kind == "lig":
            ck = {"global_label_membrane_mpnn": "global_label_membrane_mpnn_v_48_020.pt",
                  "ligand_mpnn": "ligandmpnn_v_32_020_25.pt"}[spec]
            cmd = [LIG_PY, LIG, "--model_type", spec, "--checkpoint_"+spec, os.path.join(LP, ck),
                   "--pdb_path", pdb, "--out_folder", out, "--number_of_batches", "1",
                   "--batch_size", str(N), "--temperature", "0.1", "--seed", "0"]
            extra = ["--global_transmembrane_label", "0"] if spec == "global_label_membrane_mpnn" else []
            r = run(cmd + extra)
            ok = os.path.exists(f"{out}/seqs/{t}.fa")
        else:
            cmd = [PY_, f"{ROOT}/guided_run.py", "--pdb", pdb, "--weights", spec,
                   "--n", str(N), "--out", out]
            r = run(cmd)
            ok = os.path.exists(f"{out}/seqs/{t}.fa")
        print(f"{arm:11s} {t:5s} {'ok' if ok else 'FAIL'}", flush=True)
print("BENCH30_DONE", flush=True)
