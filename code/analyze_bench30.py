import glob, os, sys, csv, itertools, collections, statistics as st
sys.path.insert(0, os.path.expanduser("~/invfold_survey"))
from dev_bench import descriptors, ca_burial, read_fasta
from sd_eval import read_ca

ROOT = os.path.expanduser("~/invfold_survey/runs/bench30")
DATA = os.path.expanduser("~/invfold_survey/data")
sc = collections.defaultdict(list)
for f in glob.glob(os.path.join(ROOT, "scTM_*.tsv")):
    for r in csv.DictReader(open(f), delimiter="\t"):
        sc[(r["model"], r["target"])].append(float(r["scTM"]))

native = {}
burial = {}
def nat(t):
    if t not in native:
        _, seq = read_ca(os.path.join(DATA, f"{t}.pdb"))
        native[t] = seq
        burial[t] = ca_burial(os.path.join(DATA, f"{t}.pdb"))
    return native[t], burial[t]

rows = []
for fa in sorted(glob.glob(os.path.join(ROOT, "*", "*", "seqs", "*.fa"))):
    arm, tgt = fa.split("/")[-4], fa.split("/")[-3]
    natseq, bur = nat(tgt)
    seqs = [(n, s) for n, s in read_fasta(fa)[1:]]
    if not seqs:
        continue
    L = min(len(s) for _, s in seqs)
    A = [list(s[:L]) for _, s in seqs]
    dv = st.mean(sum(a != b for a, b in zip(A[i], A[j])) / L
                 for i, j in itertools.combinations(range(len(A)), 2)) if len(A) > 1 else float("nan")
    for name, s in seqs:
        rec = sum(1 for a, b in zip(s, natseq) if a == b) / max(1, min(len(s), len(natseq)))
        rows.append(dict(arm=arm, target=tgt, seq=s, recovery=round(rec, 4), diversity=round(dv, 4),
                         **descriptors(s, bur)))
out = os.path.join(ROOT, "designs_analyzed.tsv")
if rows:
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)
print("WROTE", out, len(rows))
