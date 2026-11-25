import json, os, sys
ROOT = os.path.expanduser("~/invfold_survey")
TARGETS = ["1UBQ","1PGB","1FKB","3CHY","2LZM","1AKE"]
ARMS = ["vanilla","soluble","mem_global","ligand","g_surf10","g_agg10","g_charge3","g_s10_cm2","g_s10_cp2"]
def read_fa(p):
    recs=[];n=None;b=[]
    for l in open(p):
        l=l.strip()
        if l.startswith(">"):
            if n: recs.append("".join(b))
            n,b=l[1:],[]
        elif l: b.append(l)
    if n: recs.append("".join(b))
    return recs
items=[]
for a in ARMS:
    for t in TARGETS:
        fa=f"{ROOT}/runs/bench30/{a}/{t}/seqs/{t}.fa"
        if not os.path.exists(fa): continue
        seqs=read_fa(fa)[1:3]           # 2 designs per (arm,target)
        for i,s in enumerate(seqs):
            items.append({"name": f"{a}__{t}__{i}",
                          "sequences":[{"proteinChain":{"sequence":s,"count":1}}]})
json.dump(items, open(f"{ROOT}/runs/px/in.json","w"))
print("items", len(items))
