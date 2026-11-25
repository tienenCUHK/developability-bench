"""Deduplicate fold results, aggregate per arm, and emit results/AGGREGATE.md."""
import csv, glob, collections, statistics as st, numpy as np, os

sc = collections.defaultdict(dict)
dups = 0
for f in glob.glob('results/b30_scTM_*.tsv'):
    for r in csv.DictReader(open(f), delimiter='\t'):
        k = (r['model'], r['target'], int(r['idx']))
        if k in sc: dups += 1
        sc[k] = float(r['scTM'])
arm_t = collections.defaultdict(list)
for (a, t, i), v in sc.items():
    arm_t[(a, t)].append(v)
tg = sorted({t for _, t, _ in sc})

rows = list(csv.DictReader(open('results/b30_designs.tsv'), delimiter='\t'))
by = collections.defaultdict(list)
for r in rows:
    by[r['arm']].append(r)

rng = np.random.default_rng(0)
def boot(a, b):
    xs = [st.mean(arm_t[(a, t)]) for t in tg if arm_t.get((a, t)) and arm_t.get((b, t))]
    ys = [st.mean(arm_t[(b, t)]) for t in tg if arm_t.get((a, t)) and arm_t.get((b, t))]
    if len(xs) < 5: return None
    d = np.array(ys) - np.array(xs)
    idx = rng.integers(0, len(d), (10000, len(d)))
    bs = d[idx].mean(1)
    return d.mean(), np.percentile(bs, 2.5), np.percentile(bs, 97.5), len(d)

order = ['vanilla','vanilla_s1','soluble','hyper','mem_global','ligand',
         'g_surf10','g_agg10','g_charge3','g_s10_cm2','g_s10_cp2']
lines = [f"# Aggregated DevMPNN-Bench results\n",
         f"- unique (arm,target,design) folds: {len(sc)}; duplicate rows dropped: {dups}",
         f"- backbones with any fold: {len(tg)}\n",
         "| arm | n des | nT folded | scTM | surf | agg5 | charge | pI | recovery | diversity |",
         "|---|---|---|---|---|---|---|---|---|---|"]
for a in order:
    v = by.get(a, [])
    if not v: continue
    per = [st.mean(arm_t[(a, t)]) for t in tg if arm_t.get((a, t))]
    sa = [float(x['surf_apolar']) for x in v if x['surf_apolar'] != 'NA']
    lines.append(f"| {a} | {len(v)} | {len(per)} | {st.mean(per):.4f} | {st.mean(sa):.3f} | "
                 f"{st.mean(float(x['agg5frac']) for x in v):.3f} | {st.mean(float(x['charge74']) for x in v):+.2f} | "
                 f"{st.mean(float(x['pI']) for x in v):.2f} | {st.mean(float(x['recovery']) for x in v):.3f} | "
                 f"{st.mean(float(x['diversity']) for x in v):.3f} |")
lines += ["", "## Paired bootstrap vs ProteinMPNN (B=10,000, over backbones)", "",
          "| arm | d scTM | 95% CI | nT |", "|---|---|---|---|"]
for a in order:
    if a in ('vanilla',): continue
    r = boot('vanilla', a)
    if r: lines.append(f"| {a} | {r[0]:+.4f} | [{r[1]:+.4f}, {r[2]:+.4f}] | {r[3]} |")
open('results/AGGREGATE.md','w').write("\n".join(lines) + "\n")
print("\n".join(lines))
