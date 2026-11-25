# DevMPNN-Bench · benchmark card

One page, for reviewers and for the artifact README. Everything here is backed by
`results/AGGREGATE.md` and `EXPERIMENTS.md`.

## What is being measured

Whether an inverse-folding model can be **steered** along developability axes, and what that steering costs.
Three quantities are reported together, which is the point of the benchmark:

| quantity | symbol | definition | why it is needed |
|---|---|---|---|
| foldability | scTM | TM-score of the ESMFold-refolded design against the input backbone | the field's default quality gate |
| controllability | k | OLS slope of achieved property on requested control strength | separates "accurate" from "steerable" |
| leakage | Leak | sum of unrequested property drift, normalised by requested drift | catches edits that quietly break something else |

Secondary: native recovery, sequence diversity, paired ΔscTM with a bootstrap CI, and Holm-adjusted p-values
because ten arms share one baseline.

## Tasks and data

| | |
|---|---|
| backbones | 35 monomeric PDB chains, chain A, 56–415 aa (median 91) |
| sampling | ProteinMPNN-class autoregressive decoding, T = 0.1 unless stated, fixed seed |
| designs | 16 per (arm, backbone), 5,936 designs, 5,936 folds |
| properties | surface apolar fraction, aggregation-prone 5-mers, net charge at pH 7.4, pI |
| oracles | ESMFold v1 (primary, all folds), Protenix-2 single-sequence (robustness subset) |
| statistics | unit of analysis is the backbone; paired bootstrap B = 10,000; Holm correction across arms |

## Arms

| family | arms | mechanism |
|---|---|---|
| baseline | ProteinMPNN v_48_020 (+ seed-1 replicate for repeatability only) | unchanged prior |
| retrained | SolubleMPNN, HyperMPNN | prior moved by data filtering |
| conditioned | MembraneMPNN (global label, fed "soluble") | prior moved by an input label |
| context | LigandMPNN | extra atom encoder |
| guided | surf λ10, agg λ10, charge λ+3, surf10/charge−2, surf10/charge+2 | frozen weights, closed-form per-residue logit bias |

## Headline numbers

- Repeatability floor (same weights, different seed): ΔscTM = −0.0010, 95% CI [−0.0086, +0.0076].
- Paired cost of control: HyperMPNN is the only negative point estimate (−0.0635, CI [−0.1259, −0.0117]) and the
  only CI excluding zero before correction; no arm survives Holm correction. Its mean is bimodal: median backbone
  −0.005, but three backbones collapse by 0.55–0.63.
- Guided arms move surface hydrophobicity by 17% and aggregation by 46% at ΔscTM ≈ +0.007.
- Two properties are independently controllable: a 3×3 (surface × charge) grid is near-orthogonal.
- Unconstrained scalar guidance leaks: optimising GRAVY drives net charge from +0.2 to +20 while scTM stays 0.93.
- Oracle change (ESMFold → Protenix single-sequence) shifts absolute scTM by 0.15 and inverts the arm ranking
  (Spearman −0.18), so only within-oracle paired comparisons are reported.

## What this benchmark does NOT cover

Membrane proteins, complexes, and de novo backbones (single scenario). Non-linear developability objectives
(stability, real aggregation predictors). Wet-lab validation. A second sampling temperature for every arm
(only the surface arm was ablated).

## Reproduce

```
python prep_targets.py           # download + validate the 35 backbones
python bench30_grid.py           # sample all published + guided arms
python dev_fold.py <arms> <tag>  # ESMFold fold + scTM (resumable)
python analyze_bench30.py        # properties, recovery, diversity
python dedup_aggregate.py        # deduplicate, aggregate, paired bootstrap
```

Compute: one A6000 is enough; the whole v1 sweep is roughly 12 GPU-hours of folding plus a few CPU-hours of
sampling. No training is involved for any guided arm.
