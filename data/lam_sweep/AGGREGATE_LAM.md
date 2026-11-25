# Lambda sweep + global-vs-positional ablation — surface objective

- protocol: 16 designs per (lambda, backbone), ProteinMPNN v_48_020, T=0.1, seed 0, 35 backbones; T=0.1; scTM = ESMFold refold vs chain-A backbone (same folds/protocol as DevMPNN-Bench)
- lambda=10 reuses the 16-design g_surf10 cells from the main p16 run (33 backbones) plus newly sampled 1TEN/1QYS, so all lambda values cover 35 backbones
- bootstrap: B=10000, unit = backbone, seed=0; 95% percentile CI; two-sided bootstrap p with add-one correction

## Deliverable 1 — surface arm, lambda sweep (35 backbones)

| lambda | n des | nT folded | scTM | surf apolar | agg5 | charge | pI | diversity | recovery | d scTM vs vanilla [95% CI] |
|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 560 | 35 | 0.8472 | 0.368 | 0.171 | +2.08 | 7.84 | 0.282 | 0.459 | +0.0028 [-0.0006, +0.0067] |
| 5 | 560 | 35 | 0.8472 | 0.308 | 0.136 | +1.60 | 7.80 | 0.277 | 0.457 | +0.0028 [-0.0054, +0.0123] |
| 10 | 560 | 35 | 0.8507 | 0.253 | 0.114 | +1.52 | 7.79 | 0.270 | 0.450 | +0.0062 [-0.0051, +0.0207] |
| 20 | 560 | 35 | 0.8471 | 0.178 | 0.094 | +1.74 | 7.96 | 0.283 | 0.438 | +0.0027 [-0.0133, +0.0218] |

## Per-backbone controllability slope k (achieved property vs lambda)

Per backbone: OLS slope over lambda in [2.0, 5.0, 10.0, 20.0]; estimate = mean over n=35 backbones with bootstrap CI.

| quantity | k (per lambda unit) | 95% CI | median per-backbone k | mean R^2 | n backbones | p (k=0) |
|---|---|---|---|---|---|---|
| surface apolar | -0.01005 | [-0.01113, -0.00905] | -0.00994 | 0.940 | 35 | 0.0001 |
| scTM | +0.00001 | [-0.00072, +0.00083] | -0.00022 | 0.515 | 35 | 0.9804 |

Direction: surface apolar changes -0.01005 per unit lambda (negative = the guidance moves the achieved property in the intended direction); scTM changes +0.00001 per unit lambda.

## Deliverable 2 — exposure mask vs global apolar penalty

Both arms at lambda = 20; paired over backbones.

| arm | n des | nT folded | scTM | surf apolar | agg5 | charge | pI | diversity |
|---|---|---|---|---|---|---|---|---|
| mask (surf:20) | 560 | 35 | 0.8471 | 0.178 | 0.094 | +1.74 | 7.96 | 0.283 |
| global (surfglob:20) | 560 | 35 | 0.8146 | 0.159 | 0.028 | +3.28 | 8.29 | 0.295 |

- paired **d scTM (global - mask)** = -0.0325 [-0.0748, -0.0019], nT=35, raw p=0.0279
- paired **d surface apolar (global - mask)** = -0.0186 [-0.0248, -0.0130], nT=35, raw p=0.0001 (negative = global controls better)

## One-line answer

Foldability does not degrade with lambda (k_scTM = +0.00001 per lambda unit, 95% CI [-0.00072, +0.00083]); at the best lambda = 20, the exposure mask does real work: dropping it (global penalty) additionally lowers surface apolar by 0.0186 and aggregation 5-mer by 0.0662 but costs -0.0325 scTM (95% CI [-0.0748, -0.0019], p = 0.0279).

## Best lambda rule

Largest lambda whose paired d scTM vs vanilla has 95% CI lower bound >= -0.02 (no meaningful foldability cost); selected lambda = **20**.
