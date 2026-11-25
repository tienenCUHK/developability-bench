# Aggregated DevMPNN-Bench results — power-up run (16 designs per cell)

- unique (arm,target,design) folds: 5936; duplicate rows dropped: 0
- designs analyzed: 5936
- fold result files: 4
- backbones with any fold: 35
- protocol: 16 designs per (arm, backbone), T=0.1, seed 0 (vanilla_s1: seed 1), same 35 backbones as the 8-design run
- bootstrap: B=10000, unit = backbone mean, paired vs vanilla, seed=0; two-sided bootstrap p with add-one correction; Holm across the 9 comparison arms (vanilla_s1 is the seed-1 repeatability control, not an arm, and is excluded from the family)

| arm | n des | nT folded | scTM | surf | agg5 | charge | pI | recovery | diversity |
|---|---|---|---|---|---|---|---|---|---|
| vanilla | 560 | 35 | 0.8444 | 0.430 | 0.212 | +1.95 | 7.87 | 0.436 | 0.283 |
| vanilla_s1 | 560 | 35 | 0.8434 | 0.434 | 0.212 | +1.89 | 7.82 | 0.437 | 0.288 |
| soluble | 560 | 35 | 0.8457 | 0.378 | 0.168 | -1.57 | 6.64 | 0.425 | 0.283 |
| hyper | 560 | 35 | 0.7809 | 0.411 | 0.213 | +4.30 | 8.40 | 0.376 | 0.258 |
| mem_global | 528 | 33 | 0.8748 | 0.376 | 0.178 | +0.30 | 7.49 | 0.439 | 0.288 |
| ligand | 528 | 33 | 0.8717 | 0.397 | 0.176 | -1.00 | 6.64 | 0.443 | 0.263 |
| g_surf10 | 528 | 33 | 0.8703 | 0.259 | 0.120 | +1.72 | 7.87 | 0.446 | 0.272 |
| g_agg10 | 528 | 33 | 0.8718 | 0.291 | 0.124 | +1.96 | 7.97 | 0.448 | 0.280 |
| g_charge3 | 528 | 33 | 0.8618 | 0.413 | 0.206 | -8.52 | 5.35 | 0.456 | 0.277 |
| g_s10_cm2 | 528 | 33 | 0.8726 | 0.250 | 0.112 | +12.83 | 9.94 | 0.444 | 0.276 |
| g_s10_cp2 | 528 | 33 | 0.8671 | 0.258 | 0.120 | -8.70 | 5.57 | 0.448 | 0.266 |

## Paired bootstrap vs ProteinMPNN, scTM (B=10,000, over backbones)

| arm | d scTM | 95% CI | nT | p (raw) | p (Holm) |
|---|---|---|---|---|---|
| soluble | +0.0013 | [-0.0118, +0.0169] | 35 | 0.8918 | 1.0000 |
| hyper | -0.0635 | [-0.1259, -0.0117] | 35 | 0.0063 | 0.0567 |
| mem_global | +0.0207 | [-0.0013, +0.0561] | 33 | 0.0955 | 0.7639 |
| ligand | +0.0176 | [-0.0033, +0.0513] | 33 | 0.2077 | 1.0000 |
| g_surf10 | +0.0066 | [-0.0054, +0.0215] | 33 | 0.3397 | 1.0000 |
| g_agg10 | +0.0081 | [-0.0020, +0.0216] | 33 | 0.1609 | 1.0000 |
| g_charge3 | -0.0019 | [-0.0100, +0.0069] | 33 | 0.6306 | 1.0000 |
| g_s10_cm2 | +0.0089 | [-0.0057, +0.0278] | 33 | 0.3035 | 1.0000 |
| g_s10_cp2 | +0.0034 | [-0.0116, +0.0219] | 33 | 0.7358 | 1.0000 |

## Repeatability: seed-to-seed floor (vanilla_s1 vs vanilla)

| control | d scTM | 95% CI | nT | p (raw) |
|---|---|---|---|---|
| vanilla_s1 (seed 1 vs seed 0) | -0.0010 | [-0.0086, +0.0076] | 35 | 0.7716 |

This is the benchmark's resolution floor; it is not part of the Holm family.

## Paired bootstrap vs ProteinMPNN, targeted properties (per-backbone means, absolute change)

| arm | d surface apolar [95% CI] | d aggregation 5-mer [95% CI] | d net charge [95% CI] | d pI [95% CI] |
|---|---|---|---|---|
| vanilla_s1 | +0.004 [-0.004, +0.013] | -0.000 [-0.009, +0.008] | -0.06 [-0.55, +0.43] | -0.05 [-0.21, +0.11] |
| soluble | -0.052 [-0.078, -0.028] | -0.045 [-0.070, -0.022] | -3.51 [-5.02, -2.10] | -1.23 [-1.71, -0.77] |
| hyper | -0.019 [-0.044, +0.007] | +0.000 [-0.022, +0.022] | +2.35 [-0.49, +5.59] | +0.53 [-0.09, +1.17] |
| mem_global | -0.048 [-0.080, -0.019] | -0.032 [-0.063, -0.006] | -1.28 [-3.08, +0.26] | -0.29 [-0.70, +0.12] |
| ligand | -0.026 [-0.053, +0.000] | -0.034 [-0.062, -0.009] | -2.58 [-3.90, -1.05] | -1.13 [-1.51, -0.76] |
| g_surf10 | -0.171 [-0.197, -0.147] | -0.103 [-0.135, -0.077] | -0.43 [-1.40, +0.40] | -0.05 [-0.33, +0.23] |
| g_agg10 | -0.139 [-0.161, -0.119] | -0.100 [-0.126, -0.076] | -0.19 [-1.18, +0.66] | +0.05 [-0.19, +0.31] |
| g_charge3 | -0.017 [-0.032, -0.004] | -0.017 [-0.033, -0.004] | -10.67 [-13.15, -8.60] | -2.57 [-2.98, -2.14] |
| g_s10_cm2 | -0.180 [-0.208, -0.155] | -0.112 [-0.143, -0.085] | +10.69 [+7.74, +14.86] | +2.02 [+1.54, +2.53] |
| g_s10_cp2 | -0.172 [-0.199, -0.149] | -0.104 [-0.134, -0.078] | -10.84 [-15.50, -7.81] | -2.35 [-2.75, -1.92] |

Property comparisons use all designs of the arm; a property is skipped for a backbone when it is undefined for every design there (multi-chain LigandMPNN designs can have no chain-A burial map).

## Per-backbone heterogeneity of d scTM (paired vs vanilla)

| arm | nT | mean d | median d | median abs(d) | max abs(d) | max-abs backbone | #(d < -0.2) |
|---|---|---|---|---|---|---|---|
| vanilla_s1 | 35 | -0.0010 | -0.0007 | 0.0025 | 0.0949 | 2IGD | 0 |
| soluble | 35 | +0.0013 | -0.0012 | 0.0073 | 0.1491 | 2IGD | 0 |
| hyper | 35 | -0.0635 | -0.0045 | 0.0079 | 0.6318 | 1TIT | 3 |
| mem_global | 33 | +0.0207 | -0.0012 | 0.0030 | 0.4936 | 1TEN | 0 |
| ligand | 33 | +0.0176 | +0.0019 | 0.0066 | 0.4974 | 1TEN | 0 |
| g_surf10 | 33 | +0.0066 | +0.0006 | 0.0060 | 0.1790 | 2IGD | 0 |
| g_agg10 | 33 | +0.0081 | +0.0010 | 0.0054 | 0.1601 | 2IGD | 0 |
| g_charge3 | 33 | -0.0019 | -0.0006 | 0.0044 | 0.0901 | 2IGD | 0 |
| g_s10_cm2 | 33 | +0.0089 | -0.0021 | 0.0045 | 0.2275 | 2IGD | 0 |
| g_s10_cp2 | 33 | +0.0034 | -0.0019 | 0.0053 | 0.2166 | 2IGD | 0 |

- `hyper` collapses (d < -0.2): 1TIT -0.632, 3PGK -0.592, 1SRL -0.553

HyperMPNN mean-vs-median: mean d = -0.0635, median d = -0.0045, mean excluding the 3 collapse(s) = -0.0140 (n=32).

| backbone | vanilla scTM | hyper scTM | d |
|---|---|---|---|
| 1SRL | 0.7729 | 0.2198 | -0.5531 |
| 1TIT | 0.8616 | 0.2298 | -0.6318 |
| 3PGK | 0.7555 | 0.1637 | -0.5918 |

## Sign check: 16-design vs 8-design reference run (d scTM vs vanilla)

| arm | d scTM (16 designs) | d scTM (8 designs, results/AGGREGATE.md) | sign |
|---|---|---|---|
| vanilla_s1 | -0.0010 | -0.0000 | FLIP |
| soluble | +0.0013 | +0.0031 | SAME |
| hyper | -0.0635 | -0.0628 | SAME |
| mem_global | +0.0207 | +0.0248 | SAME |
| ligand | +0.0176 | +0.0207 | SAME |
| g_surf10 | +0.0066 | +0.0095 | SAME |
| g_agg10 | +0.0081 | +0.0114 | SAME |
| g_charge3 | -0.0019 | -0.0044 | SAME |
| g_s10_cm2 | +0.0089 | +0.0072 | SAME |
| g_s10_cp2 | +0.0034 | -0.0012 | FLIP |

**Sign changes: vanilla_s1 (-0.0000 -> -0.0010); g_s10_cp2 (-0.0012 -> +0.0034).**

## Caveats and sensitivity

1. **1TEN parse artifact (affects ProteinMPNN-family arms only).** ProteinMPNN emits 90-residue designs for 1TEN while chain A of data/1TEN.pdb has 89 C-alpha atoms (LigandMPNN emits 89), so the scTM alignment is shifted for vanilla/soluble/hyper/vanilla_s1 on that backbone and they score ~0.44 there, whereas mem_global/ligand score ~0.94. The artifact is identical in the 8-design reference run (vanilla 1TEN mean 0.442 there vs 0.442 here). It inflates the paired mean of the two LigandMPNN arms by about +0.015.

   Excluding 1TEN (sensitivity only, not a protocol change):

   | arm | d scTM | 95% CI | nT |
   |---|---|---|---|
   | vanilla_s1 | -0.0010 | [-0.0090, +0.0075] | 34 |
   | soluble | +0.0012 | [-0.0124, +0.0168] | 34 |
   | hyper | -0.0654 | [-0.1302, -0.0132] | 34 |
   | mem_global | +0.0059 | [-0.0029, +0.0181] | 32 |
   | ligand | +0.0026 | [-0.0050, +0.0117] | 32 |
   | g_surf10 | +0.0066 | [-0.0053, +0.0215] | 33 |
   | g_agg10 | +0.0081 | [-0.0021, +0.0214] | 33 |
   | g_charge3 | -0.0019 | [-0.0096, +0.0066] | 33 |
   | g_s10_cm2 | +0.0089 | [-0.0056, +0.0282] | 33 |
   | g_s10_cp2 | +0.0034 | [-0.0114, +0.0223] | 33 |

2. **Significance summary.** Only HyperMPNN's 95% CI excludes zero before correction (raw p = 0.0063); no arm survives Holm at alpha = 0.05 (HyperMPNN Holm p = 0.0567). In the 8-design run MembraneMPNN also excluded zero; at 16 designs its CI now contains zero (raw p = 0.0955, and its mean is dominated by the 1TEN artifact above).
3. **Seed-to-seed floor** is measured by vanilla_s1 vs vanilla at 16 designs (see the repeatability table); it is wider than the 8-design floor, which matters when reading any comparison whose CI is of that order.

## Coverage notes

- 5936 designs from 371/385 (arm,backbone) cells; the 14 missing cells are exactly the same ones missing in the 8-design reference run, so arm x backbone coverage is identical between the two runs (35/35 backbones for vanilla, soluble, hyper, vanilla_s1; 33/35 for mem_global, ligand and the five guided arms).
- `vanilla`: 35/35 backbones folded
- `vanilla_s1`: 35/35 backbones folded
- `soluble`: 35/35 backbones folded
- `hyper`: 35/35 backbones folded
- `mem_global`: 33/35 backbones folded
- `ligand`: 33/35 backbones folded
- `g_surf10`: 33/35 backbones folded
- `g_agg10`: 33/35 backbones folded
- `g_charge3`: 33/35 backbones folded
- `g_s10_cm2`: 33/35 backbones folded
- `g_s10_cp2`: 33/35 backbones folded

- missing for `mem_global`: 1AYE, 1BNI
- missing for `ligand`: 1AYE, 1BNI
- missing for `g_surf10`: 1QYS, 1TEN
- missing for `g_agg10`: 1QYS, 1TEN
- missing for `g_charge3`: 1QYS, 1TEN
- missing for `g_s10_cm2`: 1QYS, 1TEN
- missing for `g_s10_cp2`: 1QYS, 1TEN

Generation failures (identical cells to the 8-design reference run): `mem_global`/`ligand` fail on 1BNI and 1AYE with a prody shape mismatch inside LigandMPNN; the five `g_*` arms fail on 1TEN and 1QYS because the guided bias length is 89/91 while ProteinMPNN pads to 90/92. All are pre-existing tool failures, not regressions of this run.

