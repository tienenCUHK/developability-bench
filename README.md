# DevMPNN-Bench — submission artifact

Artifact for *"What Developability Benchmarks Do Not Measure: Controllability, Leakage and Bimodal Failure in
Inverse-Folding Evaluation"* (NeurIPS Evaluations & Datasets track).

The bundle contains the manuscript, the design matrix, both aggregates, the analysis code, the dataset metadata
(Croissant core + RAI) and the environment pins.

**Review mode.** This repository is hosted under the authors' account, so the submission is declared **single-blind**
on the submission form (the E&D call for papers allows this when an artifact cannot be anonymised, and this one
is a long-lived benchmark that the authors maintain). Git commits are authored as `Anonymous Authors` for the
history, but the hosting account is public.

## What is in the box

```
paper/            LaTeX source, bibliography, figures, compiled PDF, static checkers
data/
  b30_designs.tsv         per-design matrix: arm x backbone x design, with
                          sequence, native recovery, diversity, charge74, pI, GRAVY,
                          surface-apolar fraction, aggregation 5-mer fraction
  p16/                    the 16-design main run: aggregate, per-fold scTM, descriptors
  lam_sweep/              the lambda sweep (lambda in {2,5,10,20}) and the
                          exposure-mask-vs-global ablation, with aggregates
  AGGREGATE.md            the 8-design predecessor run (used only for the
                          temperature ablation in the paper)
  px_scTM.tsv             second-oracle scTM from Protenix-2 single-sequence refolding
code/             sampling, refolding, analysis, aggregation and figure scripts
DATASHEET.md      dataset documentation
BENCHMARK_CARD.md one-page summary of what is measured and what is not
croissant.json    machine-readable metadata, Croissant 1.1 core + RAI 1.0 fields
requirements.txt  runtime pins
LICENSE           MIT (code); text and figures CC BY 4.0
```

## Quick check without a GPU

The aggregate tables and the per-design matrix are plain text; every number in the paper's §4 can be recomputed
from `data/` on a laptop:

```bash
python code/p16_aggregate.py --designs data/p16/designs_analyzed.tsv \
                             --folds  data/p16/scTM_g*.tsv \
                             --out    /tmp/aggregate_check.md
```

## Full reproduction

1. `code/prep_targets.py` downloads the 35 RCSB entries listed in `data/targets30.txt` and validates chain A.
2. `code/p16_gen.py` samples 16 designs per (arm, backbone) for all 11 arms; `code/lam_gen.py` does the same for
   the lambda sweep.
3. `code/p16_fold.py` and `code/lam_fold.py` refold every design with ESMFold v1 (three recycles) and compute the
   chain-A TM-score. These steps are GPU-bound, resumable, and skip cells already present in their output TSV.
4. `code/p16_analyze.py` computes the physicochemical descriptors, native recovery and diversity;
   `code/p16_aggregate.py` produces the paired bootstrap intervals and the Holm-adjusted p-values.

Upstream models and weights (all used under their own licences, not redistributed here):

- ProteinMPNN and the SolubleMPNN weights: `github.com/dauparas/ProteinMPNN`
- LigandMPNN, including the membrane-conditioned weights: `github.com/dauparas/LigandMPNN` (client needs `prody`)
- HyperMPNN weights: `github.com/meilerlab/HyperMPNN`, file `v48_020_epoch300_hyper.pt`, loaded by ProteinMPNN's code
- ESMFold: `fair-esm` with a vendored `openfold` (the build notes are in `paper/BUILD.md`)
- Protenix-2 for the second oracle, run with `--use_msa False`

## Notes that matter when re-running

- Refold outputs append. If a job is interrupted and restarted, deduplicate with the `(arm, backbone, design)`
  key; the aggregation scripts keep the last row and report how many duplicates they dropped.
- The `surface apolar` quantity is a CA-burial proxy, not SASA. Its validation against Shrake–Rupley SASA is in
  the paper's Appendix A.
- One documented alignment artefact (backbone 1TEN: 90 designed residues against an 89-CA reference chain) is
  retained in the main table and disclosed in the paper's limitations, with an exclusion sensitivity analysis.
- Hardware used: up to four NVIDIA RTX A6000 GPUs, never more than four concurrently. No model training is
  involved; folding dominates the runtime.

## Licence

Code MIT. Manuscript text and figures CC BY 4.0. PDB structures are public-domain RCSB entries and are identified
by accession rather than redistributed.
