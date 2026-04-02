# CFG Project Submission Package

## Contents

### Required final prediction files (MEME pipeline)
These are the files matching the expected final naming convention:
- `final_predictions_meme/chr3_predictions.tsv.gz`
- `final_predictions_meme/chr10_predictions.tsv.gz`
- `final_predictions_meme/chr17_predictions.tsv.gz`

### Baseline comparison files (Markov pipeline)
These are included for internal comparison and presentation:
- `final_predictions_markov/final_chr3_predictions.tsv`
- `final_predictions_markov/final_chr10_predictions.tsv`
- `final_predictions_markov/final_chr17_predictions.tsv`

### Presentation prep
- `presentation/slide_outline.md`

## Method Summary
1. Sequence extraction from hg38 for all 200bp bins.
2. TF-specific motif discovery using MEME on bound training bins.
3. Motif scanning on unknown chromosomes via FIMO.
4. Per-bin score aggregation (max hit score) and merge into final files.

## Reproducibility notes
- MEME scripts and workflow are isolated in `analysis/meme/`.
- Final submission-style files are generated under `analysis/meme/outputs/final_scores/` before being copied here.
- Markov baseline outputs remain in `results/final_predictions/`.

## AI acknowledgment
This project used AI assistance for code organization, debugging support, and workflow scripting. All generated code and commands were reviewed and run by the team.
