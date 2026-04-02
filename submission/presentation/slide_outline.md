# Slide Deck Outline (CFG Project)

## Slide 1: Problem and Goal
Title: Predicting TF Binding Potential in K562
- Task: predict binding scores for CTCF, REST, EP300 on chr3, chr10, chr17
- Input: 200bp bins, hg38 genome, labeled training chromosomes
- Output: one numeric score per TF per bin (higher means stronger binding)

Speaker notes:
- Explain why in vivo binding is harder than motif matching alone.
- Mention ATAC availability and that sequence-first baselines were tested.

## Slide 2: Data and Constraints
- 19 labeled chromosomes for training
- 3 unknown chromosomes for final prediction
- Labels in training: B/U for each TF
- Must fetch sequence from external genome (hg38)
- No prohibited external ChIP-based labels used

Speaker notes:
- Mention that missing/repeat bins are already filtered in provided data.
- Emphasize reproducibility and fixed bin size (200bp).

## Slide 3: Method A (Markov Baseline)
- Train TF-specific positive and negative Markov models
- Score each sequence using log-likelihood ratio
- Cross-validation for model order selection
- Generate class predictions and evaluation curves

Speaker notes:
- Mention order sweep and cross-validation approach.
- State this baseline is simple and interpretable.

## Slide 4: Method B (MEME + FIMO Pipeline)
- Build TF-specific bound FASTA sets
- Discover motifs with MEME
- Scan unknown sequences with FIMO
- Aggregate hit-level motif scores to one score per 200bp bin

Speaker notes:
- Mention local MEME setup to avoid web queue bottlenecks.
- Note script fix for FIMO output format differences.

## Slide 5: Pipeline Architecture
- Input TSV + hg38 -> sequence FASTA
- MEME motifs -> FIMO raw hits
- Hit aggregation -> per-TF score files
- Merge three TF scores -> final chr*_predictions.tsv.gz

Speaker notes:
- Mention isolated workflow directory: analysis/meme.
- Mention independent processing for each TF and chromosome.

## Slide 6: Key Results Summary
Use table format:
- For each chromosome, report:
  - Markov B count per TF
  - MEME non-zero score bins per TF
  - Top-score overlap trends

Speaker notes:
- Stronger agreement seen for CTCF/EP300 on chr10 and chr17.
- REST appears weaker in motif specificity with current motifs/settings.

## Slide 7: Interpretation
- Motif-based scoring captures meaningful TF signal
- Performance is TF-dependent and chromosome-dependent
- Low-complexity motifs can inflate hit counts; filtering/tuning helps

Speaker notes:
- Discuss biological heterogeneity and sequence context limitations.
- Explain that motif presence is necessary but not always sufficient for binding.

## Slide 8: Limitations and Future Work
- Add score calibration and threshold analysis
- Combine motif score + Markov score + ATAC as features
- Try stricter motif filtering and motif quality controls
- Evaluate ROC/PR on held-out labeled chromosomes consistently

Speaker notes:
- Propose an ensemble model for final improvement.
- Mention reproducible scripts and folder layout are ready for extension.

## Slide 9: Final Deliverables
- chr3_predictions.tsv.gz
- chr10_predictions.tsv.gz
- chr17_predictions.tsv.gz
- Reproducible scripts and workflow documentation

Speaker notes:
- Point to submission folder organization.
- Confirm all outputs were generated and validated for file format.
