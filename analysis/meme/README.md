# MEME-only Analysis Workspace

This directory is intentionally separated from the Markov pipeline to avoid confusion.

## Folder layout
- `analysis/meme/scripts/`: helper scripts for FASTA prep and result assembly
- `analysis/meme/outputs/meme_input/`: bound-sequence FASTA files for MEME motif discovery
- `analysis/meme/outputs/fimo_input/`: unknown-chromosome FASTA files for FIMO scanning
- `analysis/meme/outputs/fimo_results/`: downloaded FIMO outputs and processed TF score tables
- `analysis/meme/outputs/final_scores/`: final merged score files (`chr3/10/17_predictions.tsv.gz`)

## Prerequisites
Use your existing project environment and ensure packages are installed:
- `pandas`
- `pysam`

Example:
```bash
source .venv/bin/activate
pip install pandas pysam
```

## Step 1: Build TF-specific bound FASTA for MEME
Run once per TF:
```bash
python analysis/meme/scripts/prepare_meme_bound_fasta.py --hg38 data/genome/hg38.fa --tf CTCF
python analysis/meme/scripts/prepare_meme_bound_fasta.py --hg38 data/genome/hg38.fa --tf REST
python analysis/meme/scripts/prepare_meme_bound_fasta.py --hg38 data/genome/hg38.fa --tf EP300
```
Outputs:
- `analysis/meme/outputs/meme_input/CTCF_bound.fa`
- `analysis/meme/outputs/meme_input/REST_bound.fa`
- `analysis/meme/outputs/meme_input/EP300_bound.fa`

## Step 2: Discover motifs on MEME web server
Use the online MEME tool:
- https://meme-suite.org/meme/tools/meme

For each TF:
1. Upload corresponding `*_bound.fa`
2. Alphabet: DNA
3. Keep default settings initially (or tune motif width/count later)
4. Run MEME and download motif file (`meme.txt`)

Save downloaded motif files under:
- `analysis/meme/outputs/fimo_results/motifs/`

## Step 3: Build unknown FASTA for FIMO
```bash
python analysis/meme/scripts/prepare_fimo_unknown_fasta.py --hg38 data/genome/hg38.fa
```
Outputs:
- `analysis/meme/outputs/fimo_input/chr3_unknown.fa`
- `analysis/meme/outputs/fimo_input/chr10_unknown.fa`
- `analysis/meme/outputs/fimo_input/chr17_unknown.fa`

## Step 4: Scan with FIMO web server
Use the online FIMO tool:
- https://meme-suite.org/meme/tools/fimo

For each (TF, chromosome) pair:
1. Upload TF motif file (`meme.txt` from Step 2)
2. Upload unknown FASTA (`chr*_unknown.fa`)
3. Run and download `fimo.tsv`
4. Store file using a clear name, for example:
   - `analysis/meme/outputs/fimo_results/raw/CTCF_chr3_fimo.tsv`

You will have 9 raw FIMO files (3 TFs x 3 chromosomes).

## Step 5: Convert FIMO hits to one score per bin
Run for each raw FIMO file:
```bash
python analysis/meme/scripts/process_fimo_output.py \
  --unknown_tsv data/unknown/chr3_200bp_bins_unknown.tsv \
  --fimo_tsv analysis/meme/outputs/fimo_results/raw/CTCF_chr3_fimo.tsv \
  --tf CTCF \
  --metric score \
  --out_tsv analysis/meme/outputs/fimo_results/processed/chr3_CTCF_scores.tsv
```

Notes:
- `--metric score` uses FIMO hit score directly and takes max hit score per 200bp bin
- alternatives: `--metric neglog10_p` or `--metric neglog10_q`
- bins with no motif hit get score `0.0`

## Step 6: Merge three TF files into final submission file
After you have processed files for one chromosome:
```bash
python analysis/meme/scripts/merge_meme_scores.py \
  --ctcf_tsv analysis/meme/outputs/fimo_results/processed/chr3_CTCF_scores.tsv \
  --rest_tsv analysis/meme/outputs/fimo_results/processed/chr3_REST_scores.tsv \
  --ep300_tsv analysis/meme/outputs/fimo_results/processed/chr3_EP300_scores.tsv \
  --out_tsv analysis/meme/outputs/final_scores/chr3_predictions.tsv.gz
```

Repeat for `chr10` and `chr17`.

## Deliverables from this MEME workflow
- `analysis/meme/outputs/final_scores/chr3_predictions.tsv.gz`
- `analysis/meme/outputs/final_scores/chr10_predictions.tsv.gz`
- `analysis/meme/outputs/final_scores/chr17_predictions.tsv.gz`

These match the final naming pattern required in the project description.
