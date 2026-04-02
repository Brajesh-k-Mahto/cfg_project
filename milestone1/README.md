# Milestone 1 Submission (Public Link)

This folder contains **standalone scripts** (no notebooks) and a README with exact run commands.

### Step 0: Download genome (one-time)
```bash
bash data/genome/get_hg38.sh
'''



## A) Scripts included
1. `MarkovCrossValidation.py`
   - Primary milestone script.
   - Inputs: Markov order `m`, folds `k`, single chromosome TSV, TF name.
   - Outputs: k ROC curves + k PR curves + mean ROC-AUC + mean PR-AUC.

2. `MarkovSweepChr4.py`
   - Convenience script for the form questions that require running **chr4** for **m=0..10** with `k=5`.
   - Prints 11 ROC-AUC values (3 decimals), 11 PR-AUC values (3 decimals), and 11 timing values.

3. `simplerVersion.py`
   - Takes **only 2 args**: FASTA file + order m.
   - Trains Markov model on all sequences and prints log-likelihood per sequence (one line per sequence).

## B) Requirements
* Python 3.9+
* Libraries:
  * numpy
  * pandas
  * pysam
  * scikit-learn
  * matplotlib

Install:
```bash
pip install numpy pandas pysam scikit-learn matplotlib
```

## C) Expected file structure (relative paths)
```
.
├── milestone1/
│   ├── MarkovCrossValidation.py
│   ├── MarkovSweepChr4.py
│   ├── simplerVersion.py
│   └── README.md
├── data/
│   ├── train/
│   │   └── chr4_200bp_bins.tsv   (example; any labeled chr TSV works)
│   └── genome/
│       └── hg38.fa               (FASTA file; script will create hg38.fa.fai if missing)
└── results/
```

## D) Pick a TF (required by form)
Pick one TF from: `CTCF`, `REST`, `EP300`.
You will use the same TF for the chr4 sweep questions.

## E) How to run (Milestone 1 core requirement)
Example (chr4, CTCF, m=6, k=5):
```bash
python milestone1/MarkovCrossValidation.py --order 6 --k 5 --input data/train/chr4_200bp_bins.tsv --TF CTCF
```

Outputs:
* `results/milestone1/CTCF_m6_roc_fold1.png` ... `roc_fold5.png`
* `results/milestone1/CTCF_m6_pr_fold1.png`  ... `pr_fold5.png`
* `results/milestone1/summary_CTCF_m6.txt`

Upload any one ROC and PR plot for one fold as required by the form.

## F) Run chr4 for m=0..10 with k=5 (form questions)
This prints exactly what you need to paste into the form (3 decimals), and timing values:

```bash
python milestone1/MarkovSweepChr4.py --input data/train/chr4_200bp_bins.tsv --TF CTCF --k 5
```

By default it saves ROC/PR plots for m=6 only (one set of 5 folds).
To save plots for ALL m:
```bash
python milestone1/MarkovSweepChr4.py --input data/train/chr4_200bp_bins.tsv --TF CTCF --k 5 --save_plots_for_m all
```

## G) simplerVersion requirement (2 arguments only)
```bash
python milestone1/simplerVersion.py path/to/sequences.fasta 6
```
Prints one log-likelihood per line (one per sequence).
