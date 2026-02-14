#!/usr/bin/env python3
"""
MarkovCrossValidation.py

Sequence-only Markov model classifier for TF binding (B vs U) with k-fold CV.
Generates k ROC curves and k PR curves and prints mean ROC-AUC and mean PR-AUC.

Inputs:
  --order (m): Markov order (0..10)
  --k: number of folds (3..5)
  --input: path to a single chromosome training TSV (e.g., data/train/chr1_200bp_bins.tsv)
  --TF: one of {CTCF, REST, EP300}
  --genome: reference genome FASTA (default: data/genome/hg38.fa)
  --outdir: output directory for plots and summary (default: results/milestone1)

Output:
  - <outdir>/<TF>_m<order>_roc_fold<i>.png
  - <outdir>/<TF>_m<order>_pr_fold<i>.png
  - <outdir>/summary_<TF>_m<order>.txt
"""

import os
import math
import argparse
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt

from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, precision_recall_curve, auc

ALPHABET = ("A", "C", "G", "T")


def is_valid_dna(seq: str) -> bool:
    s = seq.upper()
    return all(ch in ALPHABET for ch in s)


def train_markov(seqs, m: int):
    """
    Train an order-m Markov model with Laplace smoothing.
    Returns: dict context -> dict(base->prob)
    """
    if m < 0:
        raise ValueError("order m must be >= 0")

    counts = defaultdict(lambda: defaultdict(int))

    for s in seqs:
        s = s.upper()
        if len(s) <= m or (not is_valid_dna(s)):
            continue

        for i in range(m, len(s)):
            ctx = s[i - m:i] if m > 0 else ""
            nxt = s[i]
            counts[ctx][nxt] += 1

    probs = {}
    for ctx, next_counts in counts.items():
        denom = sum(next_counts.get(b, 0) + 1 for b in ALPHABET)  # Laplace +1
        probs[ctx] = {b: (next_counts.get(b, 0) + 1) / denom for b in ALPHABET}

    return probs


def log_prob_seq(model_probs, seq: str, m: int) -> float:
    """
    log P(seq) under Markov model. Uses uniform 1/4 for unseen contexts.
    """
    s = seq.upper()
    if len(s) <= m or (not is_valid_dna(s)):
        return float("-inf")

    logp = 0.0
    for i in range(m, len(s)):
        ctx = s[i - m:i] if m > 0 else ""
        nxt = s[i]
        if ctx in model_probs:
            p = model_probs[ctx][nxt]
        else:
            p = 0.25
        logp += math.log(p)
    return logp


def llr_score(pos_model, neg_model, seq: str, m: int) -> float:
    return log_prob_seq(pos_model, seq, m) - log_prob_seq(neg_model, seq, m)


def load_sequences_and_labels(tsv_path: str, fasta_path: str, tf_name: str):
    """
    Loads sequences and labels (B/U) from one chromosome TSV.
    Returns:
      seqs: np.ndarray of strings
      y: np.ndarray of ints (1 for B, 0 for U)
    """
    df = pd.read_csv(tsv_path, sep="\t")
    fa = pysam.FastaFile(fasta_path)

    seqs = []
    ys = []

    for _, r in df.iterrows():
        label = r[tf_name]
        if label not in ("B", "U"):
            continue

        chrom = r["chr"]
        start = int(r["start"])
        end = int(r["end"])

        seq = fa.fetch(chrom, start, end).upper()
        if len(seq) != (end - start):
            continue
        if not is_valid_dna(seq):
            continue

        y = 1 if label == "B" else 0
        seqs.append(seq)
        ys.append(y)

    return np.array(seqs, dtype=object), np.array(ys, dtype=int)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--order", type=int, required=True, help="Markov order m (0..10)")
    ap.add_argument("--k", type=int, required=True, help="Number of folds (3..5)")
    ap.add_argument("--input", required=True, help="Training chromosome TSV (single file)")
    ap.add_argument("--TF", required=True, choices=["CTCF", "REST", "EP300"], help="TF name")
    ap.add_argument("--genome", default=os.path.join("data", "genome", "hg38.fa"),
                    help="Reference genome FASTA (relative path recommended)")
    ap.add_argument("--outdir", default=os.path.join("results", "milestone1"),
                    help="Output directory for plots and summary")
    ap.add_argument("--seed", type=int, default=42, help="Random seed for CV split")
    args = ap.parse_args()

    if not (0 <= args.order <= 10):
        raise SystemExit("Error: --order must be in [0, 10]")
    if not (3 <= args.k <= 5):
        raise SystemExit("Error: --k must be in [3, 5]")

    os.makedirs(args.outdir, exist_ok=True)

    # Ensure FASTA is indexed (creates .fai if missing)
    if not os.path.exists(args.genome + ".fai"):
        pysam.faidx(args.genome)

    seqs, y = load_sequences_and_labels(args.input, args.genome, args.TF)

    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    print(f"TF={args.TF} | order(m)={args.order} | k={args.k}")
    print(f"Samples={len(y)} | Pos(B)={n_pos} | Neg(U)={n_neg}")
    if n_pos == 0 or n_neg == 0:
        raise SystemExit("Error: No B or U samples found. Check TF column in TSV.")

    skf = StratifiedKFold(n_splits=args.k, shuffle=True, random_state=args.seed)

    roc_aucs = []
    pr_aucs = []

    fold = 1
    for train_idx, test_idx in skf.split(seqs, y):
        train_seqs = seqs[train_idx]
        train_y = y[train_idx]
        test_seqs = seqs[test_idx]
        test_y = y[test_idx]

        pos_train = [s for s, yy in zip(train_seqs, train_y) if yy == 1]
        neg_train = [s for s, yy in zip(train_seqs, train_y) if yy == 0]

        pos_model = train_markov(pos_train, args.order)
        neg_model = train_markov(neg_train, args.order)

        scores = np.array([llr_score(pos_model, neg_model, s, args.order) for s in test_seqs], dtype=float)

        # ROC
        fpr, tpr, _ = roc_curve(test_y, scores)
        roc_auc = auc(fpr, tpr)
        roc_aucs.append(roc_auc)

        plt.figure()
        plt.plot(fpr, tpr)
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.title(f"ROC Fold {fold} | {args.TF} m={args.order} AUC={roc_auc:.4f}")
        roc_path = os.path.join(args.outdir, f"{args.TF}_m{args.order}_roc_fold{fold}.png")
        plt.savefig(roc_path, dpi=150, bbox_inches="tight")
        plt.close()

        # PR
        precision, recall, _ = precision_recall_curve(test_y, scores)
        pr_auc = auc(recall, precision)
        pr_aucs.append(pr_auc)

        plt.figure()
        plt.plot(recall, precision)
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.title(f"PR Fold {fold} | {args.TF} m={args.order} AUC={pr_auc:.4f}")
        pr_path = os.path.join(args.outdir, f"{args.TF}_m{args.order}_pr_fold{fold}.png")
        plt.savefig(pr_path, dpi=150, bbox_inches="tight")
        plt.close()

        print(f"Fold {fold}: ROC-AUC={roc_auc:.4f} | PR-AUC={pr_auc:.4f}")
        fold += 1

    mean_roc = float(np.mean(roc_aucs))
    mean_pr = float(np.mean(pr_aucs))

    print("\nMean ROC-AUC:", mean_roc)
    print("Mean PR-AUC:", mean_pr)

    summary_path = os.path.join(args.outdir, f"summary_{args.TF}_m{args.order}.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(f"TF={args.TF}\n")
        f.write(f"order_m={args.order}\n")
        f.write(f"k={args.k}\n")
        f.write(f"input_tsv={args.input}\n")
        f.write(f"genome_fasta={args.genome}\n")
        f.write(f"samples={len(y)}\n")
        f.write(f"pos_B={n_pos}\n")
        f.write(f"neg_U={n_neg}\n")
        f.write(f"mean_roc_auc={mean_roc}\n")
        f.write(f"mean_pr_auc={mean_pr}\n")

    print("Saved summary:", summary_path)
    print("Saved plots in:", args.outdir)


if __name__ == "__main__":
    main()
