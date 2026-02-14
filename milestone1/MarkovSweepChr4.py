#!/usr/bin/env python3
"""
MarkovSweepChr4.py

Runs MarkovCrossValidation logic for one chromosome TSV across Markov orders m=0..10
with fixed k (default 5) and one TF. Prints:

- 11 mean ROC-AUC values (m=0..10) with 3 decimals
- 11 mean PR-AUC values (m=0..10) with 3 decimals
- 11 timing values (seconds) for computing mean ROC-AUC per m

Also saves ROC/PR plots for EVERY fold for EVERY m (optional; can be heavy).
By default it saves plots only for one selected order via --save_plots_for_m.
If you want plots for all m, pass --save_plots_for_m all

Recommended for the assignment: run on chr4_200bp_bins.tsv, k=5, TF chosen by you.

Example:
python MarkovSweepChr4.py --input data/train/chr4_200bp_bins.tsv --TF CTCF --k 5

Then upload any one of the generated fold plots (ROC and PR) from results/milestone1_sweep/.
"""

import os
import math
import time
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
    if m < 0:
        raise ValueError("order m must be >= 0")
    counts = defaultdict(lambda: defaultdict(int))
    for s in seqs:
        s = s.upper()
        if len(s) <= m or (not is_valid_dna(s)):
            continue
        for i in range(m, len(s)):
            ctx = s[i-m:i] if m > 0 else ""
            nxt = s[i]
            counts[ctx][nxt] += 1
    probs = {}
    for ctx, nxt_counts in counts.items():
        denom = sum(nxt_counts.get(b, 0) + 1 for b in ALPHABET)
        probs[ctx] = {b: (nxt_counts.get(b, 0) + 1)/denom for b in ALPHABET}
    return probs

def log_prob_seq(model_probs, seq: str, m: int) -> float:
    s = seq.upper()
    if len(s) <= m or (not is_valid_dna(s)):
        return float("-inf")
    logp = 0.0
    for i in range(m, len(s)):
        ctx = s[i-m:i] if m > 0 else ""
        nxt = s[i]
        p = model_probs[ctx][nxt] if ctx in model_probs else 0.25
        logp += math.log(p)
    return logp

def llr_score(pos_model, neg_model, seq: str, m: int) -> float:
    return log_prob_seq(pos_model, seq, m) - log_prob_seq(neg_model, seq, m)

def load_sequences_and_labels(tsv_path: str, fasta_path: str, tf_name: str):
    df = pd.read_csv(tsv_path, sep="\t")
    fa = pysam.FastaFile(fasta_path)
    seqs, ys = [], []
    for _, r in df.iterrows():
        label = r[tf_name]
        if label not in ("B", "U"):
            continue
        chrom = r["chr"]
        start = int(r["start"])
        end = int(r["end"])
        seq = fa.fetch(chrom, start, end).upper()
        if len(seq) != (end-start):
            continue
        if not is_valid_dna(seq):
            continue
        ys.append(1 if label == "B" else 0)
        seqs.append(seq)
    return np.array(seqs, dtype=object), np.array(ys, dtype=int)

def run_cv_for_m(seqs, y, m: int, k: int, seed: int, outdir: str, save_plots: bool, tf: str):
    skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=seed)
    roc_aucs, pr_aucs = [], []
    fold = 1
    for train_idx, test_idx in skf.split(seqs, y):
        train_seqs = seqs[train_idx]
        train_y = y[train_idx]
        test_seqs = seqs[test_idx]
        test_y = y[test_idx]

        pos_train = [s for s, yy in zip(train_seqs, train_y) if yy == 1]
        neg_train = [s for s, yy in zip(train_seqs, train_y) if yy == 0]

        pos_model = train_markov(pos_train, m)
        neg_model = train_markov(neg_train, m)

        scores = np.array([llr_score(pos_model, neg_model, s, m) for s in test_seqs], dtype=float)

        fpr, tpr, _ = roc_curve(test_y, scores)
        roc_auc = auc(fpr, tpr)
        roc_aucs.append(roc_auc)

        precision, recall, _ = precision_recall_curve(test_y, scores)
        pr_auc = auc(recall, precision)
        pr_aucs.append(pr_auc)

        if save_plots:
            plt.figure()
            plt.plot(fpr, tpr)
            plt.xlabel("False Positive Rate")
            plt.ylabel("True Positive Rate")
            plt.title(f"ROC Fold {fold} | {tf} m={m} AUC={roc_auc:.4f}")
            plt.savefig(os.path.join(outdir, f"{tf}_m{m}_roc_fold{fold}.png"),
                        dpi=150, bbox_inches="tight")
            plt.close()

            plt.figure()
            plt.plot(recall, precision)
            plt.xlabel("Recall")
            plt.ylabel("Precision")
            plt.title(f"PR Fold {fold} | {tf} m={m} AUC={pr_auc:.4f}")
            plt.savefig(os.path.join(outdir, f"{tf}_m{m}_pr_fold{fold}.png"),
                        dpi=150, bbox_inches="tight")
            plt.close()

        fold += 1

    return float(np.mean(roc_aucs)), float(np.mean(pr_aucs))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Single chromosome training TSV, e.g. data/train/chr4_200bp_bins.tsv")
    ap.add_argument("--TF", required=True, choices=["CTCF", "REST", "EP300"], help="TF to evaluate")
    ap.add_argument("--k", type=int, default=5, help="Number of folds (default 5)")
    ap.add_argument("--genome", default=os.path.join("data", "genome", "hg38.fa"), help="Genome FASTA (relative path)")
    ap.add_argument("--outdir", default=os.path.join("results", "milestone1_sweep"), help="Output directory")
    ap.add_argument("--seed", type=int, default=42, help="CV seed")
    ap.add_argument("--save_plots_for_m", default="6",
                    help="Which m to save plots for: a number 0..10, or 'all' (default 6)")
    args = ap.parse_args()

    if not os.path.exists(args.genome + ".fai"):
        pysam.faidx(args.genome)

    os.makedirs(args.outdir, exist_ok=True)

    seqs, y = load_sequences_and_labels(args.input, args.genome, args.TF)
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    print(f"Input={args.input}")
    print(f"TF={args.TF} | k={args.k} | Samples={len(y)} | Pos(B)={n_pos} | Neg(U)={n_neg}")

    roc_list = []
    pr_list = []
    time_list = []

    for m in range(0, 11):
        save_plots = (args.save_plots_for_m.lower() == "all") or (str(m) == args.save_plots_for_m.strip())
        t0 = time.perf_counter()
        mean_roc, mean_pr = run_cv_for_m(seqs, y, m, args.k, args.seed, args.outdir, save_plots, args.TF)
        t1 = time.perf_counter()
        roc_list.append(mean_roc)
        pr_list.append(mean_pr)
        time_list.append(t1 - t0)
        print(f"m={m}: mean ROC-AUC={mean_roc:.6f} | mean PR-AUC={mean_pr:.6f} | time_sec={t1-t0:.3f}")

    # Print in the exact format needed for form entry
    print("\nROC-AUC (m=0..10), 3 decimals:")
    print(", ".join(f"{x:.3f}" for x in roc_list))
    best_m_roc = int(np.argmax(roc_list))
    print("Best ROC-AUC order:", best_m_roc)

    print("\nPR-AUC (m=0..10), 3 decimals:")
    print(", ".join(f"{x:.3f}" for x in pr_list))
    best_m_pr = int(np.argmax(pr_list))
    print("Best PR-AUC order:", best_m_pr)

    print("\nTime (seconds) for mean ROC-AUC computation per m=0..10:")
    print(", ".join(f"{t:.3f}" for t in time_list))

if __name__ == "__main__":
    main()
