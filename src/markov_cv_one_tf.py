import os
import argparse
import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt

from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, precision_recall_curve, auc

from markov_core import train_markov, llr_score, is_valid_dna

def load_sequences_and_labels(tsv_path, fasta_path, tf_name):
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
        end   = int(r["end"])
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
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--hg38", required=True)
    ap.add_argument("--tf", required=True, choices=["CTCF", "REST", "EP300"])
    ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", default="results/markov_tf_runs")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    seqs, y = load_sequences_and_labels(args.tsv, args.hg38, args.tf)
    print("TF:", args.tf, "| m:", args.m, "| k:", args.k)
    print("Total samples:", len(y))
    print("Pos (B):", int((y == 1).sum()), "| Neg (U):", int((y == 0).sum()))

    skf = StratifiedKFold(n_splits=args.k, shuffle=True, random_state=args.seed)

    roc_aucs = []
    pr_aucs = []

    fold = 1
    for train_idx, test_idx in skf.split(seqs, y):
        train_seqs = seqs[train_idx]
        train_y    = y[train_idx]
        test_seqs  = seqs[test_idx]
        test_y     = y[test_idx]

        pos_train = [s for s, yy in zip(train_seqs, train_y) if yy == 1]
        neg_train = [s for s, yy in zip(train_seqs, train_y) if yy == 0]

        pos_model = train_markov(pos_train, args.m)
        neg_model = train_markov(neg_train, args.m)

        scores = np.array([llr_score(pos_model, neg_model, s, args.m) for s in test_seqs], dtype=float)

        # ROC
        fpr, tpr, _ = roc_curve(test_y, scores)
        roc_auc = auc(fpr, tpr)
        roc_aucs.append(roc_auc)

        plt.figure()
        plt.plot(fpr, tpr)
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.title(f"ROC Fold {fold} | {args.tf} m={args.m} AUC={roc_auc:.4f}")
        plt.savefig(os.path.join(args.outdir, f"{args.tf}_m{args.m}_roc_fold{fold}.png"),
                    dpi=150, bbox_inches="tight")
        plt.close()

        # PR
        precision, recall, _ = precision_recall_curve(test_y, scores)
        pr_auc = auc(recall, precision)
        pr_aucs.append(pr_auc)

        plt.figure()
        plt.plot(recall, precision)
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.title(f"PR Fold {fold} | {args.tf} m={args.m} AUC={pr_auc:.4f}")
        plt.savefig(os.path.join(args.outdir, f"{args.tf}_m{args.m}_pr_fold{fold}.png"),
                    dpi=150, bbox_inches="tight")
        plt.close()

        print(f"Fold {fold}: ROC-AUC={roc_auc:.4f} | PR-AUC={pr_auc:.4f}")
        fold += 1

    mean_roc = float(np.mean(roc_aucs))
    mean_pr = float(np.mean(pr_aucs))

    print("\nMean ROC-AUC:", mean_roc)
    print("Mean PR-AUC:", mean_pr)

    # Save a one-line summary
    summary_path = os.path.join(args.outdir, f"summary_{args.tf}_m{args.m}.txt")
    with open(summary_path, "w") as f:
        f.write(f"TF={args.tf}\n")
        f.write(f"m={args.m}\n")
        f.write(f"k={args.k}\n")
        f.write(f"mean_roc_auc={mean_roc}\n")
        f.write(f"mean_pr_auc={mean_pr}\n")
    print("Saved:", summary_path)

if __name__ == "__main__":
    main()

