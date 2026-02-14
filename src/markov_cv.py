import os
import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt

from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, precision_recall_curve, auc

from markov_core import train_markov, llr_score, is_valid_dna

TSV = "data/train/chr1_200bp_bins.tsv"   # adjust if needed
FA  = "data/genome/hg38.fa"

TF = "CTCF"
M  = 3
K  = 5
SEED = 42

OUT_DIR = "results/markov_cv"
os.makedirs(OUT_DIR, exist_ok=True)

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
    seqs, y = load_sequences_and_labels(TSV, FA, TF)
    print("Total samples:", len(y))
    print("Pos (B):", int((y == 1).sum()), "| Neg (U):", int((y == 0).sum()))

    skf = StratifiedKFold(n_splits=K, shuffle=True, random_state=SEED)

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

        pos_model = train_markov(pos_train, M)
        neg_model = train_markov(neg_train, M)

        scores = np.array([llr_score(pos_model, neg_model, s, M) for s in test_seqs], dtype=float)

        # ROC
        fpr, tpr, _ = roc_curve(test_y, scores)
        roc_auc = auc(fpr, tpr)
        roc_aucs.append(roc_auc)

        plt.figure()
        plt.plot(fpr, tpr)
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.title(f"ROC Fold {fold} | TF={TF} m={M} AUC={roc_auc:.4f}")
        roc_path = os.path.join(OUT_DIR, f"roc_fold{fold}.png")
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
        plt.title(f"PR Fold {fold} | TF={TF} m={M} AUC={pr_auc:.4f}")
        pr_path = os.path.join(OUT_DIR, f"pr_fold{fold}.png")
        plt.savefig(pr_path, dpi=150, bbox_inches="tight")
        plt.close()

        print(f"Fold {fold}: ROC-AUC={roc_auc:.4f} | PR-AUC={pr_auc:.4f}")
        fold += 1

    print("\nMean ROC-AUC:", float(np.mean(roc_aucs)))
    print("Mean PR-AUC:", float(np.mean(pr_aucs)))
    print("Saved plots in:", OUT_DIR)

if __name__ == "__main__":
    main()

