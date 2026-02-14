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
K = 5
SEED = 42
M_VALUES = list(range(0, 11))  # 0..10

OUT_DIR = "results/markov_m_sweep"
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

def cv_auc_for_m(seqs, y, m, k=5, seed=42):
    skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=seed)
    roc_aucs = []
    pr_aucs = []

    for train_idx, test_idx in skf.split(seqs, y):
        train_seqs = seqs[train_idx]
        train_y    = y[train_idx]
        test_seqs  = seqs[test_idx]
        test_y     = y[test_idx]

        pos_train = [s for s, yy in zip(train_seqs, train_y) if yy == 1]
        neg_train = [s for s, yy in zip(train_seqs, train_y) if yy == 0]

        pos_model = train_markov(pos_train, m)
        neg_model = train_markov(neg_train, m)

        scores = np.array([llr_score(pos_model, neg_model, s, m) for s in test_seqs], dtype=float)

        fpr, tpr, _ = roc_curve(test_y, scores)
        roc_aucs.append(auc(fpr, tpr))

        precision, recall, _ = precision_recall_curve(test_y, scores)
        pr_aucs.append(auc(recall, precision))

    return float(np.mean(roc_aucs)), float(np.mean(pr_aucs))

def main():
    seqs, y = load_sequences_and_labels(TSV, FA, TF)
    print("Total samples:", len(y))
    print("Pos (B):", int((y == 1).sum()), "| Neg (U):", int((y == 0).sum()))

    rows = []
    for m in M_VALUES:
        mean_roc, mean_pr = cv_auc_for_m(seqs, y, m, k=K, seed=SEED)
        print(f"m={m}: Mean ROC-AUC={mean_roc:.6f} | Mean PR-AUC={mean_pr:.6f}")
        rows.append({"m": m, "mean_roc_auc": mean_roc, "mean_pr_auc": mean_pr})

    res = pd.DataFrame(rows)
    out_csv = os.path.join(OUT_DIR, f"m_sweep_{TF}.csv")
    res.to_csv(out_csv, index=False)
    print("Saved table:", out_csv)

    # Plot ROC-AUC vs m
    plt.figure()
    plt.plot(res["m"], res["mean_roc_auc"])
    plt.xlabel("Markov order (m)")
    plt.ylabel("Mean ROC-AUC")
    plt.title(f"Markov order sweep | TF={TF}")
    plt.savefig(os.path.join(OUT_DIR, f"roc_auc_vs_m_{TF}.png"), dpi=150, bbox_inches="tight")
    plt.close()

    # Plot PR-AUC vs m
    plt.figure()
    plt.plot(res["m"], res["mean_pr_auc"])
    plt.xlabel("Markov order (m)")
    plt.ylabel("Mean PR-AUC")
    plt.title(f"Markov order sweep | TF={TF}")
    plt.savefig(os.path.join(OUT_DIR, f"pr_auc_vs_m_{TF}.png"), dpi=150, bbox_inches="tight")
    plt.close()

    best_roc_row = res.iloc[res["mean_roc_auc"].idxmax()]
    best_pr_row  = res.iloc[res["mean_pr_auc"].idxmax()]

    print("\nBest by ROC-AUC:", dict(best_roc_row))
    print("Best by PR-AUC :", dict(best_pr_row))
    print("Saved plots in:", OUT_DIR)

if __name__ == "__main__":
    main()

