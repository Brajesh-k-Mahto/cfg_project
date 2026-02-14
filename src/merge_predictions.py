import os
import pandas as pd

PRED_DIR = "results/markov_predictions"
UNKNOWN_DIR = "data/unknown"
OUT_DIR = "results/final_predictions"
os.makedirs(OUT_DIR, exist_ok=True)

def load_pred(chrom, tf):
    path = os.path.join(PRED_DIR, f"{chrom}_{tf}_pred.tsv")
    df = pd.read_csv(path, sep="\t")
    return df[["chr", "start", "end", "pred"]].rename(columns={"pred": f"{tf}_pred"})

def merge_one(chrom):
    unk_path = os.path.join(UNKNOWN_DIR, f"{chrom}_200bp_bins_unknown.tsv")
    base = pd.read_csv(unk_path, sep="\t")  # chr, start, end, ATAC

    out = base.copy()
    for tf in ["CTCF", "REST", "EP300"]:
        p = load_pred(chrom, tf)
        out = out.merge(p, on=["chr", "start", "end"], how="left")

    out_path = os.path.join(OUT_DIR, f"final_{chrom}_predictions.tsv")
    out.to_csv(out_path, sep="\t", index=False)
    print("Wrote:", out_path)

def main():
    for chrom in ["chr3", "chr10", "chr17"]:
        merge_one(chrom)

if __name__ == "__main__":
    main()

