import os
import pandas as pd

# CHANGE THIS to your actual chr1 file name/path
INPUT_TSV = "./data/train/chr1_200bp_bins.tsv"

# Output folder + file
OUT_DIR = "results"
OUT_TSV = os.path.join(OUT_DIR, "preview_first20.tsv")

def main():
    # 1) Read file
    df = pd.read_csv(INPUT_TSV, sep="\t")

    # 2) Print number of rows
    print("Rows (bins):", len(df))

    # 3) Print column names
    print("Columns:", list(df.columns))

    # 4) Count B and U for one TF (CTCF)
    if "CTCF" not in df.columns:
        raise ValueError("CTCF column not found in this file. Check the file/columns.")

    counts = df["CTCF"].value_counts(dropna=False)
    b_count = int(counts.get("B", 0))
    u_count = int(counts.get("U", 0))
    print("CTCF counts -> B:", b_count, ", U:", u_count)

    # 5) Write first 20 rows
    os.makedirs(OUT_DIR, exist_ok=True)
    df.head(20).to_csv(OUT_TSV, sep="\t", index=False)
    print("Wrote first 20 rows to:", OUT_TSV)

if __name__ == "__main__":
    main()

