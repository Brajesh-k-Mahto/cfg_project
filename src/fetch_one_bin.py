import pandas as pd
import pysam

TSV = "data/train/chr1_200bp_bins.tsv"   # change if your filename differs
FA  = "data/genome/hg38.fa"

# Pick a row number to test (0 = first data row)
ROW_INDEX = 0

df = pd.read_csv(TSV, sep="\t")
row = df.iloc[ROW_INDEX]

chrom = row["chr"]
start = int(row["start"])
end   = int(row["end"])

fa = pysam.FastaFile(FA)

# pysam uses 0-based, half-open coordinates: [start, end)
# Your TSV is typically already 0-based bins; we use as-is.
seq = fa.fetch(chrom, start, end).upper()

print("Bin:", chrom, start, end)
print("Length:", len(seq))
print("Seq (first 60):", seq[:60])
print("Has_N:", "N" in seq)

