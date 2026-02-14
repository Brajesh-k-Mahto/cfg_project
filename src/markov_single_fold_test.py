import pandas as pd
import pysam
from markov_core import train_markov, llr_score, is_valid_dna

TSV = "data/train/chr1_200bp_bins.tsv"      # adjust if needed
FA  = "data/genome/hg38.fa"
TF  = "CTCF"
M   = 3   # Markov order (change later)

def main():
    df = pd.read_csv(TSV, sep="\t")
    fa = pysam.FastaFile(FA)

    pos_seqs = []
    neg_seqs = []

    for _, r in df.iterrows():
        chrom = r["chr"]
        start = int(r["start"])
        end   = int(r["end"])
        label = r[TF]

        seq = fa.fetch(chrom, start, end).upper()

        if len(seq) != (end - start):
            continue
        if not is_valid_dna(seq):
            continue

        if label == "B":
            pos_seqs.append(seq)
        elif label == "U":
            neg_seqs.append(seq)

    print("TF:", TF, "| m:", M)
    print("Pos (B) sequences:", len(pos_seqs))
    print("Neg (U) sequences:", len(neg_seqs))

    if len(pos_seqs) == 0 or len(neg_seqs) == 0:
        raise RuntimeError("No sequences found for B or U. Check labels and TF name.")

    pos_model = train_markov(pos_seqs, M)
    neg_model = train_markov(neg_seqs, M)

    # Score a few examples
    print("\nSample scores (higher => more likely Bound):")
    for i in range(5):
        s = pos_seqs[i]
        print("B score:", llr_score(pos_model, neg_model, s, M))

    for i in range(5):
        s = neg_seqs[i]
        print("U score:", llr_score(pos_model, neg_model, s, M))

if __name__ == "__main__":
    main()

