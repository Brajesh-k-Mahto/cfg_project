import argparse
import glob
import os

import pandas as pd
import pysam

ALPHABET = {"A", "C", "G", "T"}


def is_valid_dna(seq: str) -> bool:
    return all(base in ALPHABET for base in seq)


def fetch_bound_sequences(train_files, fasta_path, tf_name):
    fa = pysam.FastaFile(fasta_path)
    out = []

    for fp in train_files:
        df = pd.read_csv(fp, sep="\t")
        for _, row in df.iterrows():
            label = row[tf_name]
            if label != "B":
                continue

            chrom = row["chr"]
            start = int(row["start"])
            end = int(row["end"])
            seq = fa.fetch(chrom, start, end).upper()

            if len(seq) != (end - start):
                continue
            if not is_valid_dna(seq):
                continue

            header = f"{chrom}:{start}-{end}"
            out.append((header, seq))

    return out


def write_fasta(records, out_fasta):
    os.makedirs(os.path.dirname(out_fasta), exist_ok=True)
    with open(out_fasta, "w", encoding="utf-8") as handle:
        for header, seq in records:
            handle.write(f">{header}\n{seq}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Prepare TF-specific bound FASTA for MEME motif discovery."
    )
    parser.add_argument("--train_dir", default="data/train", help="Directory with labeled chr TSV files")
    parser.add_argument("--hg38", required=True, help="Path to hg38 FASTA")
    parser.add_argument("--tf", required=True, choices=["CTCF", "REST", "EP300"], help="Target TF")
    parser.add_argument(
        "--out_fasta",
        default=None,
        help="Output FASTA path; default uses analysis/meme/outputs/meme_input/{TF}_bound.fa",
    )
    args = parser.parse_args()

    train_files = sorted(glob.glob(os.path.join(args.train_dir, "chr*_200bp_bins.tsv")))
    if not train_files:
        raise FileNotFoundError(f"No training TSV files found under: {args.train_dir}")

    out_fasta = args.out_fasta or os.path.join(
        "analysis", "meme", "outputs", "meme_input", f"{args.tf}_bound.fa"
    )

    records = fetch_bound_sequences(train_files, args.hg38, args.tf)
    if not records:
        raise RuntimeError(f"No valid bound sequences found for TF={args.tf}")

    write_fasta(records, out_fasta)
    print(f"TF: {args.tf}")
    print(f"Training files scanned: {len(train_files)}")
    print(f"Bound FASTA records: {len(records)}")
    print(f"Wrote: {out_fasta}")


if __name__ == "__main__":
    main()
