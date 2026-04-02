import argparse
import glob
import os

import pandas as pd
import pysam

ALPHABET = {"A", "C", "G", "T"}


def is_valid_dna(seq: str) -> bool:
    return all(base in ALPHABET for base in seq)


def write_unknown_fasta(unknown_tsv, fasta_path, out_fasta):
    fa = pysam.FastaFile(fasta_path)
    df = pd.read_csv(unknown_tsv, sep="\t")

    written = 0
    with open(out_fasta, "w", encoding="utf-8") as handle:
        for _, row in df.iterrows():
            chrom = row["chr"]
            start = int(row["start"])
            end = int(row["end"])
            seq = fa.fetch(chrom, start, end).upper()

            if len(seq) != (end - start):
                continue
            if not is_valid_dna(seq):
                continue

            # Keep header deterministic so it can be matched in FIMO output.
            header = f"{chrom}:{start}-{end}"
            handle.write(f">{header}\n{seq}\n")
            written += 1

    return written


def main():
    parser = argparse.ArgumentParser(
        description="Prepare unknown chromosome FASTA files for FIMO scanning."
    )
    parser.add_argument("--unknown_dir", default="data/unknown", help="Directory with unknown chr TSV files")
    parser.add_argument("--hg38", required=True, help="Path to hg38 FASTA")
    parser.add_argument(
        "--out_dir",
        default=os.path.join("analysis", "meme", "outputs", "fimo_input"),
        help="Output directory for unknown FASTA files",
    )
    args = parser.parse_args()

    unknown_files = sorted(glob.glob(os.path.join(args.unknown_dir, "chr*_200bp_bins_unknown.tsv")))
    if not unknown_files:
        raise FileNotFoundError(f"No unknown TSV files found under: {args.unknown_dir}")

    os.makedirs(args.out_dir, exist_ok=True)

    for uf in unknown_files:
        chrom = os.path.basename(uf).replace("_200bp_bins_unknown.tsv", "")
        out_fasta = os.path.join(args.out_dir, f"{chrom}_unknown.fa")
        written = write_unknown_fasta(uf, args.hg38, out_fasta)
        print(f"{chrom}: wrote {written} FASTA records -> {out_fasta}")


if __name__ == "__main__":
    main()
