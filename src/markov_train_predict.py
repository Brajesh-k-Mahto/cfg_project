import os
import glob
import argparse
import pandas as pd
import pysam

from markov_core import train_markov, llr_score, is_valid_dna

def load_labeled_seqs(train_files, fasta_path, tf):
    fa = pysam.FastaFile(fasta_path)
    pos, neg = [], []

    for fp in train_files:
        df = pd.read_csv(fp, sep="\t")
        for _, r in df.iterrows():
            label = r[tf]
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

            if label == "B":
                pos.append(seq)
            else:
                neg.append(seq)

    return pos, neg

def predict_file(unknown_tsv, fasta_path, pos_model, neg_model, tf, m, threshold, out_path):
    fa = pysam.FastaFile(fasta_path)
    df = pd.read_csv(unknown_tsv, sep="\t")

    out_rows = []
    for _, r in df.iterrows():
        chrom = r["chr"]
        start = int(r["start"])
        end   = int(r["end"])

        seq = fa.fetch(chrom, start, end).upper()
        if len(seq) != (end - start) or (not is_valid_dna(seq)):
            score = float("nan")
            pred = "U"  # conservative default
        else:
            score = llr_score(pos_model, neg_model, seq, m)
            pred = "B" if score >= threshold else "U"

        out_rows.append({
            "chr": chrom,
            "start": start,
            "end": end,
            "TF": tf,
            "score": score,
            "pred": pred
        })

    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(out_path, sep="\t", index=False)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train_dir", default="data/train")
    ap.add_argument("--unknown_dir", default="data/unknown")
    ap.add_argument("--hg38", required=True)
    ap.add_argument("--tf", required=True, choices=["CTCF", "REST", "EP300"])
    ap.add_argument("--m", type=int, default=6)
    ap.add_argument("--threshold", type=float, default=0.0)
    ap.add_argument("--outdir", default="results/markov_predictions")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    train_files = sorted(glob.glob(os.path.join(args.train_dir, "chr*_200bp_bins.tsv")))
    unknown_files = sorted(glob.glob(os.path.join(args.unknown_dir, "chr*_200bp_bins_unknown.tsv")))

    print("TF:", args.tf, "| m:", args.m, "| threshold:", args.threshold)
    print("Training files:", len(train_files))
    print("Unknown files:", len(unknown_files))

    # 1) Train models on all labelled data
    pos_seqs, neg_seqs = load_labeled_seqs(train_files, args.hg38, args.tf)
    print("Training sequences -> B:", len(pos_seqs), "| U:", len(neg_seqs))

    pos_model = train_markov(pos_seqs, args.m)
    neg_model = train_markov(neg_seqs, args.m)

    # 2) Predict each unknown chromosome file
    for uf in unknown_files:
        base = os.path.basename(uf).replace("_200bp_bins_unknown.tsv", "")
        out_path = os.path.join(args.outdir, f"{base}_{args.tf}_pred.tsv")
        predict_file(uf, args.hg38, pos_model, neg_model, args.tf, args.m, args.threshold, out_path)
        print("Wrote:", out_path)

if __name__ == "__main__":
    main()

