import argparse
import os

import pandas as pd


REQUIRED_TFS = ["CTCF", "REST", "EP300"]


def read_tf_score_file(path, tf):
    df = pd.read_csv(path, sep="\t")
    expected_cols = {"chr", "start", "end", "ATAC", tf}
    missing = expected_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
    return df[["chr", "start", "end", "ATAC", tf]]


def main():
    parser = argparse.ArgumentParser(
        description="Merge TF-specific MEME/FIMO scores into final submission-style TSV(.gz)."
    )
    parser.add_argument("--ctcf_tsv", required=True, help="Processed CTCF score TSV")
    parser.add_argument("--rest_tsv", required=True, help="Processed REST score TSV")
    parser.add_argument("--ep300_tsv", required=True, help="Processed EP300 score TSV")
    parser.add_argument("--out_tsv", required=True, help="Output file path; use .tsv.gz to gzip-compress")
    args = parser.parse_args()

    ctcf = read_tf_score_file(args.ctcf_tsv, "CTCF")
    rest = read_tf_score_file(args.rest_tsv, "REST")
    ep300 = read_tf_score_file(args.ep300_tsv, "EP300")

    out = ctcf.merge(rest[["chr", "start", "end", "REST"]], on=["chr", "start", "end"], how="inner")
    out = out.merge(ep300[["chr", "start", "end", "EP300"]], on=["chr", "start", "end"], how="inner")

    # Keep expected final-column order.
    out = out[["chr", "start", "end", "ATAC", "CTCF", "REST", "EP300"]]

    os.makedirs(os.path.dirname(args.out_tsv), exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)

    print(f"Merged rows: {len(out)}")
    print(f"Wrote: {args.out_tsv}")


if __name__ == "__main__":
    main()
