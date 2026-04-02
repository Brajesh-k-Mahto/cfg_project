import argparse
import math
import os

import pandas as pd


def transform_score(row, metric):
    if metric == "score":
        return float(row["score"])
    if metric == "neglog10_p":
        p = max(float(row["p-value"]), 1e-300)
        return -math.log10(p)
    if metric == "neglog10_q":
        q = max(float(row["q-value"]), 1e-300)
        return -math.log10(q)
    raise ValueError(f"Unsupported metric: {metric}")


def load_fimo(fimo_tsv):
    if not os.path.exists(fimo_tsv):
        raise FileNotFoundError(f"FIMO file not found: {fimo_tsv}")

    # FIMO output usually contains comment lines starting with '#'.
    rows = []
    with open(fimo_tsv, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            rows.append(line)

    if not rows:
        return pd.DataFrame(
            columns=["motif_id", "motif_alt_id", "sequence_name", "start", "stop", "strand", "score", "p-value", "q-value", "matched_sequence"]
        )

    from io import StringIO

    return pd.read_csv(StringIO("".join(rows)), sep="\t")


def build_bin_index(unknown_df):
    # Keyed by (chrom, bin_start) for constant-time lookup of valid bins.
    index = {}
    for _, row in unknown_df.iterrows():
        index[(row["chr"], int(row["start"]))] = f"{row['chr']}:{int(row['start'])}-{int(row['end'])}"
    return index


def hit_to_sequence_name(row, bin_index):
    seq_name = str(row["sequence_name"])

    # Case 1: FASTA header preserved (e.g., chr3:12200-12400).
    if ":" in seq_name and "-" in seq_name:
        return seq_name

    # Case 2: FIMO parsed genomic coordinates (sequence_name=chr, start/stop are genomic).
    chrom = seq_name
    hit_start = int(row["start"])
    bin_start = (hit_start // 200) * 200
    return bin_index.get((chrom, bin_start))


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate FIMO hits into one numeric TF score per 200bp bin."
    )
    parser.add_argument("--unknown_tsv", required=True, help="Unknown chromosome TSV (chr/start/end/ATAC)")
    parser.add_argument("--fimo_tsv", required=True, help="FIMO TSV output for one TF and one chromosome")
    parser.add_argument("--tf", required=True, choices=["CTCF", "REST", "EP300"], help="Target TF")
    parser.add_argument(
        "--metric",
        default="score",
        choices=["score", "neglog10_p", "neglog10_q"],
        help="How to convert hit-level FIMO results to numeric hit score before max aggregation",
    )
    parser.add_argument("--out_tsv", required=True, help="Output TSV path")
    args = parser.parse_args()

    unknown = pd.read_csv(args.unknown_tsv, sep="\t")
    unknown["sequence_name"] = unknown.apply(
        lambda r: f"{r['chr']}:{int(r['start'])}-{int(r['end'])}", axis=1
    )

    fimo = load_fimo(args.fimo_tsv)

    score_col = args.tf
    if fimo.empty:
        score_map = {}
    else:
        bin_index = build_bin_index(unknown)
        fimo["agg_score"] = fimo.apply(lambda r: transform_score(r, args.metric), axis=1)
        fimo["sequence_name"] = fimo.apply(lambda r: hit_to_sequence_name(r, bin_index), axis=1)
        fimo = fimo.dropna(subset=["sequence_name"])
        fimo = fimo[["sequence_name", "agg_score"]]
        fimo = fimo.groupby("sequence_name", as_index=False).agg({"agg_score": "max"})
        score_map = dict(zip(fimo["sequence_name"], fimo["agg_score"]))

    unknown[score_col] = unknown["sequence_name"].map(score_map).fillna(0.0)
    out = unknown[["chr", "start", "end", "ATAC", score_col]]

    os.makedirs(os.path.dirname(args.out_tsv), exist_ok=True)
    out.to_csv(args.out_tsv, sep="\t", index=False)

    print(f"TF: {args.tf}")
    print(f"Metric: {args.metric}")
    print(f"Rows in unknown file: {len(out)}")
    print(f"Rows with non-zero score: {(out[score_col] > 0).sum()}")
    print(f"Wrote: {args.out_tsv}")


if __name__ == "__main__":
    main()
