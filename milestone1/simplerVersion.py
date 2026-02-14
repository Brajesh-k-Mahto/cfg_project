#!/usr/bin/env python3
"""
simplerVersion.py

Takes only TWO arguments:
  1) FASTA file containing DNA sequences
  2) Markov order m (integer >= 0)

Trains an order-m Markov model (Laplace-smoothed) on ALL sequences in the FASTA file,
then prints the log-likelihood of EACH sequence under the trained model (one number per line).

Usage:
python simplerVersion.py <sequences.fasta> <m>

Notes:
- FASTA records may have mixed case; only A/C/G/T are used.
- Any sequence containing other letters (including N) is skipped from training AND scored as -inf.
"""

import sys
import math
from collections import defaultdict

ALPHABET = ("A", "C", "G", "T")

def is_valid_dna(seq: str) -> bool:
    s = seq.upper()
    return all(ch in ALPHABET for ch in s)

def read_fasta(path: str):
    header = None
    seq_chunks = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_chunks)
                header = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header is not None:
            yield header, "".join(seq_chunks)

def train_markov(seqs, m: int):
    counts = defaultdict(lambda: defaultdict(int))
    for s in seqs:
        s = s.upper()
        if len(s) <= m or (not is_valid_dna(s)):
            continue
        for i in range(m, len(s)):
            ctx = s[i-m:i] if m > 0 else ""
            nxt = s[i]
            counts[ctx][nxt] += 1

    probs = {}
    for ctx, nxt_counts in counts.items():
        denom = sum(nxt_counts.get(b, 0) + 1 for b in ALPHABET)
        probs[ctx] = {b: (nxt_counts.get(b, 0) + 1)/denom for b in ALPHABET}
    return probs

def log_prob_seq(model_probs, seq: str, m: int) -> float:
    s = seq.upper()
    if len(s) <= m or (not is_valid_dna(s)):
        return float("-inf")
    logp = 0.0
    for i in range(m, len(s)):
        ctx = s[i-m:i] if m > 0 else ""
        nxt = s[i]
        p = model_probs[ctx][nxt] if ctx in model_probs else 0.25
        logp += math.log(p)
    return logp

def main():
    if len(sys.argv) != 3:
        print("Usage: python simplerVersion.py <sequences.fasta> <m>", file=sys.stderr)
        sys.exit(1)

    fasta_path = sys.argv[1]
    m = int(sys.argv[2])
    if m < 0:
        raise SystemExit("m must be >= 0")

    seqs = [seq for _, seq in read_fasta(fasta_path)]
    model = train_markov(seqs, m)

    for _, seq in read_fasta(fasta_path):
        print(log_prob_seq(model, seq, m))

if __name__ == "__main__":
    main()
