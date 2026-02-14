import math
from collections import defaultdict

ALPHABET = ("A", "C", "G", "T")

def is_valid_dna(seq: str) -> bool:
    seq = seq.upper()
    return all(ch in ALPHABET for ch in seq)

def train_markov(seqs, m: int):
    """
    Train an order-m Markov model with Laplace smoothing.
    Returns: dict mapping (context:str -> dict(base->prob))
    """
    if m < 0:
        raise ValueError("m must be >= 0")

    # counts[context][base] = count
    counts = defaultdict(lambda: defaultdict(int))

    for s in seqs:
        s = s.upper()
        if len(s) <= m:
            continue
        if not is_valid_dna(s):
            continue

        for i in range(m, len(s)):
            ctx = s[i-m:i] if m > 0 else ""
            nxt = s[i]
            counts[ctx][nxt] += 1

    # convert to probabilities with Laplace smoothing
    probs = {}
    for ctx, next_counts in counts.items():
        # Laplace: +1 for each base
        denom = sum(next_counts.get(b, 0) + 1 for b in ALPHABET)
        probs[ctx] = {b: (next_counts.get(b, 0) + 1) / denom for b in ALPHABET}

    # Also handle unseen contexts at scoring time:
    # we will use uniform 1/4 if context never seen
    return probs

def log_prob_seq(model_probs, seq: str, m: int) -> float:
    """
    Compute log P(seq) under an order-m Markov model.
    Uses uniform 1/4 for unseen contexts.
    """
    seq = seq.upper()
    if len(seq) <= m:
        return float("-inf")
    if not is_valid_dna(seq):
        return float("-inf")

    logp = 0.0
    for i in range(m, len(seq)):
        ctx = seq[i-m:i] if m > 0 else ""
        nxt = seq[i]

        if ctx in model_probs:
            p = model_probs[ctx][nxt]
        else:
            p = 0.25  # unseen context fallback

        logp += math.log(p)
    return logp

def llr_score(pos_model, neg_model, seq: str, m: int) -> float:
    return log_prob_seq(pos_model, seq, m) - log_prob_seq(neg_model, seq, m)

