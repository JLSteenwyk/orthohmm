"""Multi-sequence profile HMM construction from MSAs.

Builds position-specific scoring profiles (PSSMs) from multiple
sequence alignments, for iterative profile-based search.

Pipeline:
  1. Align orthogroup members with the built-in center-star MSA
  2. Compute position frequency matrix with BLOSUM pseudocounts
  3. Convert to log-odds emission scores
  4. Use with existing batch_viterbi_c infrastructure
"""

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

from .matrices import ALPHABET, ALPHABET_SIZE, FULL_CHAR_TO_IDX
from .msa_center_star import center_star_msa


@dataclass
class MSAProfile:
    """MSA-derived position-specific profile."""
    length: int                      # number of match states (non-gap-majority columns)
    match_emissions: np.ndarray      # (L, 20) int8 log-odds emission scores
    insert_emissions: np.ndarray     # (20,) int8 background scores
    transitions: np.ndarray          # (7,) int32
    consensus: np.ndarray            # uint8 consensus sequence for prefilter use
    member_ids: List[str] = field(default_factory=list)  # source genes


def encode_alignment(aligned_seqs: List[str]) -> np.ndarray:
    """Encode aligned sequences as int8 array with 20=gap.

    Returns (n_seqs, alignment_length) int8.
    """
    if not aligned_seqs:
        return np.empty((0, 0), dtype=np.int8)

    L = len(aligned_seqs[0])
    n = len(aligned_seqs)
    enc = np.full((n, L), ALPHABET_SIZE, dtype=np.int8)  # 20 = gap/unknown

    for i, seq in enumerate(aligned_seqs):
        for j, c in enumerate(seq):
            if c == "-" or c == ".":
                enc[i, j] = 20  # gap
            else:
                enc[i, j] = FULL_CHAR_TO_IDX.get(ord(c), 20)
    return enc


def compute_pssm(
    encoded_msa: np.ndarray,
    sub_matrix: np.ndarray,
    bg_freqs: np.ndarray,
    pseudocount_weight: float = 5.0,
    gap_fraction_threshold: float = 0.5,
) -> tuple:
    """Compute position-specific scoring matrix from encoded MSA.

    Parameters
    ----------
    encoded_msa : (n_seqs, L) int8, values 0-19 for amino acids, 20 for gap
    sub_matrix : (20, 20) int8 substitution matrix (for pseudocounts)
    bg_freqs : (20,) float64 background frequencies
    pseudocount_weight : BLOSUM pseudocount weight
    gap_fraction_threshold : drop columns with more gaps than this

    Returns
    -------
    pssm : (L_match, 20) int8 log-odds scores
    match_mask : (L,) bool, True for match columns (retained in PSSM)
    consensus : (L_match,) uint8 consensus sequence
    """
    n, L = encoded_msa.shape

    # Compute gap fraction per column
    gap_counts = (encoded_msa == 20).sum(axis=0)
    gap_frac = gap_counts / n
    match_mask = gap_frac < gap_fraction_threshold

    L_match = int(match_mask.sum())
    if L_match == 0:
        return (np.empty((0, 20), dtype=np.int8),
                match_mask,
                np.empty(0, dtype=np.uint8))

    pssm = np.zeros((L_match, 20), dtype=np.int8)
    consensus = np.zeros(L_match, dtype=np.uint8)

    # Convert sub_matrix to conditional probability-like pseudocount:
    # For each observed AA a, pseudocount for AA b is 2^(sub_matrix[a,b]/2) * bg[b]
    # This is BLOSUM pseudocount formulation
    sub_pseudo = np.zeros((20, 20), dtype=np.float64)
    for a in range(20):
        for b in range(20):
            sub_pseudo[a, b] = (2.0 ** (sub_matrix[a, b] / 2.0)) * bg_freqs[b]
        # Normalize row
        s = sub_pseudo[a].sum()
        if s > 0:
            sub_pseudo[a] /= s

    match_idx = 0
    for col in range(L):
        if not match_mask[col]:
            continue

        # Position frequency: count of each AA at this column
        pfm = np.zeros(20, dtype=np.float64)
        for a in range(20):
            pfm[a] = (encoded_msa[:, col] == a).sum()

        n_aa = pfm.sum()  # non-gap count
        if n_aa == 0:
            # All gaps — use background
            probs = bg_freqs.copy()
        else:
            # Empirical probability
            emp_prob = pfm / n_aa

            # Pseudocount: mix observed with substitution-matrix-derived distribution
            pseudo_prob = np.zeros(20, dtype=np.float64)
            for a in range(20):
                if pfm[a] > 0:
                    weight = pfm[a] / n_aa
                    pseudo_prob += weight * sub_pseudo[a]
            if pseudo_prob.sum() == 0:
                pseudo_prob = bg_freqs.copy()

            # Combine with pseudocount weight (alpha=n, beta=pseudocount_weight)
            # p(b) = (n*emp_prob(b) + alpha*pseudo_prob(b)) / (n + alpha)
            alpha = pseudocount_weight
            probs = (n_aa * emp_prob + alpha * pseudo_prob) / (n_aa + alpha)
            probs = probs / probs.sum()  # normalize

        # Log-odds vs background
        for b in range(20):
            if probs[b] > 0 and bg_freqs[b] > 0:
                score = 2.0 * np.log2(probs[b] / bg_freqs[b])
                pssm[match_idx, b] = np.clip(round(score), -128, 127)
            else:
                pssm[match_idx, b] = -10

        # Consensus = most frequent non-gap AA
        consensus[match_idx] = int(np.argmax(probs))
        match_idx += 1

    return pssm, match_mask, consensus


def build_msa_profile(
    sequences: List[str],
    ids: List[str],
    sub_matrix: np.ndarray,
    bg_freqs: np.ndarray,
    gap_open: int = -12,
    gap_extend: int = -3,
    align_threads: int = 1,
) -> Optional[MSAProfile]:
    """Build a multi-sequence profile from a list of sequences.

    Uses the built-in center-star MSA (banded Needleman-Wunsch + projection)
    — no external aligner dependency.

    Returns None if profile can't be built (too few sequences, etc.).
    """
    if len(sequences) < 2:
        return None

    aligned = center_star_msa(sequences, sub_matrix,
                              gap_open=gap_open, gap_extend=gap_extend,
                              n_threads=align_threads)
    if aligned is None:
        return None

    # Encode
    enc = encode_alignment(aligned)
    if enc.size == 0:
        return None

    # Compute PSSM
    pssm, match_mask, consensus = compute_pssm(enc, sub_matrix, bg_freqs)
    if pssm.shape[0] == 0:
        return None

    # Standard HMM parameters
    insert_emissions = np.full(ALPHABET_SIZE, -1, dtype=np.int8)
    gap_close = max(gap_extend + 1, -1)
    transitions = np.array([
        0, gap_open, gap_open, gap_close, gap_extend, gap_close, gap_extend
    ], dtype=np.int32)

    return MSAProfile(
        length=pssm.shape[0],
        match_emissions=pssm,
        insert_emissions=insert_emissions,
        transitions=transitions,
        consensus=consensus,
        member_ids=list(ids),
    )
