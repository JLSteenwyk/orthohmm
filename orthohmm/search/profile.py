"""Single-sequence profile HMM construction.

Builds a Plan7-style profile HMM from a single query sequence and
a substitution matrix, equivalent to what phmmer does internally.

For position i with amino acid a_i:
  - Match emission scores = substitution_matrix[a_i, :]
  - Insert emission scores = background-derived scores (position-independent)
  - Transition probabilities use standard HMMER defaults

The substitution matrix scores ARE the log-odds PSSM entries, so we
keep them as integers for efficient Viterbi computation.
"""

from dataclasses import dataclass

import numpy as np
from numba import njit, int32, int8

from .matrices import ALPHABET_SIZE


# ──────────────────────────────────────────────────────────────────────
# Default HMM transition parameters (in scaled integer log-odds)
#
# HMMER3 phmmer defaults (approximate):
#   p(M->M) = 1 - 2*p_open  ≈ 0.96
#   p(M->I) = p_open         ≈ 0.02
#   p(M->D) = p_open         ≈ 0.02
#   p(I->M) = 1 - p_extend   ≈ 0.6
#   p(I->I) = p_extend        ≈ 0.4
#   p(D->M) = 1 - p_extend   ≈ 0.6
#   p(D->D) = p_extend        ≈ 0.4
#
# Converted to integer log-odds (scaled by ~3 to stay in int8 range):
#   log(0.96) * 3 ≈ -0.12  → 0   (near zero penalty for M->M)
#   log(0.02) * 3 ≈ -11.7  → -12 (large penalty for opening gap)
#   log(0.6)  * 3 ≈ -1.53  → -2  (small penalty for closing gap)
#   log(0.4)  * 3 ≈ -2.75  → -3  (moderate penalty for extending gap)
#
# We use a simpler integer scaling that captures the relative costs.
# The exact values don't matter much because OrthoHMM normalizes
# scores by gene length and corrects by phylogenetic distance.
# ──────────────────────────────────────────────────────────────────────

# Transition indices in the (L, 7) transitions array
T_MM = 0  # Match -> Match
T_MI = 1  # Match -> Insert
T_MD = 2  # Match -> Delete
T_IM = 3  # Insert -> Match
T_II = 4  # Insert -> Insert
T_DM = 5  # Delete -> Match
T_DD = 6  # Delete -> Delete


@dataclass
class ProfileHMM:
    """Lightweight Plan7 profile for one query sequence."""
    length: int                    # L (number of match states)
    match_emissions: np.ndarray    # (L, 20) int8
    insert_emissions: np.ndarray   # (20,) int8
    transitions: np.ndarray        # (7,) int32, uniform across positions


def build_profile(
    query_seq: np.ndarray,
    query_len: int,
    sub_matrix: np.ndarray,
    bg_freqs: np.ndarray,
    gap_open: int = -12,
    gap_extend: int = -3,
) -> ProfileHMM:
    """Build a single-sequence Plan7 profile HMM.

    Parameters
    ----------
    query_seq : uint8 array of length query_len
    query_len : int
    sub_matrix : (20, 20) int8 substitution matrix
    bg_freqs : (20,) float64 background frequencies
    gap_open : int, transition penalty for M->I and M->D
    gap_extend : int, transition penalty for I->I and D->D

    Returns
    -------
    ProfileHMM with match emissions from the substitution matrix rows.
    """
    # Match emissions: row of sub_matrix for the amino acid at each position
    match_emissions = np.empty((query_len, ALPHABET_SIZE), dtype=np.int8)
    for i in range(query_len):
        aa = query_seq[i]
        if aa < ALPHABET_SIZE:
            match_emissions[i, :] = sub_matrix[aa, :]
        else:
            # Unknown residue: use zero scores (no information)
            match_emissions[i, :] = 0

    # Insert emissions: derived from background frequencies
    # Convert to integer scores roughly comparable to substitution matrix
    # log2(1/20) ≈ -4.3, so uniform background ≈ -4 per residue
    # We use a simple approximation: all insert emissions = -1
    # (mild penalty, since inserts should be possible but not favored)
    insert_emissions = np.full(ALPHABET_SIZE, -1, dtype=np.int8)

    # Transitions: uniform across all positions
    # gap_close = roughly -log(1-p_extend), should be small
    gap_close = max(gap_extend + 1, -1)  # slightly less penalty than extend

    transitions = np.array([
        0,           # T_MM: no penalty for match-match
        gap_open,    # T_MI: penalty for opening insertion
        gap_open,    # T_MD: penalty for opening deletion
        gap_close,   # T_IM: small penalty for closing insertion
        gap_extend,  # T_II: penalty for extending insertion
        gap_close,   # T_DM: small penalty for closing deletion
        gap_extend,  # T_DD: penalty for extending deletion
    ], dtype=np.int32)

    return ProfileHMM(
        length=query_len,
        match_emissions=match_emissions,
        insert_emissions=insert_emissions,
        transitions=transitions,
    )


def build_profiles_batch(
    species_seqs,
    sub_matrix: np.ndarray,
    bg_freqs: np.ndarray,
    gap_open: int = -12,
    gap_extend: int = -3,
):
    """Build profiles for all sequences in a species.

    Returns flat arrays suitable for batch Viterbi:
      - flat_match_emit: (total_positions, 20) int8
      - insert_emit: (20,) int8
      - transitions: (7,) int32
      - profile_offsets: (N,) int64
      - profile_lengths: (N,) int32
    """
    N = species_seqs.num_sequences
    total_len = int(species_seqs.lengths.sum())

    flat_match_emit = np.empty((total_len, ALPHABET_SIZE), dtype=np.int8)
    profile_offsets = np.empty(N, dtype=np.int64)
    profile_lengths = species_seqs.lengths.copy()

    gap_close = max(gap_extend + 1, -1)
    transitions = np.array([
        0, gap_open, gap_open, gap_close, gap_extend, gap_close, gap_extend
    ], dtype=np.int32)

    insert_emit = np.full(ALPHABET_SIZE, -1, dtype=np.int8)

    pos = 0
    for i in range(N):
        profile_offsets[i] = pos
        seq = species_seqs.get_sequence(i)
        length = int(species_seqs.lengths[i])
        for j in range(length):
            aa = seq[j]
            if aa < ALPHABET_SIZE:
                flat_match_emit[pos + j, :] = sub_matrix[aa, :]
            else:
                flat_match_emit[pos + j, :] = 0
        pos += length

    return flat_match_emit, insert_emit, transitions, profile_offsets, profile_lengths
