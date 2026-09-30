"""Single-sequence profile HMM construction.

Builds an integer scoring profile for local match/insert/delete Viterbi
alignment. This is not an exact implementation of Plan7 or phmmer scoring.

For position i with amino acid a_i:
  - Match emission scores = substitution_matrix[a_i, :]
  - Insert emission scores = uniform -1 (position-independent)
  - Transition costs = uniform additive integer penalties

The supplied substitution matrix rows are used directly as integer match
scores. The background-frequency argument does not alter these emissions.
"""

from dataclasses import dataclass

import numpy as np
from numba import njit, int32, int8

from .matrices import ALPHABET_SIZE


# ──────────────────────────────────────────────────────────────────────
# Default additive costs in MM, MI, MD, IM, II, DM, DD order:
# [0, -12, -12, -1, -3, -1, -3]. These are not normalized probabilities.
# Insertions also pay the emission cost of -1 per residue; deletions do not.
# Raw scores determine the E-value gate before downstream normalization,
# so changing these penalties changes search decisions.
# ──────────────────────────────────────────────────────────────────────

# Transition indices in the position-independent (7,) transitions array
T_MM = 0  # Match -> Match
T_MI = 1  # Match -> Insert
T_MD = 2  # Match -> Delete
T_IM = 3  # Insert -> Match
T_II = 4  # Insert -> Insert
T_DM = 5  # Delete -> Match
T_DD = 6  # Delete -> Delete


@dataclass
class ProfileHMM:
    """Integer match/insert/delete scoring profile for one query sequence."""
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
    """Build a single-sequence match/insert/delete scoring profile.

    Parameters
    ----------
    query_seq : uint8 array of length query_len
    query_len : int
    sub_matrix : (20, 20) int8 substitution matrix
    bg_freqs : (20,) float64 background frequencies, retained for API compatibility
        but not used to construct single-sequence emissions
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

    # Uniform insert penalty, independent of bg_freqs.
    insert_emissions = np.full(ALPHABET_SIZE, -1, dtype=np.int8)

    # Transitions: uniform across all positions
    # The default gap_extend=-3 gives gap_close=-1.
    gap_close = max(gap_extend + 1, -1)

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

    flat_match_emit = np.zeros((total_len, ALPHABET_SIZE), dtype=np.int8)
    profile_offsets = species_seqs.offsets.copy()
    profile_lengths = species_seqs.lengths.copy()

    gap_close = max(gap_extend + 1, -1)
    transitions = np.array([
        0, gap_open, gap_open, gap_close, gap_extend, gap_close, gap_extend
    ], dtype=np.int32)

    insert_emit = np.full(ALPHABET_SIZE, -1, dtype=np.int8)

    np.take(
        sub_matrix,
        species_seqs.flat_sequences,
        axis=0,
        out=flat_match_emit,
        mode="clip",
    )
    flat_match_emit[species_seqs.flat_sequences >= ALPHABET_SIZE] = 0

    return flat_match_emit, insert_emit, transitions, profile_offsets, profile_lengths
