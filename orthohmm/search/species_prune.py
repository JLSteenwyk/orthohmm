"""Phylogeny-aware species pruning.

Precompute a cheap k-mer signature per species, use it to estimate
pairwise species similarity, then skip the bottom-ranked target species
for each query species.

The default strategy (`top_k=None`) computes signatures but does no
pruning; pass `top_k=8` (for 12 species) to keep each query species's
8 most similar target species.

Signature: presence/absence of the top-4096 most frequent reduced-alphabet
5-mers. Stored as a dense uint8 array (4096 bytes per species). Pairwise
distance = 1 - Jaccard of present k-mers.
"""

from __future__ import annotations

import numpy as np
from typing import Dict, List, Tuple
from numba import njit, int64, int32, uint8

from .prefilter import REDUCED_ALPHA, REDUCED_ALPHA_SIZE, _remap_flat


SIG_K = 5
SIG_BITS = 4096  # 4 KB per species signature


@njit(cache=True)
def _top_kmers(flat_seq, offsets, lengths, k, alpha_size, top_n):
    """Return top-N most frequent k-mers in a species.

    Returns (kmer_codes_sorted, counts_sorted) both int64, sorted desc
    by count.
    """
    # Max kmer code value = alpha_size**k; for 10^5 = 1e5 kmers
    max_kmers = int64(1)
    for _ in range(k):
        max_kmers *= alpha_size
    counts = np.zeros(max_kmers, dtype=np.int64)

    n_seqs = len(offsets) - 1  # if offsets are [start0, start1, ..., end]
    # But our offsets layout has N entries, not N+1. Use lengths.
    n_seqs = len(lengths)

    for s in range(n_seqs):
        off = offsets[s]
        L = lengths[s]
        if L < k:
            continue
        # Rolling hash
        code = int64(0)
        mult = int64(1)
        for _ in range(k - 1):
            mult *= alpha_size
        # prime first window
        for j in range(k):
            code = code * alpha_size + flat_seq[off + j]
        counts[code] += 1
        for j in range(k, L):
            old_letter = flat_seq[off + j - k]
            code -= int64(old_letter) * mult
            code = code * alpha_size + flat_seq[off + j]
            counts[code] += 1

    # Argpartition to get top N
    if top_n > max_kmers:
        top_n = max_kmers
    idx = np.argpartition(-counts, top_n - 1)[:top_n]
    # sort by count desc
    sorted_idx = idx[np.argsort(-counts[idx])]
    return sorted_idx, counts[sorted_idx]


def build_species_signature(species_seqs) -> np.ndarray:
    """Build a 4096-byte presence signature for one species.

    Returns uint8 array of length SIG_BITS where entry[i] == 1 if one of
    the top-SIG_BITS reduced-alphabet 5-mers is present in this species.
    """
    reduced = _remap_flat(
        species_seqs.flat_sequences, REDUCED_ALPHA,
        len(species_seqs.flat_sequences),
    )
    top_kmers, _ = _top_kmers(
        reduced, species_seqs.offsets, species_seqs.lengths,
        int32(SIG_K), int32(REDUCED_ALPHA_SIZE), int64(SIG_BITS),
    )
    sig = np.zeros(SIG_BITS, dtype=np.uint8)
    sig[: len(top_kmers)] = 1  # all present by definition of "top"
    # But we want signatures to be comparable across species. Store the
    # actual top-kmer-codes; two species with overlapping top-K kmers
    # indicate similarity.
    return top_kmers  # int64 array, length SIG_BITS


def species_similarity_matrix(signatures: Dict[str, np.ndarray]
                              ) -> Tuple[List[str], np.ndarray]:
    """Pairwise Jaccard similarity of species signatures.

    Returns (species_list, sim_matrix) where sim_matrix[i,j] is |A∩B|/|A∪B|
    of the top-k k-mer sets.
    """
    names = sorted(signatures.keys())
    n = len(names)
    # Convert each signature to a set for fast intersection
    sets = {name: set(signatures[name].tolist()) for name in names}

    M = np.ones((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            a, b = sets[names[i]], sets[names[j]]
            inter = len(a & b)
            union = len(a | b)
            sim = inter / union if union > 0 else 0.0
            M[i, j] = sim
            M[j, i] = sim
    return names, M


def top_k_target_species(sim_matrix: np.ndarray, names: List[str],
                         top_k: int) -> Dict[str, List[str]]:
    """For each query species, return the top-K most similar target species.

    The query species itself is always included. top_k counts OTHER species
    (so top_k=8 for a query gives self + 8 others = 9 total targets).
    """
    result: Dict[str, List[str]] = {}
    n = len(names)
    for i, name in enumerate(names):
        # rank others by similarity desc
        others = [(sim_matrix[i, j], names[j]) for j in range(n) if j != i]
        others.sort(reverse=True)
        keep = [name]  # self first
        keep.extend([o[1] for o in others[:top_k]])
        result[name] = keep
    return result


def prune_species_pairs(all_pairs: List[Tuple[str, str]],
                        sim_matrix: np.ndarray,
                        names: List[str],
                        top_k: int) -> List[Tuple[str, str]]:
    """Given N² species pairs, return the subset using top-K target pruning.

    For 12 species with top_k=8 this returns ~12 * 9 = 108 pairs (out of 144).
    """
    keep_map = top_k_target_species(sim_matrix, names, top_k)
    kept = []
    for qf, tf in all_pairs:
        if tf in keep_map.get(qf, [qf]):
            kept.append((qf, tf))
    return kept
