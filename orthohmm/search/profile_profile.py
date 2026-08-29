"""Calibrated profile-profile evidence for repairing split orthogroups."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
from numba import get_num_threads, njit, prange, set_num_threads

from ..accuracy import deduplicate_undirected_edges
from ..helpers import IndexedEdges
from .msa_profile import MSAProfile
from .prefilter import prefilter_candidates
from .sequences import SpeciesSequences


@dataclass(frozen=True)
class ProfilePairResult:
    """Profile-pair edges plus bounded-search diagnostics."""

    edges: IndexedEdges
    candidate_pairs: int
    reciprocal_pairs: int


def _pack_profile_distributions(
    profiles: Mapping[int, MSAProfile],
    background: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    profile_ids = np.asarray(sorted(profiles), dtype=np.int32)
    lengths = np.fromiter(
        (profiles[int(profile_id)].length for profile_id in profile_ids),
        dtype=np.int32,
        count=len(profile_ids),
    )
    offsets = np.zeros(len(profile_ids), dtype=np.int64)
    if len(profile_ids) > 1:
        offsets[1:] = np.cumsum(lengths[:-1], dtype=np.int64)
    if not len(profile_ids):
        empty = np.empty((0, 20), dtype=np.float32)
        return profile_ids, lengths, offsets, empty, empty.copy()
    emissions = np.vstack([
        profiles[int(profile_id)].match_emissions
        for profile_id in profile_ids
    ]).astype(np.float32)
    probabilities = background[None, :] * np.exp2(emissions / 2.0)
    probabilities /= probabilities.sum(axis=1, keepdims=True)
    return profile_ids, lengths, offsets, emissions, probabilities.astype(
        np.float32, copy=False
    )


def _candidate_pairs(
    profiles: Mapping[int, MSAProfile],
    profile_ids: np.ndarray,
    lengths: np.ndarray,
    offsets: np.ndarray,
    cpu: int,
    kmer_k: int,
    max_candidates_per_profile: int,
) -> np.ndarray:
    if not len(profile_ids):
        return np.empty((0, 2), dtype=np.int32)
    consensus = np.concatenate([
        profiles[int(profile_id)].consensus for profile_id in profile_ids
    ]).astype(np.uint8, copy=False)
    database = SpeciesSequences(
        "profile_consensus",
        [str(profile_id) for profile_id in profile_ids],
        consensus,
        offsets,
        lengths,
    )
    candidates, candidate_offsets = prefilter_candidates(
        database,
        database,
        k=kmer_k,
        use_reduced_alphabet=False,
        min_total_hits=2,
        min_diag_hits=1,
        max_candidates_per_query=max_candidates_per_profile,
        n_threads=cpu,
    )
    sources = np.repeat(
        np.arange(len(profile_ids), dtype=np.int32),
        np.diff(candidate_offsets),
    )
    keep = sources != candidates
    if not keep.any():
        return np.empty((0, 2), dtype=np.int32)
    left = np.minimum(sources[keep], candidates[keep])
    right = np.maximum(sources[keep], candidates[keep])
    encoded = np.unique(
        left.astype(np.int64) * len(profile_ids) + right.astype(np.int64)
    )
    return np.column_stack((
        encoded // len(profile_ids), encoded % len(profile_ids)
    )).astype(np.int32, copy=False)


@njit(cache=True, parallel=True)
def _batch_local_profile_scores(
    emissions: np.ndarray,
    probabilities: np.ndarray,
    offsets: np.ndarray,
    lengths: np.ndarray,
    pairs: np.ndarray,
    band_width: int,
    gap_penalty: float,
) -> np.ndarray:
    scores = np.zeros(len(pairs), dtype=np.float32)
    for pair_index in prange(len(pairs)):
        profile_a = pairs[pair_index, 0]
        profile_b = pairs[pair_index, 1]
        length_a = lengths[profile_a]
        length_b = lengths[profile_b]
        offset_a = offsets[profile_a]
        offset_b = offsets[profile_b]
        previous = np.zeros(length_b + 1, dtype=np.float32)
        best = 0.0
        for i in range(1, length_a + 1):
            current = np.zeros(length_b + 1, dtype=np.float32)
            lower = max(1, i - band_width)
            upper = min(length_b, i + band_width)
            for j in range(lower, upper + 1):
                row_a = offset_a + i - 1
                row_b = offset_b + j - 1
                match_score = 0.0
                for amino_acid in range(20):
                    match_score += 0.5 * (
                        probabilities[row_a, amino_acid]
                        * emissions[row_b, amino_acid]
                        + probabilities[row_b, amino_acid]
                        * emissions[row_a, amino_acid]
                    )
                value = max(
                    0.0,
                    previous[j - 1] + match_score,
                    previous[j] + gap_penalty,
                    current[j - 1] + gap_penalty,
                )
                current[j] = value
                if value > best:
                    best = value
            previous = current
        scores[pair_index] = best
    return scores


def profile_profile_similarities(
    profiles: Mapping[int, MSAProfile],
    background: np.ndarray,
    cpu: int,
    band_width: int = 64,
    gap_penalty: float = -5.0,
    kmer_k: int = 4,
    max_candidates_per_profile: int = 30,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return profile IDs, candidate pairs, and self-normalized similarities."""
    profile_ids, lengths, offsets, emissions, probabilities = (
        _pack_profile_distributions(profiles, background)
    )
    pairs = _candidate_pairs(
        profiles,
        profile_ids,
        lengths,
        offsets,
        cpu,
        kmer_k,
        max_candidates_per_profile,
    )
    if not len(pairs):
        return profile_ids, pairs, np.empty(0, dtype=np.float32)
    self_pairs = np.column_stack((
        np.arange(len(profile_ids), dtype=np.int32),
        np.arange(len(profile_ids), dtype=np.int32),
    ))
    previous_threads = get_num_threads()
    set_num_threads(max(1, min(int(cpu), previous_threads)))
    try:
        self_scores = _batch_local_profile_scores(
            emissions,
            probabilities,
            offsets,
            lengths,
            self_pairs,
            band_width,
            gap_penalty,
        )
        pair_scores = _batch_local_profile_scores(
            emissions,
            probabilities,
            offsets,
            lengths,
            pairs,
            band_width,
            gap_penalty,
        )
    finally:
        set_num_threads(previous_threads)
    denominators = np.sqrt(self_scores[pairs[:, 0]] * self_scores[pairs[:, 1]])
    similarities = np.divide(
        pair_scores,
        denominators,
        out=np.zeros_like(pair_scores),
        where=denominators > 0.0,
    )
    return profile_ids, pairs, similarities


def build_profile_profile_edges(
    profiles: Mapping[int, MSAProfile],
    clusters: Sequence[Sequence[int]],
    gene_names: Sequence[str],
    background: np.ndarray,
    cpu: int,
    similarity_threshold: float = 0.4,
    max_combined_genes: int = 80,
) -> ProfilePairResult:
    """Connect mutually closest profiles that clear a calibrated similarity."""
    if not 0.0 < similarity_threshold <= 1.0:
        raise ValueError("similarity_threshold must be in (0, 1]")
    if max_combined_genes < 2:
        raise ValueError("max_combined_genes must be >= 2")
    profile_ids, pairs, similarities = profile_profile_similarities(
        profiles, background, cpu
    )
    if not len(pairs):
        empty = deduplicate_undirected_edges(gene_names, [], [], [])
        return ProfilePairResult(empty, 0, 0)

    best_scores = np.full(len(profile_ids), -np.inf, dtype=np.float32)
    best_partners = np.full(len(profile_ids), -1, dtype=np.int32)
    for pair, similarity in sorted(
        zip(pairs, similarities),
        key=lambda item: (-float(item[1]), int(item[0][0]), int(item[0][1])),
    ):
        left, right = int(pair[0]), int(pair[1])
        if similarity > best_scores[left]:
            best_scores[left] = similarity
            best_partners[left] = right
        if similarity > best_scores[right]:
            best_scores[right] = similarity
            best_partners[right] = left

    selected = []
    for local_profile, partner in enumerate(best_partners):
        partner = int(partner)
        if partner <= local_profile or partner < 0:
            continue
        if best_partners[partner] != local_profile:
            continue
        similarity = float(best_scores[local_profile])
        if similarity < similarity_threshold:
            continue
        left_cluster = int(profile_ids[local_profile])
        right_cluster = int(profile_ids[partner])
        if len(clusters[left_cluster]) + len(clusters[right_cluster]) > max_combined_genes:
            continue
        selected.append((left_cluster, right_cluster, similarity))

    edge_sources = []
    edge_targets = []
    edge_weights = []
    for left_cluster, right_cluster, similarity in selected:
        for source in clusters[left_cluster]:
            for target in clusters[right_cluster]:
                edge_sources.append(int(source))
                edge_targets.append(int(target))
                edge_weights.append(similarity)
    edges = deduplicate_undirected_edges(
        gene_names, edge_sources, edge_targets, edge_weights
    )
    return ProfilePairResult(edges, len(pairs), len(selected))
