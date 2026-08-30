"""Strict MSA-profile expansion for high-sensitivity orthogroup inference."""

from __future__ import annotations

import multiprocessing
import os
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import dataclass
from typing import Collection, Mapping, Sequence

import numpy as np

from ..accuracy import combine_edges, deduplicate_undirected_edges
from ..helpers import IndexedEdges
from .evalue import batch_estimate_evalues
from .matrices import (
    ALPHABET,
    get_background_freqs,
    get_ka_params,
    get_matrix,
)
from .msa_profile import MSAProfile, build_msa_profile
from .prefilter import prefilter_candidates
from .profile_profile import build_profile_profile_edges
from .sequences import SequenceStore, SpeciesSequences
from .viterbi import (
    batch_viterbi_c,
    batch_viterbi_multipair_c,
    batch_viterbi_score,
    is_c_available,
)

try:
    from .viterbi_cuda_ctypes import is_cuda_available, batch_viterbi_cuda
except Exception:  # pragma: no cover - CUDA is optional
    def is_cuda_available():
        return False

    batch_viterbi_cuda = None


@dataclass(frozen=True)
class ProfileHits:
    """Significant gene-to-MSA-profile scores using global integer IDs.

    Raw scores remain the source for E-values, profile ranking, and graph edge
    weights.  Per-target-residue scores are an optional decision statistic for
    comparing candidates with member-derived thresholds.
    """

    profile_cluster_ids: np.ndarray
    gene_ids: np.ndarray
    scores: np.ndarray
    evalues: np.ndarray
    candidate_count: int
    scores_per_target_residue: np.ndarray | None = None


@dataclass(frozen=True)
class ProfileExpansionResult:
    """Strict profile edges and counts used by benchmark provenance."""

    edges: IndexedEdges
    profiles_built: int
    profile_candidates: int
    significant_profile_hits: int
    strict_profile_edges: int = 0
    reciprocal_profile_pairs: int = 0
    reciprocal_profile_edges: int = 0
    calibrated_profiles: int = 0
    profile_pair_candidates: int = 0
    profile_pair_merges: int = 0
    profile_pair_edges: int = 0
    direct_profile_fallback_edges: int = 0
    species_completion_profiles: int = 0


def load_global_sequence_database(
    fasta_directory: str,
    files: Sequence[str],
    gene_names: Sequence[str],
) -> SpeciesSequences:
    """Load input FASTAs into one sequence database aligned to the gene table."""
    store = SequenceStore.from_fasta_directory(fasta_directory, list(files))
    sequence_by_id = {}
    for filename in files:
        species = store.species[filename]
        for index, gene in enumerate(species.ids):
            gene = str(gene)
            if gene in sequence_by_id:
                raise ValueError(f"duplicate gene ID across FASTA files: {gene}")
            sequence_by_id[gene] = species.get_sequence(index)
    missing = [gene for gene in gene_names if gene not in sequence_by_id]
    if missing or len(sequence_by_id) != len(gene_names):
        raise ValueError("FASTA sequences do not match the global gene table")
    ids = list(gene_names)
    sequences = [sequence_by_id[gene] for gene in ids]
    lengths = np.fromiter(
        (len(sequence) for sequence in sequences),
        dtype=np.int32,
        count=len(sequences),
    )
    offsets = np.zeros(len(sequences), dtype=np.int64)
    if len(sequences) > 1:
        offsets[1:] = np.cumsum(lengths[:-1], dtype=np.int64)
    flat = (
        np.concatenate(sequences).astype(np.uint8, copy=False)
        if sequences
        else np.empty(0, dtype=np.uint8)
    )
    return SpeciesSequences("all_genes", ids, flat, offsets, lengths)


def _decode_sequence(sequence: np.ndarray) -> str:
    return "".join(
        ALPHABET[int(amino_acid)] if amino_acid < 20 else "X"
        for amino_acid in sequence
    )


def _build_profile_worker(task):
    cluster_id, member_ids, sequences, sub_matrix, background = task
    try:
        profile = build_msa_profile(
            [_decode_sequence(sequence) for sequence in sequences],
            member_ids,
            sub_matrix,
            background,
        )
    except Exception:
        profile = None
    return cluster_id, profile


def _execute_profile_build_tasks(tasks, task_count: int, cpu: int):
    """Build a bounded stream of profiles serially or with spawn workers."""
    task_iter = iter(tasks)
    profiles = {}
    workers = max(1, min(int(cpu), task_count)) if task_count else 1
    if workers == 1:
        for task in task_iter:
            cluster_id, profile = _build_profile_worker(task)
            if profile is not None:
                profiles[cluster_id] = profile
        return profiles

    max_pending = workers * 2
    context = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(
        max_workers=workers,
        mp_context=context,
    ) as executor:
        pending = {
            executor.submit(_build_profile_worker, next(task_iter))
            for _ in range(min(max_pending, task_count))
        }
        while pending:
            completed, pending = wait(
                pending, return_when=FIRST_COMPLETED
            )
            for future in completed:
                cluster_id, profile = future.result()
                if profile is not None:
                    profiles[cluster_id] = profile
                try:
                    task = next(task_iter)
                except StopIteration:
                    continue
                pending.add(executor.submit(_build_profile_worker, task))
    return profiles


def build_cluster_profiles(
    clusters: Sequence[Sequence[int]],
    gene_names: Sequence[str],
    sequence_database: SpeciesSequences,
    matrix_name: str,
    cpu: int,
    min_cluster_size: int = 3,
    max_cluster_size: int = 200,
    gene_to_species=None,
    min_species_count: int = 1,
    allowed_cluster_ids: Collection[int] | None = None,
) -> dict[int, MSAProfile]:
    """Build center-star MSA profiles for informative cluster sizes."""
    if min_species_count < 1:
        raise ValueError("min_species_count must be >= 1")
    if min_species_count > 1 and gene_to_species is None:
        raise ValueError(
            "gene_to_species is required when min_species_count > 1"
        )
    species = (
        np.asarray(gene_to_species, dtype=np.int32)
        if gene_to_species is not None
        else None
    )
    if species is not None and len(species) != len(gene_names):
        raise ValueError("gene_to_species must match gene_names")
    sub_matrix = get_matrix(matrix_name)
    background = get_background_freqs(matrix_name)
    allowed = (
        set(int(cluster_id) for cluster_id in allowed_cluster_ids)
        if allowed_cluster_ids is not None
        else None
    )
    cluster_ids = [
        cluster_id
        for cluster_id, cluster in enumerate(clusters)
        if min_cluster_size <= len(cluster) <= max_cluster_size
        and (allowed is None or cluster_id in allowed)
        and (
            min_species_count == 1
            or len(np.unique(species[np.asarray(cluster, dtype=np.int32)]))
            >= min_species_count
        )
    ]

    def tasks():
        for cluster_id in cluster_ids:
            cluster = clusters[cluster_id]
            member_ids = [gene_names[int(gene)] for gene in cluster]
            sequences = [
                sequence_database.get_sequence(int(gene)) for gene in cluster
            ]
            yield (
                cluster_id, member_ids, sequences, sub_matrix, background
            )

    return _execute_profile_build_tasks(tasks(), len(cluster_ids), cpu)


def select_species_completion_cluster_ids(
    clusters: Sequence[Sequence[int]],
    gene_to_species,
    min_cluster_size: int = 3,
    max_cluster_size: int = 200,
) -> set[int]:
    """Select incomplete seed clusters containing at most one gene per species."""
    species = np.asarray(gene_to_species, dtype=np.int32)
    total_species = len(np.unique(species))
    selected = set()
    for cluster_id, cluster in enumerate(clusters):
        if not min_cluster_size <= len(cluster) <= max_cluster_size:
            continue
        represented = np.unique(species[np.asarray(cluster, dtype=np.int32)])
        if len(represented) == len(cluster) and len(represented) < total_species:
            selected.add(cluster_id)
    return selected


def select_single_copy_profile_ids(
    profiles: Mapping[int, MSAProfile],
    clusters: Sequence[Sequence[int]],
    gene_to_species,
) -> set[int]:
    """Return profiles whose source cluster has at most one gene per species."""
    if gene_to_species is None:
        raise ValueError("gene_to_species is required for single-copy profiles")
    species = np.asarray(gene_to_species, dtype=np.int32)
    selected = set()
    for profile_id in profiles:
        cluster = np.asarray(clusters[int(profile_id)], dtype=np.int32)
        if len(cluster) and len(np.unique(species[cluster])) == len(cluster):
            selected.add(int(profile_id))
    return selected


def _pack_profiles(profiles: Mapping[int, MSAProfile]):
    profile_ids = np.asarray(sorted(profiles), dtype=np.int32)
    if len(profile_ids) == 0:
        return profile_ids, None
    lengths = np.fromiter(
        (profiles[int(profile_id)].length for profile_id in profile_ids),
        dtype=np.int32,
        count=len(profile_ids),
    )
    offsets = np.zeros(len(profile_ids), dtype=np.int64)
    if len(profile_ids) > 1:
        offsets[1:] = np.cumsum(lengths[:-1], dtype=np.int64)
    flat_emissions = np.vstack([
        profiles[int(profile_id)].match_emissions
        for profile_id in profile_ids
    ]).astype(np.int8, copy=False)
    consensus = np.concatenate([
        profiles[int(profile_id)].consensus for profile_id in profile_ids
    ]).astype(np.uint8, copy=False)
    first = profiles[int(profile_ids[0])]
    packed = (
        lengths,
        offsets,
        flat_emissions,
        consensus,
        first.insert_emissions,
        first.transitions,
    )
    return profile_ids, packed


def _score_profile_pairs(
    flat_emissions,
    insert_emissions,
    transitions,
    profile_offsets,
    profile_lengths,
    sequence_database,
    pairs,
    band_width,
    cpu,
):
    """Route profile pairs through CUDA, AVX2, scalar C, or Numba."""
    def score_cpu(cpu_pairs):
        if not is_c_available():
            return batch_viterbi_score(
                flat_emissions,
                insert_emissions,
                transitions,
                profile_offsets,
                profile_lengths,
                sequence_database.flat_sequences,
                sequence_database.offsets,
                sequence_database.lengths,
                cpu_pairs,
                band_width,
            )
        try:
            return batch_viterbi_multipair_c(
                flat_emissions,
                insert_emissions,
                transitions,
                profile_offsets,
                profile_lengths,
                sequence_database.flat_sequences,
                sequence_database.offsets,
                sequence_database.lengths,
                cpu_pairs,
                band_width,
                n_threads=cpu,
            )
        except RuntimeError:
            return batch_viterbi_c(
                flat_emissions,
                insert_emissions,
                transitions,
                profile_offsets,
                profile_lengths,
                sequence_database.flat_sequences,
                sequence_database.offsets,
                sequence_database.lengths,
                cpu_pairs,
                band_width,
                n_threads=cpu,
            )

    scores = np.empty(len(pairs), dtype=np.int32)
    gpu_mask = np.zeros(len(pairs), dtype=bool)
    if is_cuda_available():
        gpu_mask = sequence_database.lengths[pairs[:, 1]] <= 1998
    if gpu_mask.any():
        try:
            scores[gpu_mask] = batch_viterbi_cuda(
                flat_emissions,
                insert_emissions,
                transitions,
                profile_offsets,
                profile_lengths,
                sequence_database.flat_sequences,
                sequence_database.offsets,
                sequence_database.lengths,
                pairs[gpu_mask],
                band_width,
            )
        except RuntimeError:
            gpu_mask.fill(False)
    if (~gpu_mask).any():
        scores[~gpu_mask] = score_cpu(pairs[~gpu_mask])
    return scores


def search_genes_against_profiles(
    profiles: Mapping[int, MSAProfile],
    sequence_database: SpeciesSequences,
    matrix_name: str,
    cpu: int,
    band_width: int = 64,
    kmer_k: int = 4,
    max_candidates_per_gene: int = 20,
    evalue_threshold: float = 1e-4,
) -> ProfileHits:
    """Prefilter and Viterbi-score all genes against MSA profiles."""
    profile_ids, packed = _pack_profiles(profiles)
    if packed is None:
        empty_i = np.empty(0, dtype=np.int32)
        empty_f = np.empty(0, dtype=np.float64)
        return ProfileHits(empty_i, empty_i.copy(), empty_f, empty_f.copy(), 0)
    (
        profile_lengths,
        profile_offsets,
        flat_emissions,
        consensus,
        insert_emissions,
        transitions,
    ) = packed
    profile_database = SpeciesSequences(
        "msa_profiles",
        [str(profile_id) for profile_id in profile_ids],
        consensus,
        profile_offsets,
        profile_lengths,
    )
    sub_matrix = get_matrix(matrix_name)
    candidate_profiles, candidate_offsets = prefilter_candidates(
        sequence_database,
        profile_database,
        k=kmer_k,
        use_reduced_alphabet=False,
        min_total_hits=4,
        min_diag_hits=1,
        max_candidates_per_query=max_candidates_per_gene,
        n_threads=cpu,
        sub_matrix=sub_matrix,
    )
    candidate_count = len(candidate_profiles)
    if candidate_count == 0:
        empty_i = np.empty(0, dtype=np.int32)
        empty_f = np.empty(0, dtype=np.float64)
        return ProfileHits(
            empty_i, empty_i.copy(), empty_f, empty_f.copy(), 0
        )

    gene_ids = np.repeat(
        np.arange(sequence_database.num_sequences, dtype=np.int32),
        np.diff(candidate_offsets),
    )
    pairs = np.column_stack((candidate_profiles, gene_ids)).astype(
        np.int32, copy=False
    )
    scores = _score_profile_pairs(
        flat_emissions,
        insert_emissions,
        transitions,
        profile_offsets,
        profile_lengths,
        sequence_database,
        pairs,
        band_width,
        cpu,
    )
    lam, k_parameter = get_ka_params(matrix_name)
    query_lengths = profile_lengths[pairs[:, 0]]
    evalues = batch_estimate_evalues(
        scores,
        query_lengths,
        int(sequence_database.lengths.sum()),
        lam,
        k_parameter,
    )
    significant = evalues < evalue_threshold
    significant_gene_ids = pairs[significant, 1].copy()
    scores_per_target_residue = scores[significant].astype(np.float64)
    scores_per_target_residue /= np.maximum(
        sequence_database.lengths[significant_gene_ids], 1
    )
    return ProfileHits(
        profile_ids[pairs[significant, 0]],
        significant_gene_ids,
        scores[significant].astype(np.float64),
        evalues[significant],
        candidate_count,
        scores_per_target_residue,
    )


def compute_profile_self_thresholds(
    profiles: Mapping[int, MSAProfile],
    gene_names: Sequence[str],
    sequence_database: SpeciesSequences,
    cpu: int,
    band_width: int = 64,
    matrix_name: str | None = None,
    calibrate_weakest_member: bool = False,
    calibration_profile_ids: Collection[int] | None = None,
    score_per_target_residue: bool = False,
    pair_profile_threshold_ratio: float = 1.0,
) -> dict[int, float]:
    """Return member thresholds, optionally calibrated by one jackknife.

    The default reproduces the historical in-sample minimum.  Calibration
    rebuilds each profile without its weakest in-sample member, scores that
    member against the held-out profile, and keeps the lower score.  This
    removes profile-construction optimism without using benchmark labels.
    Per-target-residue mode applies the same normalization to members and
    candidates so partial proteins are not compared with full-length totals.
    """
    if calibrate_weakest_member and not matrix_name:
        raise ValueError(
            "matrix_name is required for weakest-member calibration"
        )
    if not 0.0 < pair_profile_threshold_ratio <= 1.0:
        raise ValueError("pair_profile_threshold_ratio must be in (0, 1]")
    profile_ids, packed = _pack_profiles(profiles)
    if packed is None:
        return {}
    (
        profile_lengths,
        profile_offsets,
        flat_emissions,
        _consensus,
        insert_emissions,
        transitions,
    ) = packed
    gene_to_id = {gene: idx for idx, gene in enumerate(gene_names)}
    local_profile_ids = {int(profile_id): idx for idx, profile_id in enumerate(profile_ids)}
    pairs = []
    pair_clusters = []
    for profile_id in profile_ids:
        cluster_id = int(profile_id)
        for member in profiles[cluster_id].member_ids:
            gene_id = gene_to_id.get(member)
            if gene_id is not None:
                pairs.append((local_profile_ids[cluster_id], gene_id))
                pair_clusters.append(cluster_id)
    if not pairs:
        return {}
    pairs = np.asarray(pairs, dtype=np.int32)
    scores = _score_profile_pairs(
        flat_emissions,
        insert_emissions,
        transitions,
        profile_offsets,
        profile_lengths,
        sequence_database,
        pairs,
        band_width,
        cpu,
    )
    decision_scores = scores.astype(np.float64)
    if score_per_target_residue:
        decision_scores /= np.maximum(
            sequence_database.lengths[pairs[:, 1]], 1
        )
    thresholds = {}
    weakest_members = {}
    for cluster_id, pair, score in zip(
        pair_clusters, pairs, decision_scores
    ):
        score = float(score)
        if score < thresholds.get(cluster_id, np.inf):
            thresholds[cluster_id] = score
            weakest_members[cluster_id] = int(pair[1])
    if not calibrate_weakest_member:
        return thresholds

    calibrated_ids = (
        set(int(profile_id) for profile_id in calibration_profile_ids)
        if calibration_profile_ids is not None
        else set(weakest_members)
    )
    calibrated_ids.intersection_update(weakest_members)
    if not calibrated_ids:
        return thresholds

    sub_matrix = get_matrix(matrix_name)
    background = get_background_freqs(matrix_name)

    pair_tasks = []
    pair_task_metadata = []
    pair_task_id = 0
    for cluster_id in sorted(calibrated_ids):
        member_ids = profiles[cluster_id].member_ids
        if len(member_ids) != 2:
            continue
        for held_index in range(2):
            retained_name = member_ids[1 - held_index]
            held_name = member_ids[held_index]
            retained_gene = gene_to_id.get(retained_name)
            held_gene = gene_to_id.get(held_name)
            if retained_gene is None or held_gene is None:
                continue
            retained_sequence = sequence_database.get_sequence(retained_gene)
            pair_tasks.append((
                pair_task_id,
                [retained_name, retained_name],
                [retained_sequence, retained_sequence],
                sub_matrix,
                background,
            ))
            pair_task_metadata.append((cluster_id, held_gene))
            pair_task_id += 1
    pair_profiles = _execute_profile_build_tasks(
        pair_tasks, len(pair_tasks), cpu
    )
    if pair_profiles:
        pair_ids, pair_packed = _pack_profiles(pair_profiles)
        if pair_packed is not None:
            (
                pair_lengths,
                pair_offsets,
                pair_emissions,
                _pair_consensus,
                pair_inserts,
                pair_transitions,
            ) = pair_packed
            pair_pairs = np.asarray([
                (local_id, pair_task_metadata[int(task_id)][1])
                for local_id, task_id in enumerate(pair_ids)
            ], dtype=np.int32)
            pair_scores = _score_profile_pairs(
                pair_emissions,
                pair_inserts,
                pair_transitions,
                pair_offsets,
                pair_lengths,
                sequence_database,
                pair_pairs,
                band_width,
                cpu,
            ).astype(np.float64)
            if score_per_target_residue:
                pair_scores /= np.maximum(
                    sequence_database.lengths[pair_pairs[:, 1]], 1
                )
            for task_id, score in zip(pair_ids, pair_scores):
                cluster_id = pair_task_metadata[int(task_id)][0]
                thresholds[cluster_id] = min(
                    thresholds[cluster_id],
                    float(score) * pair_profile_threshold_ratio,
                )

    def jackknife_tasks():
        for cluster_id in sorted(calibrated_ids):
            weakest_gene_id = weakest_members[cluster_id]
            weakest_name = gene_names[weakest_gene_id]
            retained_names = [
                member
                for member in profiles[cluster_id].member_ids
                if member != weakest_name
            ]
            if len(retained_names) < 2:
                continue
            retained_sequences = [
                sequence_database.get_sequence(gene_to_id[member])
                for member in retained_names
                if member in gene_to_id
            ]
            if len(retained_sequences) != len(retained_names):
                continue
            yield (
                cluster_id,
                retained_names,
                retained_sequences,
                sub_matrix,
                background,
            )

    jackknife_cluster_ids = [
        cluster_id
        for cluster_id, weakest_gene_id in weakest_members.items()
        if cluster_id in calibrated_ids
        if len(profiles[cluster_id].member_ids) >= 3
        and gene_names[weakest_gene_id] in profiles[cluster_id].member_ids
    ]
    jackknife_profiles = _execute_profile_build_tasks(
        jackknife_tasks(), len(jackknife_cluster_ids), cpu
    )
    jackknife_ids, jackknife_packed = _pack_profiles(jackknife_profiles)
    if jackknife_packed is None:
        return thresholds
    (
        jackknife_lengths,
        jackknife_offsets,
        jackknife_emissions,
        _jackknife_consensus,
        jackknife_inserts,
        jackknife_transitions,
    ) = jackknife_packed
    jackknife_pairs = np.asarray([
        (local_profile_id, weakest_members[int(cluster_id)])
        for local_profile_id, cluster_id in enumerate(jackknife_ids)
    ], dtype=np.int32)
    jackknife_scores = _score_profile_pairs(
        jackknife_emissions,
        jackknife_inserts,
        jackknife_transitions,
        jackknife_offsets,
        jackknife_lengths,
        sequence_database,
        jackknife_pairs,
        band_width,
        cpu,
    )
    jackknife_decision_scores = jackknife_scores.astype(np.float64)
    if score_per_target_residue:
        jackknife_decision_scores /= np.maximum(
            sequence_database.lengths[jackknife_pairs[:, 1]], 1
        )
    for cluster_id, score in zip(
        jackknife_ids, jackknife_decision_scores
    ):
        cluster_id = int(cluster_id)
        thresholds[cluster_id] = min(
            thresholds[cluster_id], float(score)
        )
    return thresholds


def build_strict_profile_edges(
    gene_names: Sequence[str],
    clusters: Sequence[Sequence[int]],
    profile_hits: ProfileHits,
    self_thresholds: Mapping[int, float],
    hit_queries,
    hit_targets,
    hit_scores,
    score_per_target_residue: bool = False,
    candidate_missing_species_only: bool = False,
    gene_to_species=None,
) -> IndexedEdges:
    """Create one sequence-supported anchor for each strict profile candidate."""
    n_genes = len(gene_names)
    cluster_by_gene = np.full(n_genes, -1, dtype=np.int32)
    member_rank = np.full(n_genes, np.iinfo(np.int32).max, dtype=np.int32)
    for cluster_id, cluster in enumerate(clusters):
        members = np.asarray(cluster, dtype=np.int32)
        cluster_by_gene[members] = cluster_id
        member_rank[members] = np.arange(len(members), dtype=np.int32)

    profile_ids = np.asarray(profile_hits.profile_cluster_ids, dtype=np.int32)
    candidate_genes = np.asarray(profile_hits.gene_ids, dtype=np.int32)
    profile_scores = np.asarray(profile_hits.scores, dtype=np.float64)
    if len(profile_ids) == 0:
        return deduplicate_undirected_edges(gene_names, [], [], [])
    decision_scores = profile_scores
    if score_per_target_residue:
        if profile_hits.scores_per_target_residue is None:
            raise ValueError(
                "per-residue profile scores are required for normalized "
                "thresholding"
            )
        decision_scores = np.asarray(
            profile_hits.scores_per_target_residue, dtype=np.float64
        )
    thresholds = np.fromiter(
        (self_thresholds.get(int(profile_id), np.inf) for profile_id in profile_ids),
        dtype=np.float64,
        count=len(profile_ids),
    )
    eligible = (
        (cluster_by_gene[candidate_genes] != profile_ids)
        & (decision_scores >= thresholds)
    )
    if candidate_missing_species_only:
        if gene_to_species is None:
            raise ValueError(
                "gene_to_species is required for missing-species candidates"
            )
        species = np.asarray(gene_to_species, dtype=np.int32)
        if len(species) != n_genes:
            raise ValueError("gene_to_species must match gene_names")
        species_count = int(species.max()) + 1 if len(species) else 0
        represented = np.zeros(
            (len(clusters), species_count), dtype=np.bool_
        )
        for cluster_id, cluster in enumerate(clusters):
            members = np.asarray(cluster, dtype=np.int32)
            represented[cluster_id, species[members]] = True
        eligible &= ~represented[
            profile_ids, species[candidate_genes]
        ]
    if not eligible.any():
        return deduplicate_undirected_edges(gene_names, [], [], [])

    eligible_rows = np.flatnonzero(eligible)
    best_scores = np.full(n_genes, -np.inf, dtype=np.float64)
    np.maximum.at(
        best_scores,
        candidate_genes[eligible_rows],
        profile_scores[eligible_rows],
    )
    winning_rows = eligible_rows[
        profile_scores[eligible_rows]
        == best_scores[candidate_genes[eligible_rows]]
    ]
    winning_genes, first = np.unique(
        candidate_genes[winning_rows], return_index=True
    )
    selected_profile = np.full(n_genes, -1, dtype=np.int32)
    selected_profile[winning_genes] = profile_ids[winning_rows[first]]

    queries = np.asarray(hit_queries, dtype=np.int32)
    targets = np.asarray(hit_targets, dtype=np.int32)
    scores = np.asarray(hit_scores, dtype=np.float64)
    outgoing = (
        (selected_profile[queries] >= 0)
        & (cluster_by_gene[targets] == selected_profile[queries])
        & (queries != targets)
        & (scores > 0.0)
    )
    reverse = (
        (selected_profile[targets] >= 0)
        & (cluster_by_gene[queries] == selected_profile[targets])
        & (queries != targets)
        & (scores > 0.0)
    )
    outgoing_keys = (
        queries[outgoing].astype(np.int64) * n_genes + targets[outgoing]
    )
    reverse_rows = np.flatnonzero(reverse)
    reverse_keys = (
        targets[reverse_rows].astype(np.int64) * n_genes
        + queries[reverse_rows]
    )
    if len(outgoing_keys) and len(reverse_keys):
        reverse_rows = reverse_rows[
            ~np.isin(reverse_keys, np.unique(outgoing_keys))
        ]

    anchor_genes = np.concatenate((queries[outgoing], targets[reverse_rows]))
    anchor_members = np.concatenate((targets[outgoing], queries[reverse_rows]))
    anchor_scores = np.concatenate((scores[outgoing], scores[reverse_rows]))
    if len(anchor_genes) == 0:
        return deduplicate_undirected_edges(gene_names, [], [], [])
    best_anchor_scores = np.full(n_genes, -np.inf, dtype=np.float64)
    np.maximum.at(best_anchor_scores, anchor_genes, anchor_scores)
    tied_anchor_rows = np.flatnonzero(
        anchor_scores == best_anchor_scores[anchor_genes]
    )
    # Historical strict expansion iterated cluster members in file order and
    # retained the first member on equal sequence-hit scores. Preserve that
    # deterministic tie rule instead of depending on hit-array order.
    order = np.lexsort((
        member_rank[anchor_members[tied_anchor_rows]],
        anchor_genes[tied_anchor_rows],
    ))
    ordered_rows = tied_anchor_rows[order]
    selected_genes, first = np.unique(
        anchor_genes[ordered_rows], return_index=True
    )
    selected_rows = ordered_rows[first]
    return deduplicate_undirected_edges(
        gene_names,
        selected_genes,
        anchor_members[selected_rows],
        anchor_scores[selected_rows],
    )


def build_direct_profile_fallback_edges(
    gene_names: Sequence[str],
    clusters: Sequence[Sequence[int]],
    profile_hits: ProfileHits,
    self_thresholds: Mapping[int, float],
    hit_queries,
    hit_targets,
    hit_scores,
    allowed_profile_ids: Collection[int],
    score_per_target_residue: bool = False,
) -> IndexedEdges:
    """Anchor calibrated profile winners lacking a pairwise sequence hit."""
    n_genes = len(gene_names)
    cluster_by_gene = np.full(n_genes, -1, dtype=np.int32)
    representatives = np.full(len(clusters), -1, dtype=np.int32)
    for cluster_id, cluster in enumerate(clusters):
        members = np.asarray(cluster, dtype=np.int32)
        if len(members):
            cluster_by_gene[members] = cluster_id
            representatives[cluster_id] = int(members[0])

    profile_ids = np.asarray(profile_hits.profile_cluster_ids, dtype=np.int32)
    candidate_genes = np.asarray(profile_hits.gene_ids, dtype=np.int32)
    raw_scores = np.asarray(profile_hits.scores, dtype=np.float64)
    if not len(profile_ids) or not allowed_profile_ids:
        return deduplicate_undirected_edges(gene_names, [], [], [])
    decision_scores = raw_scores
    if score_per_target_residue:
        if profile_hits.scores_per_target_residue is None:
            raise ValueError(
                "per-residue profile scores are required for normalized "
                "thresholding"
            )
        decision_scores = np.asarray(
            profile_hits.scores_per_target_residue, dtype=np.float64
        )
    thresholds = np.fromiter(
        (
            self_thresholds.get(int(profile_id), np.inf)
            for profile_id in profile_ids
        ),
        dtype=np.float64,
        count=len(profile_ids),
    )
    allowed = np.isin(
        profile_ids,
        np.fromiter(allowed_profile_ids, dtype=np.int32),
    )
    eligible = (
        allowed
        & (cluster_by_gene[candidate_genes] != profile_ids)
        & np.isfinite(decision_scores)
        & np.isfinite(thresholds)
        & (thresholds > 0.0)
        & (decision_scores >= thresholds)
    )
    if not eligible.any():
        return deduplicate_undirected_edges(gene_names, [], [], [])

    eligible_rows = np.flatnonzero(eligible)
    best_scores = np.full(n_genes, -np.inf, dtype=np.float64)
    np.maximum.at(
        best_scores,
        candidate_genes[eligible_rows],
        raw_scores[eligible_rows],
    )
    winning_rows = eligible_rows[
        raw_scores[eligible_rows]
        == best_scores[candidate_genes[eligible_rows]]
    ]
    winning_genes, first = np.unique(
        candidate_genes[winning_rows], return_index=True
    )
    selected_rows = winning_rows[first]
    selected_profiles = profile_ids[selected_rows]

    queries = np.asarray(hit_queries, dtype=np.int32)
    targets = np.asarray(hit_targets, dtype=np.int32)
    scores = np.asarray(hit_scores, dtype=np.float64)
    anchored = np.zeros(n_genes, dtype=np.bool_)
    selected_profile_by_gene = np.full(n_genes, -1, dtype=np.int32)
    selected_profile_by_gene[winning_genes] = selected_profiles
    outgoing = (
        (selected_profile_by_gene[queries] >= 0)
        & (
            cluster_by_gene[targets]
            == selected_profile_by_gene[queries]
        )
        & (queries != targets)
        & (scores > 0.0)
    )
    reverse = (
        (selected_profile_by_gene[targets] >= 0)
        & (
            cluster_by_gene[queries]
            == selected_profile_by_gene[targets]
        )
        & (queries != targets)
        & (scores > 0.0)
    )
    anchored[queries[outgoing]] = True
    anchored[targets[reverse]] = True

    fallback = ~anchored[winning_genes]
    fallback_genes = winning_genes[fallback]
    fallback_rows = selected_rows[fallback]
    fallback_profiles = profile_ids[fallback_rows]
    fallback_members = representatives[fallback_profiles]
    valid = fallback_members >= 0
    confidence = np.clip(
        decision_scores[fallback_rows] / thresholds[fallback_rows],
        1.0,
        5.0,
    )
    return deduplicate_undirected_edges(
        gene_names,
        fallback_genes[valid],
        fallback_members[valid],
        confidence[valid],
    )


def build_reciprocal_profile_edges(
    gene_names: Sequence[str],
    clusters: Sequence[Sequence[int]],
    profile_hits: ProfileHits,
    self_thresholds: Mapping[int, float],
    gene_to_species,
    threshold_ratio: float = 0.7,
    min_support: int = 2,
    score_per_target_residue: bool = False,
) -> tuple[IndexedEdges, int]:
    """Connect complementary split clusters with reciprocal profile evidence.

    A directed observation is a member of one cluster scoring against another
    cluster's profile.  Edges are emitted only when both directions have the
    requested number of unique supporting genes and the two clusters have no
    species in common.  Edge weights are the weaker of the two profile scores
    after calibration by each profile's weakest-member score.
    """
    if not 0.0 < threshold_ratio <= 1.0:
        raise ValueError("threshold_ratio must be in (0, 1]")
    if min_support < 1:
        raise ValueError("min_support must be >= 1")
    if gene_to_species is None:
        raise ValueError("gene_to_species is required for reciprocal profiles")
    species = np.asarray(gene_to_species, dtype=np.int32)
    if len(species) != len(gene_names):
        raise ValueError("gene_to_species must match gene_names")

    n_genes = len(gene_names)
    cluster_by_gene = np.full(n_genes, -1, dtype=np.int32)
    cluster_species = []
    for cluster_id, cluster in enumerate(clusters):
        members = np.asarray(cluster, dtype=np.int32)
        cluster_by_gene[members] = cluster_id
        cluster_species.append(set(int(item) for item in species[members]))

    profile_ids = np.asarray(profile_hits.profile_cluster_ids, dtype=np.int32)
    candidate_genes = np.asarray(profile_hits.gene_ids, dtype=np.int32)
    profile_scores = np.asarray(profile_hits.scores, dtype=np.float64)
    if len(profile_ids) == 0:
        return deduplicate_undirected_edges(gene_names, [], [], []), 0
    if score_per_target_residue:
        if profile_hits.scores_per_target_residue is None:
            raise ValueError(
                "per-residue profile scores are required for normalized "
                "thresholding"
            )
        profile_scores = np.asarray(
            profile_hits.scores_per_target_residue, dtype=np.float64
        )
    if not (
        len(profile_ids) == len(candidate_genes) == len(profile_scores)
    ):
        raise ValueError("profile-hit arrays must align")
    valid_profiles = (profile_ids >= 0) & (profile_ids < len(clusters))
    valid_genes = (candidate_genes >= 0) & (candidate_genes < n_genes)
    thresholds = np.fromiter(
        (
            self_thresholds.get(int(profile_id), np.inf)
            for profile_id in profile_ids
        ),
        dtype=np.float64,
        count=len(profile_ids),
    )
    source_clusters = np.full(len(candidate_genes), -1, dtype=np.int32)
    source_clusters[valid_genes] = cluster_by_gene[candidate_genes[valid_genes]]
    eligible = (
        valid_profiles
        & valid_genes
        & (source_clusters >= 0)
        & (source_clusters != profile_ids)
        & np.isfinite(thresholds)
        & (thresholds > 0.0)
        & np.isfinite(profile_scores)
        & (profile_scores >= thresholds * threshold_ratio)
    )
    if not eligible.any():
        return deduplicate_undirected_edges(gene_names, [], [], []), 0

    directed: dict[tuple[int, int], dict[int, float]] = {}
    for row in np.flatnonzero(eligible):
        source = int(source_clusters[row])
        target = int(profile_ids[row])
        gene = int(candidate_genes[row])
        ratio = float(profile_scores[row] / thresholds[row])
        support = directed.setdefault((source, target), {})
        support[gene] = max(support.get(gene, -np.inf), ratio)

    edge_sources = []
    edge_targets = []
    edge_weights = []
    reciprocal_pairs = 0
    for source, target in sorted(directed):
        if source >= target:
            continue
        forward = directed[(source, target)]
        reverse = directed.get((target, source))
        if reverse is None:
            continue
        if len(forward) < min_support or len(reverse) < min_support:
            continue
        if cluster_species[source] & cluster_species[target]:
            continue
        reciprocal_pairs += 1
        for source_gene, forward_ratio in sorted(forward.items()):
            for target_gene, reverse_ratio in sorted(reverse.items()):
                edge_sources.append(source_gene)
                edge_targets.append(target_gene)
                edge_weights.append(min(forward_ratio, reverse_ratio))

    return (
        deduplicate_undirected_edges(
            gene_names, edge_sources, edge_targets, edge_weights
        ),
        reciprocal_pairs,
    )


def expand_profiles(
    clusters: Sequence[Sequence[int]],
    gene_names: Sequence[str],
    fasta_directory: str,
    files: Sequence[str],
    matrix_name: str,
    cpu: int,
    evalue_threshold: float,
    hit_queries,
    hit_targets,
    hit_scores,
    calibrate_weakest_member: bool = False,
    gene_to_species=None,
    min_profile_species: int = 1,
    calibrate_single_copy_profiles: bool = False,
    reciprocal_profile_merges: bool = False,
    reciprocal_profile_threshold_ratio: float = 0.7,
    reciprocal_profile_min_support: int = 2,
    score_per_target_residue: bool = False,
    profile_profile_merges: bool = False,
    profile_profile_similarity_threshold: float = 0.6,
    profile_profile_max_combined_genes: int = 80,
    direct_profile_fallback: bool = False,
    profile_min_cluster_size: int = 3,
    pair_profile_threshold_ratio: float = 1.0,
    species_completion_only: bool = False,
) -> ProfileExpansionResult:
    """Run the complete strict profile-expansion stage."""
    os.environ["OMP_NUM_THREADS"] = str(cpu)
    sequence_database = load_global_sequence_database(
        fasta_directory, files, gene_names
    )
    completion_profile_ids = None
    if species_completion_only:
        if gene_to_species is None:
            raise ValueError(
                "gene_to_species is required for species completion"
            )
        completion_profile_ids = select_species_completion_cluster_ids(
            clusters,
            gene_to_species,
            min_cluster_size=profile_min_cluster_size,
        )
    profiles = build_cluster_profiles(
        clusters,
        gene_names,
        sequence_database,
        matrix_name,
        cpu,
        min_cluster_size=profile_min_cluster_size,
        gene_to_species=gene_to_species,
        min_species_count=min_profile_species,
        allowed_cluster_ids=completion_profile_ids,
    )
    profile_hits = search_genes_against_profiles(
        profiles,
        sequence_database,
        matrix_name,
        cpu,
        evalue_threshold=evalue_threshold,
    )
    calibration_profile_ids = None
    calibrate_profiles = calibrate_weakest_member
    if calibrate_single_copy_profiles and not calibrate_weakest_member:
        calibration_profile_ids = select_single_copy_profile_ids(
            profiles, clusters, gene_to_species
        )
        calibrate_profiles = True
    self_thresholds = compute_profile_self_thresholds(
        profiles,
        gene_names,
        sequence_database,
        cpu,
        matrix_name=matrix_name,
        calibrate_weakest_member=calibrate_profiles,
        calibration_profile_ids=calibration_profile_ids,
        score_per_target_residue=score_per_target_residue,
        pair_profile_threshold_ratio=pair_profile_threshold_ratio,
    )
    edges = build_strict_profile_edges(
        gene_names,
        clusters,
        profile_hits,
        self_thresholds,
        hit_queries,
        hit_targets,
        hit_scores,
        score_per_target_residue=score_per_target_residue,
        candidate_missing_species_only=species_completion_only,
        gene_to_species=gene_to_species,
    )
    strict_edge_count = len(edges)
    direct_fallback_edge_count = 0
    if direct_profile_fallback:
        fallback_edges = build_direct_profile_fallback_edges(
            gene_names,
            clusters,
            profile_hits,
            self_thresholds,
            hit_queries,
            hit_targets,
            hit_scores,
            calibration_profile_ids or set(),
            score_per_target_residue=score_per_target_residue,
        )
        direct_fallback_edge_count = len(fallback_edges)
        edges = combine_edges(edges, fallback_edges)
    reciprocal_pairs = 0
    reciprocal_edge_count = 0
    if reciprocal_profile_merges:
        reciprocal_edges, reciprocal_pairs = build_reciprocal_profile_edges(
            gene_names,
            clusters,
            profile_hits,
            self_thresholds,
            gene_to_species,
            threshold_ratio=reciprocal_profile_threshold_ratio,
            min_support=reciprocal_profile_min_support,
            score_per_target_residue=score_per_target_residue,
        )
        reciprocal_edge_count = len(reciprocal_edges)
        edges = combine_edges(edges, reciprocal_edges)
    profile_pair_candidates = 0
    profile_pair_merges = 0
    profile_pair_edge_count = 0
    if profile_profile_merges:
        profile_pair_result = build_profile_profile_edges(
            profiles,
            clusters,
            gene_names,
            get_background_freqs(matrix_name),
            cpu,
            similarity_threshold=profile_profile_similarity_threshold,
            max_combined_genes=profile_profile_max_combined_genes,
        )
        profile_pair_candidates = profile_pair_result.candidate_pairs
        profile_pair_merges = profile_pair_result.reciprocal_pairs
        profile_pair_edge_count = len(profile_pair_result.edges)
        edges = combine_edges(edges, profile_pair_result.edges)
    return ProfileExpansionResult(
        edges=edges,
        profiles_built=len(profiles),
        profile_candidates=profile_hits.candidate_count,
        significant_profile_hits=len(profile_hits.scores),
        strict_profile_edges=strict_edge_count,
        reciprocal_profile_pairs=reciprocal_pairs,
        reciprocal_profile_edges=reciprocal_edge_count,
        calibrated_profiles=(
            len(profiles)
            if calibrate_weakest_member
            else len(calibration_profile_ids or ())
        ),
        profile_pair_candidates=profile_pair_candidates,
        profile_pair_merges=profile_pair_merges,
        profile_pair_edges=profile_pair_edge_count,
        direct_profile_fallback_edges=direct_fallback_edge_count,
        species_completion_profiles=(
            len(profiles) if species_completion_only else 0
        ),
    )
