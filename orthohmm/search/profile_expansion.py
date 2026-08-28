"""Strict MSA-profile expansion for high-sensitivity orthogroup inference."""

from __future__ import annotations

import multiprocessing
import os
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from ..accuracy import deduplicate_undirected_edges
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
    """Significant gene-to-MSA-profile scores using global integer IDs."""

    profile_cluster_ids: np.ndarray
    gene_ids: np.ndarray
    scores: np.ndarray
    evalues: np.ndarray
    candidate_count: int


@dataclass(frozen=True)
class ProfileExpansionResult:
    """Strict profile edges and counts used by benchmark provenance."""

    edges: IndexedEdges
    profiles_built: int
    profile_candidates: int
    significant_profile_hits: int


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


def build_cluster_profiles(
    clusters: Sequence[Sequence[int]],
    gene_names: Sequence[str],
    sequence_database: SpeciesSequences,
    matrix_name: str,
    cpu: int,
    min_cluster_size: int = 3,
    max_cluster_size: int = 200,
) -> dict[int, MSAProfile]:
    """Build center-star MSA profiles for informative cluster sizes."""
    sub_matrix = get_matrix(matrix_name)
    background = get_background_freqs(matrix_name)
    cluster_ids = [
        cluster_id
        for cluster_id, cluster in enumerate(clusters)
        if min_cluster_size <= len(cluster) <= max_cluster_size
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

    profiles = {}
    workers = max(1, min(int(cpu), len(cluster_ids))) if cluster_ids else 1
    if workers == 1:
        for task in tasks():
            cluster_id, profile = _build_profile_worker(task)
            if profile is not None:
                profiles[cluster_id] = profile
        return profiles

    task_iter = iter(tasks())
    max_pending = workers * 2
    context = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(
        max_workers=workers,
        mp_context=context,
    ) as executor:
        pending = {
            executor.submit(_build_profile_worker, next(task_iter))
            for _ in range(min(max_pending, len(cluster_ids)))
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
    return ProfileHits(
        profile_ids[pairs[significant, 0]],
        pairs[significant, 1].copy(),
        scores[significant].astype(np.float64),
        evalues[significant],
        candidate_count,
    )


def compute_profile_self_thresholds(
    profiles: Mapping[int, MSAProfile],
    gene_names: Sequence[str],
    sequence_database: SpeciesSequences,
    cpu: int,
    band_width: int = 64,
) -> dict[int, float]:
    """Return each profile's weakest existing-member Viterbi score."""
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
    thresholds = {}
    for cluster_id, score in zip(pair_clusters, scores):
        thresholds[cluster_id] = min(
            thresholds.get(cluster_id, np.inf), float(score)
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
    thresholds = np.fromiter(
        (self_thresholds.get(int(profile_id), np.inf) for profile_id in profile_ids),
        dtype=np.float64,
        count=len(profile_ids),
    )
    eligible = (
        (cluster_by_gene[candidate_genes] != profile_ids)
        & (profile_scores >= thresholds)
    )
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
) -> ProfileExpansionResult:
    """Run the complete strict profile-expansion stage."""
    os.environ["OMP_NUM_THREADS"] = str(cpu)
    sequence_database = load_global_sequence_database(
        fasta_directory, files, gene_names
    )
    profiles = build_cluster_profiles(
        clusters,
        gene_names,
        sequence_database,
        matrix_name,
        cpu,
    )
    profile_hits = search_genes_against_profiles(
        profiles,
        sequence_database,
        matrix_name,
        cpu,
        evalue_threshold=evalue_threshold,
    )
    self_thresholds = compute_profile_self_thresholds(
        profiles,
        gene_names,
        sequence_database,
        cpu,
    )
    edges = build_strict_profile_edges(
        gene_names,
        clusters,
        profile_hits,
        self_thresholds,
        hit_queries,
        hit_targets,
        hit_scores,
    )
    return ProfileExpansionResult(
        edges=edges,
        profiles_built=len(profiles),
        profile_candidates=profile_hits.candidate_count,
        significant_profile_hits=len(profile_hits.scores),
    )
