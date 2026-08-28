"""Accuracy-oriented search and graph-inference profiles.

The high-sensitivity path keeps the standard pipeline unchanged while
productionizing the broader search and multipass graph recipe validated on
OrthoBench. All graph operations use global integer gene IDs so benchmark-scale
datasets do not require Python dictionaries for every hit or edge.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .helpers import IndexedEdges


DEFAULT_ACCURACY_PROFILE = "standard"


@dataclass(frozen=True)
class AccuracyProfile:
    """Parameters that change sequence-search sensitivity and graph inference."""

    name: str
    kmer_k: int
    max_candidates_per_query: int
    multipass_graph: bool
    profile_expansion: bool
    leiden_seed: int


_ACCURACY_PROFILES = {
    "standard": AccuracyProfile(
        name="standard",
        kmer_k=5,
        max_candidates_per_query=0,
        multipass_graph=False,
        profile_expansion=False,
        leiden_seed=0,
    ),
    "high_sensitivity": AccuracyProfile(
        name="high_sensitivity",
        kmer_k=4,
        max_candidates_per_query=100,
        multipass_graph=True,
        profile_expansion=True,
        leiden_seed=4,
    ),
}


def resolve_accuracy_profile(name: str) -> AccuracyProfile:
    """Resolve a named accuracy profile without exposing mutable state."""
    profile_name = (
        name.lower() if isinstance(name, str) else DEFAULT_ACCURACY_PROFILE
    )
    try:
        return _ACCURACY_PROFILES[profile_name]
    except KeyError as exc:
        choices = ", ".join(sorted(_ACCURACY_PROFILES))
        raise ValueError(
            f"Unknown accuracy profile {name!r}; expected one of: {choices}"
        ) from exc


def _empty_edges(gene_names: Sequence[str]) -> IndexedEdges:
    return IndexedEdges(
        list(gene_names),
        np.empty(0, dtype=np.int32),
        np.empty(0, dtype=np.int32),
        np.empty(0, dtype=np.float64),
    )


def deduplicate_undirected_edges(
    gene_names: Sequence[str],
    sources,
    targets,
    weights,
) -> IndexedEdges:
    """Canonicalize undirected edges and retain the maximum duplicate weight."""
    sources = np.asarray(sources, dtype=np.int32)
    targets = np.asarray(targets, dtype=np.int32)
    weights = np.asarray(weights, dtype=np.float64)
    if not (len(sources) == len(targets) == len(weights)):
        raise ValueError("edge source, target, and weight arrays must align")
    if len(sources) == 0:
        return _empty_edges(gene_names)

    low = np.minimum(sources, targets)
    high = np.maximum(sources, targets)
    keep = (low != high) & np.isfinite(weights)
    if not keep.any():
        return _empty_edges(gene_names)

    low = low[keep]
    high = high[keep]
    weights = weights[keep]
    packed = (low.astype(np.int64) << 32) | high.astype(np.int64)
    order = np.argsort(packed, kind="stable")
    packed = packed[order]
    weights = weights[order]
    starts = np.concatenate(
        (np.array([0], dtype=np.int64), np.flatnonzero(np.diff(packed)) + 1)
    )
    unique = packed[starts]
    return IndexedEdges(
        list(gene_names),
        (unique >> 32).astype(np.int32),
        (unique & 0xFFFFFFFF).astype(np.int32),
        np.maximum.reduceat(weights, starts),
    )


def combine_edges(*edge_sets: IndexedEdges) -> IndexedEdges:
    """Merge edge collections that share the same global gene table."""
    if not edge_sets:
        return _empty_edges([])
    gene_names = edge_sets[0].gene_names
    if any(edges.gene_names != gene_names for edges in edge_sets[1:]):
        raise ValueError("cannot merge edges with different gene tables")
    return deduplicate_undirected_edges(
        gene_names,
        np.concatenate([edges.sources for edges in edge_sets]),
        np.concatenate([edges.targets for edges in edge_sets]),
        np.concatenate([edges.weights for edges in edge_sets]),
    )


def build_rbnh_edges(
    gene_names: Sequence[str],
    gene_to_species,
    hit_queries,
    hit_targets,
    hit_scores,
) -> IndexedEdges:
    """Build reciprocal-best normalized-hit edges from significant hits.

    A gene's edge threshold is its weakest reciprocal best hit across target
    species. A directed hit becomes an undirected edge when it satisfies at
    least one endpoint's threshold, matching the validated multipass recipe.
    """
    species = np.asarray(gene_to_species, dtype=np.int32)
    queries = np.asarray(hit_queries, dtype=np.int32)
    targets = np.asarray(hit_targets, dtype=np.int32)
    scores = np.asarray(hit_scores, dtype=np.float64)
    n_genes = len(gene_names)
    if len(species) != n_genes:
        raise ValueError("gene species IDs must align with the gene table")
    if not (len(queries) == len(targets) == len(scores)):
        raise ValueError("hit query, target, and score arrays must align")
    if len(queries) == 0:
        return _empty_edges(gene_names)
    if (
        queries.min() < 0
        or targets.min() < 0
        or queries.max() >= n_genes
        or targets.max() >= n_genes
    ):
        raise ValueError("hit IDs fall outside the gene table")

    n_species = int(species.max()) + 1 if len(species) else 0
    keys = queries.astype(np.int64) * n_species + species[targets]
    n_slots = n_genes * n_species
    best_scores = np.full(n_slots, -np.inf, dtype=np.float64)
    np.maximum.at(best_scores, keys, scores)

    # Preserve the first target on exact ties, matching the historical strict
    # greater-than update used by the validated benchmark implementation.
    winning_rows = np.flatnonzero(scores == best_scores[keys])
    winning_keys, first = np.unique(keys[winning_rows], return_index=True)
    selected_rows = winning_rows[first]
    best_targets = np.full(n_slots, -1, dtype=np.int32)
    best_targets[winning_keys] = targets[selected_rows]

    source_genes = (winning_keys // n_species).astype(np.int32)
    target_genes = best_targets[winning_keys]
    reverse_keys = (
        target_genes.astype(np.int64) * n_species + species[source_genes]
    )
    reciprocal = best_targets[reverse_keys] == source_genes
    reciprocal_sources = source_genes[reciprocal]
    reciprocal_targets = target_genes[reciprocal]
    reciprocal_keys = winning_keys[reciprocal]
    reverse_keys = reverse_keys[reciprocal]
    reciprocal_scores = (
        best_scores[reciprocal_keys] + best_scores[reverse_keys]
    ) / 2.0

    thresholds = np.full(n_genes, np.inf, dtype=np.float64)
    np.minimum.at(thresholds, reciprocal_sources, reciprocal_scores)
    np.minimum.at(thresholds, reciprocal_targets, reciprocal_scores)
    keep = (
        (queries != targets)
        & np.isfinite(scores)
        & (scores >= np.minimum(thresholds[queries], thresholds[targets]))
    )
    return deduplicate_undirected_edges(
        gene_names,
        queries[keep],
        targets[keep],
        scores[keep],
    )


def read_index_clusters(path: str, gene_names: Sequence[str]) -> list[list[int]]:
    """Read a whitespace-delimited cluster file into global gene IDs."""
    gene_to_id = {gene: idx for idx, gene in enumerate(gene_names)}
    clusters = []
    with open(path, "r") as handle:
        for line in handle:
            members = [gene_to_id[gene] for gene in line.split()]
            if members:
                clusters.append(members)
    return clusters


def build_singleton_assignment_edges(
    gene_names: Sequence[str],
    clusters: Sequence[Sequence[int]],
    hit_queries,
    hit_targets,
    hit_scores,
) -> IndexedEdges:
    """Attach each singleton to its best-supported non-singleton cluster."""
    n_genes = len(gene_names)
    cluster_by_gene = np.full(n_genes, -1, dtype=np.int32)
    singleton = np.ones(n_genes, dtype=bool)
    for cluster_id, cluster in enumerate(clusters):
        members = np.asarray(cluster, dtype=np.int32)
        if len(members) >= 2:
            cluster_by_gene[members] = cluster_id
            singleton[members] = False

    queries = np.asarray(hit_queries, dtype=np.int32)
    targets = np.asarray(hit_targets, dtype=np.int32)
    scores = np.asarray(hit_scores, dtype=np.float64)
    if not (len(queries) == len(targets) == len(scores)):
        raise ValueError("hit query, target, and score arrays must align")
    if len(queries) == 0:
        return _empty_edges(gene_names)

    eligible = (
        singleton[queries]
        & (cluster_by_gene[targets] >= 0)
        & np.isfinite(scores)
        & (scores > 0.0)
    )
    if not eligible.any():
        return _empty_edges(gene_names)

    eligible_rows = np.flatnonzero(eligible)
    best_scores = np.full(n_genes, -np.inf, dtype=np.float64)
    np.maximum.at(best_scores, queries[eligible_rows], scores[eligible_rows])
    winning_rows = eligible_rows[
        scores[eligible_rows] == best_scores[queries[eligible_rows]]
    ]
    winning_queries, first = np.unique(
        queries[winning_rows], return_index=True
    )
    best_clusters = np.full(n_genes, -1, dtype=np.int32)
    best_clusters[winning_queries] = cluster_by_gene[
        targets[winning_rows[first]]
    ]

    keep = eligible & (
        cluster_by_gene[targets] == best_clusters[queries]
    )
    return deduplicate_undirected_edges(
        gene_names,
        queries[keep],
        targets[keep],
        scores[keep],
    )
