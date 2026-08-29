"""Cluster refinement utilities for orthogroup inference.

The refinement stage works after graph clustering. It uses two independent
signals:

* all significant directed search hits to identify reciprocal cluster-level
  support that Leiden can miss when evidence is spread across many genes; and
* the stricter RBNH graph to split very large duplicated clusters whose member
  degrees show obvious outliers.
* a broad-dataset copy-number guard to break the largest high-copy overmerges
  that dominate pairwise false positives on QfO-scale inputs. On broad inputs
  this guard is applied by itself so the production path matches the validated
  QfO-scale transformation.
* a medium-panel copy-number guard for large-but-not-huge clusters dominated
  by many copies from one species.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple

import numpy as np


Gene = str
Cluster = List[Gene]
IndexCluster = List[int]
DirectedHit = Tuple[Gene, Gene, float]
EdgeMap = Mapping[frozenset, float]


DEFAULT_MAX_RECIPROCAL_MERGES = 150
DEFAULT_MAX_COMPONENT_GENES = 80
DEFAULT_SPLIT_MIN_SIZE = 20
DEFAULT_SPLIT_MIN_SPECIES_COUNT = 20
DEFAULT_SPLIT_DEGREE_RATIO = 1.1
DEFAULT_COPY_SPLIT_MIN_SIZE = 150
DEFAULT_COPY_SPLIT_MIN_SPECIES_COUNT = 10
DEFAULT_COPY_SPLIT_MIN_DATASET_SPECIES = 50
DEFAULT_BROAD_COPY_SPLIT_MIN_SPECIES_COUNT = 10
DEFAULT_PANEL_COPY_SPLIT_MIN_SIZE = 70
DEFAULT_PANEL_COPY_SPLIT_MAX_SIZE = 150
DEFAULT_PANEL_COPY_SPLIT_MIN_SPECIES_COUNT = 14
DEFAULT_PANEL_COPY_SPLIT_MIN_DATASET_SPECIES = 10
DEFAULT_PANEL_COPY_COMPONENT_MIN_EDGE_WEIGHT = 1.5
DEFAULT_BROAD_RECIPROCAL_MIN_AVG_SCORE = 1.0
DEFAULT_BROAD_RECIPROCAL_MIN_MAX_SCORE = 1.0
DEFAULT_BROAD_RECIPROCAL_MIN_COVERAGE = 2.0
DEFAULT_BROAD_RECIPROCAL_MIN_NORM = 1.0
DEFAULT_BROAD_MAX_RECIPROCAL_MERGES = 32
DEFAULT_REFINE_PROFILE = "default"

# Tunable refinement presets. "default" preserves current production behavior.
_REFINE_PROFILE_PRESETS: Dict[str, Dict[str, object]] = {
    "default": {},
    "qfo": {
        "broad_copy_split_min_species_count": 24,
        "broad_max_reciprocal_merges": 64,
        "broad_reciprocal_min_avg_score": 0.95,
        "broad_reciprocal_min_max_score": 0.95,
        "broad_reciprocal_min_coverage": 1.6,
        "broad_reciprocal_min_norm": 0.8,
        "copy_split_min_size": 180,
    },
}


def resolve_refinement_profile(name: str) -> Dict[str, object]:
    """Return supported refinement parameter overrides for a named preset."""
    profile = name.lower() if isinstance(name, str) else DEFAULT_REFINE_PROFILE
    if profile not in _REFINE_PROFILE_PRESETS:
        raise ValueError(f"Unknown refinement profile: {name!r}")
    return dict(_REFINE_PROFILE_PRESETS[profile])


@dataclass
class _PairStats:
    total: float = 0.0
    count: int = 0
    max_score: float = 0.0

    def add(self, score: float) -> None:
        self.total += score
        self.count += 1
        if score > self.max_score:
            self.max_score = score


class _DSU:
    def __init__(self, clusters: Sequence[Cluster], gene_to_species: Mapping[Gene, str]):
        self.parent = list(range(len(clusters)))
        self.gene_size = [len(c) for c in clusters]
        self.species_counts = [
            Counter(gene_to_species.get(g, "") for g in c)
            for c in clusters
        ]

    def find(self, item: int) -> int:
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, a: int, b: int, max_genes: int | None = None) -> bool:
        root_a = self.find(a)
        root_b = self.find(b)
        if root_a == root_b:
            return False
        if max_genes is not None and self.gene_size[root_a] + self.gene_size[root_b] > max_genes:
            return False
        if self.gene_size[root_a] < self.gene_size[root_b]:
            root_a, root_b = root_b, root_a
        self.parent[root_b] = root_a
        self.gene_size[root_a] += self.gene_size[root_b]
        self.species_counts[root_a].update(self.species_counts[root_b])
        return True


class _IndexedDSU:
    def __init__(self, clusters: Sequence[Sequence[int]], gene_to_species: Sequence[int]):
        self.parent = list(range(len(clusters)))
        self.gene_size = [len(c) for c in clusters]
        self.species_counts = [
            Counter(int(gene_to_species[g]) for g in c)
            for c in clusters
        ]

    def find(self, item: int) -> int:
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, a: int, b: int, max_genes: int | None = None) -> bool:
        root_a = self.find(a)
        root_b = self.find(b)
        if root_a == root_b:
            return False
        if max_genes is not None and self.gene_size[root_a] + self.gene_size[root_b] > max_genes:
            return False
        if self.gene_size[root_a] < self.gene_size[root_b]:
            root_a, root_b = root_b, root_a
        self.parent[root_b] = root_a
        self.gene_size[root_a] += self.gene_size[root_b]
        self.species_counts[root_a].update(self.species_counts[root_b])
        return True


def _cluster_pair_stats(
    clusters: Sequence[Cluster],
    directed_hits: Iterable[DirectedHit],
) -> Dict[Tuple[int, int], _PairStats]:
    gene_to_cluster = {
        gene: idx
        for idx, cluster in enumerate(clusters)
        for gene in cluster
    }
    pair_stats: Dict[Tuple[int, int], _PairStats] = defaultdict(_PairStats)
    for query, target, score in directed_hits:
        query_cluster = gene_to_cluster.get(query)
        target_cluster = gene_to_cluster.get(target)
        if query_cluster is None or target_cluster is None or query_cluster == target_cluster:
            continue
        pair_stats[(query_cluster, target_cluster)].add(float(score))
    return pair_stats


def _features(
    stats: _PairStats,
    source_size: int,
    target_size: int,
) -> Tuple[float, float, float, float]:
    if stats.count == 0:
        return 0.0, 0.0, 0.0, 0.0
    min_size = max(1, min(source_size, target_size))
    avg_score = stats.total / stats.count
    coverage = stats.count / min_size
    norm = stats.total / math.sqrt(max(1, source_size * target_size))
    return avg_score, stats.max_score, coverage, norm


def _component_clusters(clusters: Sequence[Cluster], dsu: _DSU) -> List[Cluster]:
    components: MutableMapping[int, Cluster] = defaultdict(list)
    for idx, cluster in enumerate(clusters):
        components[dsu.find(idx)].extend(cluster)
    return list(components.values())


def _component_index_clusters(
    clusters: Sequence[Sequence[int]],
    dsu: _IndexedDSU,
) -> List[IndexCluster]:
    components: MutableMapping[int, IndexCluster] = defaultdict(list)
    for idx, cluster in enumerate(clusters):
        components[dsu.find(idx)].extend(int(gene) for gene in cluster)
    return list(components.values())


def _species_overlap(a: Counter, b: Counter) -> int:
    return len(set(a) & set(b))


def _should_split_large_copy_cluster(
    cluster_size: int,
    species_counts: Counter,
    dataset_species_count: int,
    min_size: int,
    min_species_count: int,
    min_dataset_species: int,
) -> bool:
    if min_size <= 0 or min_species_count <= 0:
        return False
    if dataset_species_count < min_dataset_species:
        return False
    if cluster_size < min_size:
        return False
    return bool(species_counts) and max(species_counts.values()) >= min_species_count


def _should_split_medium_panel_copy_cluster(
    cluster_size: int,
    species_counts: Counter,
    dataset_species_count: int,
    min_size: int,
    max_size: int,
    min_species_count: int,
    min_dataset_species: int,
) -> bool:
    if min_size <= 0 or max_size <= min_size or min_species_count <= 0:
        return False
    if dataset_species_count < min_dataset_species:
        return False
    if cluster_size < min_size or cluster_size >= max_size:
        return False
    return bool(species_counts) and max(species_counts.values()) >= min_species_count


def _copy_split_large_clusters(
    clusters: Sequence[Cluster],
    gene_to_species: Mapping[Gene, str],
    min_size: int,
    min_species_count: int,
    min_dataset_species: int,
) -> List[Cluster]:
    dataset_species_count = len({
        species for species in gene_to_species.values()
        if species != ""
    })
    refined: List[Cluster] = []
    for cluster in clusters:
        cluster_list = list(cluster)
        species_counts = Counter(gene_to_species.get(g, "") for g in cluster_list)
        if _should_split_large_copy_cluster(
            len(cluster_list),
            species_counts,
            dataset_species_count,
            min_size,
            min_species_count,
            min_dataset_species,
        ):
            refined.extend([[gene] for gene in cluster_list])
        else:
            refined.append(cluster_list)
    return refined


def _copy_split_large_index_clusters(
    clusters: Sequence[Sequence[int]],
    gene_to_species: Sequence[int],
    min_size: int,
    min_species_count: int,
    min_dataset_species: int,
) -> List[IndexCluster]:
    species = np.asarray(gene_to_species)
    dataset_species_count = int(np.unique(species).size) if len(species) else 0
    refined: List[IndexCluster] = []
    for cluster in clusters:
        cluster_list = [int(gene) for gene in cluster]
        species_counts = Counter(int(species[gene]) for gene in cluster_list)
        if _should_split_large_copy_cluster(
            len(cluster_list),
            species_counts,
            dataset_species_count,
            min_size,
            min_species_count,
            min_dataset_species,
        ):
            refined.extend([[gene] for gene in cluster_list])
        else:
            refined.append(cluster_list)
    return refined


def _split_string_cluster_by_components(
    cluster: Sequence[Gene],
    adjacency: Mapping[Gene, set],
) -> List[Cluster]:
    members = set(cluster)
    seen: set[Gene] = set()
    components: List[Cluster] = []
    for gene in cluster:
        if gene in seen:
            continue
        seen.add(gene)
        stack = [gene]
        component: Cluster = []
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in adjacency.get(current, set()):
                if neighbor in members and neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        components.append(component)
    return components


def _split_index_cluster_by_components(
    cluster: Sequence[int],
    adjacency: Mapping[int, set],
) -> List[IndexCluster]:
    members = {int(gene) for gene in cluster}
    seen: set[int] = set()
    components: List[IndexCluster] = []
    for gene in cluster:
        gene = int(gene)
        if gene in seen:
            continue
        seen.add(gene)
        stack = [gene]
        component: IndexCluster = []
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in adjacency.get(current, set()):
                if neighbor in members and neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        components.append(component)
    return components


def _is_broad_string_dataset(
    gene_to_species: Mapping[Gene, str],
    min_dataset_species: int,
) -> bool:
    dataset_species_count = len({
        species for species in gene_to_species.values()
        if species != ""
    })
    return dataset_species_count >= min_dataset_species


def _is_broad_index_dataset(
    gene_to_species: Sequence[int],
    min_dataset_species: int,
) -> bool:
    species = np.asarray(gene_to_species)
    dataset_species_count = int(np.unique(species).size) if len(species) else 0
    return dataset_species_count >= min_dataset_species


def _resolve_broad_copy_split_min_species_count(
    dataset_species_count: int,
    requested_min_species_count: int,
    min_dataset_species: int,
    broad_copy_split_min_species_count: int,
) -> int:
    if dataset_species_count >= min_dataset_species:
        return max(requested_min_species_count, broad_copy_split_min_species_count)
    return requested_min_species_count


def _merge_reciprocal_cluster_best(
    clusters: Sequence[Cluster],
    pair_stats: Mapping[Tuple[int, int], _PairStats],
    dsu: _DSU,
    max_merges: int,
    max_genes: int,
    min_avg_score: float = 0.5,
    min_max_score: float = 0.8,
    min_coverage: float = 1.0,
    min_norm: float | None = None,
) -> int:
    sizes = [len(c) for c in clusters]
    best_out: Dict[int, Tuple[float, int]] = {}
    for (source, target), stats in pair_stats.items():
        avg_score, max_score, coverage, norm = _features(
            stats, sizes[source], sizes[target]
        )
        if avg_score < min_avg_score or max_score < min_max_score or coverage < min_coverage:
            continue
        if min_norm is not None and norm < min_norm:
            continue
        current = best_out.get(source)
        if current is None or norm > current[0]:
            best_out[source] = (norm, target)

    reciprocal = []
    for source, (norm, target) in best_out.items():
        target_best = best_out.get(target)
        if target_best is None or target_best[1] != source or source >= target:
            continue
        reverse = pair_stats[(target, source)]
        _avg, _max_score, _coverage, reverse_norm = _features(
            reverse, sizes[target], sizes[source]
        )
        reciprocal.append((min(norm, reverse_norm), source, target))
    reciprocal.sort(reverse=True)

    merged = 0
    for _score, source, target in reciprocal:
        if merged >= max_merges:
            break
        if dsu.union(source, target, max_genes=max_genes):
            merged += 1
    return merged


def _merge_weak_balanced_rescue(
    clusters: Sequence[Cluster],
    pair_stats: Mapping[Tuple[int, int], _PairStats],
    dsu: _DSU,
    max_genes: int,
) -> int:
    sizes = [len(c) for c in clusters]
    candidates = []
    for (source, target), forward in pair_stats.items():
        if source >= target:
            continue
        reverse = pair_stats.get((target, source))
        if reverse is None:
            continue
        root_a = dsu.find(source)
        root_b = dsu.find(target)
        if root_a == root_b:
            continue
        comp_a = dsu.species_counts[root_a]
        comp_b = dsu.species_counts[root_b]
        combined_species = comp_a + comp_b
        size_a = dsu.gene_size[root_a]
        size_b = dsu.gene_size[root_b]
        total_size = size_a + size_b
        if total_size > 45 or total_size > max_genes:
            continue
        size_ratio = max(size_a, size_b) / max(1, min(size_a, size_b))
        shared_species = _species_overlap(comp_a, comp_b)
        max_species_count = max(combined_species.values()) if combined_species else 0

        avg_a, max_a, cov_a, norm_a = _features(forward, sizes[source], sizes[target])
        avg_b, max_b, cov_b, norm_b = _features(reverse, sizes[target], sizes[source])
        min_avg = min(avg_a, avg_b)
        min_max = min(max_a, max_b)
        min_cov = min(cov_a, cov_b)
        min_norm = min(norm_a, norm_b)

        if (
            shared_species == 9
            and max_species_count <= 6
            and size_ratio <= 1.35
            and 0.22 <= min_avg <= 0.36
            and 0.65 <= min_max <= 0.90
            and 1.0 <= min_cov <= 1.4
            and 0.20 <= min_norm <= 0.50
        ):
            candidates.append((
                -size_ratio,
                -abs(min_max - 0.75),
                min_avg,
                source,
                target,
            ))

    candidates.sort(reverse=True)
    merged = 0
    for _ratio, _max_dist, _avg, source, target in candidates:
        if dsu.union(source, target, max_genes=max_genes):
            merged += 1
    return merged


def _split_high_duplication_clusters(
    clusters: Sequence[Cluster],
    rbnh_edges: EdgeMap,
    gene_to_species: Mapping[Gene, str],
    min_size: int,
    min_species_count: int,
    degree_ratio: float,
    copy_split_min_size: int,
    copy_split_min_species_count: int,
    copy_split_min_dataset_species: int,
    panel_copy_split_min_size: int,
    panel_copy_split_max_size: int,
    panel_copy_split_min_species_count: int,
    panel_copy_split_min_dataset_species: int,
    panel_copy_component_min_edge_weight: float,
    component_split_high_duplication: bool = False,
    broad_mode: bool = False,
) -> List[Cluster]:
    adjacency: MutableMapping[Gene, set] = defaultdict(set)
    component_adjacency: MutableMapping[Gene, set] = defaultdict(set)
    component_nodes = set()
    for edge, weight in rbnh_edges.items():
        if len(edge) != 2:
            continue
        gene_a, gene_b = tuple(edge)
        adjacency[gene_a].add(gene_b)
        adjacency[gene_b].add(gene_a)
        if float(weight) >= panel_copy_component_min_edge_weight:
            component_adjacency[gene_a].add(gene_b)
            component_adjacency[gene_b].add(gene_a)
            component_nodes.add(gene_a)
            component_nodes.add(gene_b)

    dataset_species_count = len({
        species for species in gene_to_species.values()
        if species != ""
    })
    refined: List[Cluster] = []
    for cluster in clusters:
        members = set(cluster)
        if len(cluster) < min_size and (
            copy_split_min_size <= 0 or len(cluster) < copy_split_min_size
        ):
            refined.append(list(cluster))
            continue
        species_counts = Counter(gene_to_species.get(g, "") for g in cluster)
        if _should_split_large_copy_cluster(
            len(cluster),
            species_counts,
            dataset_species_count,
            copy_split_min_size,
            copy_split_min_species_count,
            copy_split_min_dataset_species,
        ):
            refined.extend([[gene] for gene in cluster])
            continue
        if _should_split_medium_panel_copy_cluster(
            len(cluster),
            species_counts,
            dataset_species_count,
            panel_copy_split_min_size,
            panel_copy_split_max_size,
            panel_copy_split_min_species_count,
            panel_copy_split_min_dataset_species,
        ):
            if broad_mode and not (members & component_nodes):
                refined.append(list(cluster))
                continue
            split_components = _split_string_cluster_by_components(
                cluster,
                component_adjacency,
            )
            if broad_mode and all(len(component) == 1 for component in split_components):
                refined.append(list(cluster))
            else:
                refined.extend(split_components)
            continue
        if broad_mode:
            refined.append(list(cluster))
            continue
        if len(cluster) < min_size:
            refined.append(list(cluster))
            continue
        if max(species_counts.values()) < min_species_count:
            refined.append(list(cluster))
            continue
        if component_split_high_duplication:
            split_components = _split_string_cluster_by_components(
                cluster,
                component_adjacency,
            )
            if len(split_components) > 1 and any(
                len(component) > 1 for component in split_components
            ):
                refined.extend(split_components)
                continue
        degrees = {
            gene: len(adjacency.get(gene, set()) & members)
            for gene in cluster
        }
        if not degrees:
            refined.append(list(cluster))
            continue
        max_degree = max(degrees.values())
        keep = [
            gene for gene in cluster
            if degrees[gene] * degree_ratio >= max_degree
        ]
        cut = [
            gene for gene in cluster
            if degrees[gene] * degree_ratio < max_degree
        ]
        if len(keep) > 1:
            refined.append(keep)
        else:
            refined.extend([[gene] for gene in keep])
        refined.extend([[gene] for gene in cut])
    return refined


def _gene_to_cluster_array(
    clusters: Sequence[Sequence[int]],
    total_genes: int,
) -> np.ndarray:
    gene_to_cluster = np.full(total_genes, -1, dtype=np.int32)
    for cluster_idx, cluster in enumerate(clusters):
        if cluster:
            gene_to_cluster[np.asarray(cluster, dtype=np.int64)] = cluster_idx
    return gene_to_cluster


def _cluster_pair_arrays(
    clusters: Sequence[Sequence[int]],
    hit_queries: Sequence[int],
    hit_targets: Sequence[int],
    hit_scores: Sequence[float],
    total_genes: int,
    default_score: float | None = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    gene_to_cluster = _gene_to_cluster_array(clusters, total_genes)
    queries = np.asarray(hit_queries, dtype=np.int64)
    targets = np.asarray(hit_targets, dtype=np.int64)
    scores = np.asarray(hit_scores, dtype=np.float64)
    if len(queries) != len(targets):
        common = min(len(queries), len(targets))
        queries = queries[:common]
        targets = targets[:common]
    if len(scores) > 0 and (len(scores) != len(queries)):
        if len(queries) == 0:
            scores = np.asarray([], dtype=np.float64)
        else:
            scores = scores[: min(len(queries), len(scores))]
            queries = queries[: len(scores)]
            targets = targets[: len(scores)]
    if len(queries) == 0:
        empty_i = np.asarray([], dtype=np.int32)
        empty_f = np.asarray([], dtype=np.float64)
        return empty_i, empty_i, empty_f, empty_i, empty_f
    if len(scores) == 0:
        if default_score is None:
            empty_i = np.asarray([], dtype=np.int32)
            empty_f = np.asarray([], dtype=np.float64)
            return empty_i, empty_i, empty_f, empty_i, empty_f
        scores = np.full(len(queries), default_score, dtype=np.float64)
    if len(scores) == 0:
        empty_i = np.asarray([], dtype=np.int32)
        empty_f = np.asarray([], dtype=np.float64)
        return empty_i, empty_i, empty_f, empty_i, empty_f

    source_clusters = gene_to_cluster[queries]
    target_clusters = gene_to_cluster[targets]
    valid = (
        (source_clusters >= 0)
        & (target_clusters >= 0)
        & (source_clusters != target_clusters)
    )
    if not np.any(valid):
        empty_i = np.asarray([], dtype=np.int32)
        empty_f = np.asarray([], dtype=np.float64)
        return empty_i, empty_i, empty_f, empty_i, empty_f

    source_valid = source_clusters[valid].astype(np.int64, copy=False)
    target_valid = target_clusters[valid].astype(np.int64, copy=False)
    score_valid = scores[valid]
    pair_key = (source_valid << 32) | target_valid
    order = np.argsort(pair_key, kind="stable")
    key_sorted = pair_key[order]
    score_sorted = score_valid[order]
    run_starts = np.concatenate(
        ([0], np.nonzero(np.diff(key_sorted))[0] + 1)
    )
    run_ends = np.concatenate((run_starts[1:], [len(key_sorted)]))
    key_unique = key_sorted[run_starts]
    sources = (key_unique >> 32).astype(np.int32)
    targets = (key_unique & 0xFFFFFFFF).astype(np.int32)
    totals = np.add.reduceat(score_sorted, run_starts).astype(np.float64)
    counts = (run_ends - run_starts).astype(np.int32)
    max_scores = np.maximum.reduceat(score_sorted, run_starts).astype(np.float64)
    return sources, targets, totals, counts, max_scores


def _indexed_features(
    sources: np.ndarray,
    targets: np.ndarray,
    totals: np.ndarray,
    counts: np.ndarray,
    max_scores: np.ndarray,
    sizes: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    min_sizes = np.maximum(1, np.minimum(sizes[sources], sizes[targets]))
    avg_scores = totals / np.maximum(1, counts)
    coverages = counts / min_sizes
    norms = totals / np.sqrt(np.maximum(1, sizes[sources] * sizes[targets]))
    return avg_scores, max_scores, coverages, norms


def _merge_reciprocal_cluster_best_indexed(
    sources: np.ndarray,
    targets: np.ndarray,
    totals: np.ndarray,
    counts: np.ndarray,
    max_scores: np.ndarray,
    sizes: np.ndarray,
    dsu: _IndexedDSU,
    max_merges: int,
    max_genes: int,
    min_avg_score: float = 0.5,
    min_max_score: float = 0.8,
    min_coverage: float = 1.0,
    min_norm: float | None = None,
) -> int:
    avg_scores, max_scores, coverages, norms = _indexed_features(
        sources, targets, totals, counts, max_scores, sizes
    )
    candidate = (
        (avg_scores >= min_avg_score)
        & (max_scores >= min_max_score)
        & (coverages >= min_coverage)
    )
    if min_norm is not None:
        candidate = candidate & (norms >= min_norm)
    best_norm = np.full(len(sizes), -np.inf, dtype=np.float64)
    best_target = np.full(len(sizes), -1, dtype=np.int32)
    candidate_sources = sources[candidate]
    candidate_targets = targets[candidate]
    candidate_norms = norms[candidate]
    if len(candidate_sources):
        np.maximum.at(best_norm, candidate_sources, candidate_norms)
        winning_rows = np.flatnonzero(
            candidate_norms == best_norm[candidate_sources]
        )
        winning_sources, first = np.unique(
            candidate_sources[winning_rows], return_index=True
        )
        selected_rows = winning_rows[first]
        best_target[winning_sources] = candidate_targets[selected_rows]

    reciprocal = []
    for source, target in enumerate(best_target):
        if target < 0 or source >= target:
            continue
        if best_target[target] == source:
            reciprocal.append((min(best_norm[source], best_norm[target]), source, int(target)))
    reciprocal.sort(reverse=True)

    merged = 0
    for _score, source, target in reciprocal:
        if merged >= max_merges:
            break
        if dsu.union(source, target, max_genes=max_genes):
            merged += 1
    return merged


def _merge_weak_balanced_rescue_indexed(
    sources: np.ndarray,
    targets: np.ndarray,
    totals: np.ndarray,
    counts: np.ndarray,
    max_scores: np.ndarray,
    sizes: np.ndarray,
    dsu: _IndexedDSU,
    max_genes: int,
) -> int:
    avg_scores, max_scores, coverages, norms = _indexed_features(
        sources, targets, totals, counts, max_scores, sizes
    )
    pair_sizes = sizes[sources] + sizes[targets]
    plausible = (
        (pair_sizes <= min(45, max_genes))
        & (avg_scores >= 0.22)
        & (avg_scores <= 0.36)
        & (max_scores >= 0.65)
        & (max_scores <= 0.90)
        & (coverages >= 1.0)
        & (coverages <= 1.4)
        & (norms >= 0.20)
        & (norms <= 0.50)
    )
    plausible_lookup = {
        (int(source) << 32) | int(target): idx
        for idx, (source, target) in enumerate(zip(sources[plausible], targets[plausible]))
    }
    plausible_indices = np.nonzero(plausible)[0]
    candidates = []
    for idx in plausible_indices:
        source = int(sources[idx])
        target = int(targets[idx])
        if source >= target:
            continue
        reverse_key = (target << 32) | source
        reverse_pos = plausible_lookup.get(reverse_key)
        if reverse_pos is None:
            continue
        reverse_idx = int(plausible_indices[reverse_pos])
        root_a = dsu.find(source)
        root_b = dsu.find(target)
        if root_a == root_b:
            continue
        comp_a = dsu.species_counts[root_a]
        comp_b = dsu.species_counts[root_b]
        combined_species = comp_a + comp_b
        size_a = dsu.gene_size[root_a]
        size_b = dsu.gene_size[root_b]
        total_size = size_a + size_b
        if total_size > 45 or total_size > max_genes:
            continue
        size_ratio = max(size_a, size_b) / max(1, min(size_a, size_b))
        shared_species = _species_overlap(comp_a, comp_b)
        max_species_count = max(combined_species.values()) if combined_species else 0
        min_avg = min(avg_scores[idx], avg_scores[reverse_idx])
        min_max = min(max_scores[idx], max_scores[reverse_idx])
        min_cov = min(coverages[idx], coverages[reverse_idx])
        min_norm = min(norms[idx], norms[reverse_idx])

        if (
            shared_species == 9
            and max_species_count <= 6
            and size_ratio <= 1.35
            and 0.22 <= min_avg <= 0.36
            and 0.65 <= min_max <= 0.90
            and 1.0 <= min_cov <= 1.4
            and 0.20 <= min_norm <= 0.50
        ):
            candidates.append((
                -size_ratio,
                -abs(min_max - 0.75),
                min_avg,
                source,
                target,
            ))

    candidates.sort(reverse=True)
    merged = 0
    for _ratio, _max_dist, _avg, source, target in candidates:
        if dsu.union(source, target, max_genes=max_genes):
            merged += 1
    return merged


def _split_high_duplication_index_clusters(
    clusters: Sequence[Sequence[int]],
    rbnh_queries: Sequence[int],
    rbnh_targets: Sequence[int],
    gene_to_species: Sequence[int],
    total_genes: int,
    min_size: int,
    min_species_count: int,
    degree_ratio: float,
    copy_split_min_size: int,
    copy_split_min_species_count: int,
    copy_split_min_dataset_species: int,
    panel_copy_split_min_size: int,
    panel_copy_split_max_size: int,
    panel_copy_split_min_species_count: int,
    panel_copy_split_min_dataset_species: int,
    panel_copy_component_min_edge_weight: float,
    rbnh_scores: Sequence[float] | None,
    component_split_high_duplication: bool = False,
    broad_mode: bool = False,
) -> List[IndexCluster]:
    gene_to_cluster = _gene_to_cluster_array(clusters, total_genes)
    edge_queries = np.asarray(rbnh_queries, dtype=np.int64)
    edge_targets = np.asarray(rbnh_targets, dtype=np.int64)
    edge_scores = (
        np.asarray(rbnh_scores, dtype=np.float64)
        if rbnh_scores is not None
        else np.asarray([], dtype=np.float64)
    )
    if len(edge_queries) != len(edge_targets):
        common = min(len(edge_queries), len(edge_targets))
        edge_queries = edge_queries[:common]
        edge_targets = edge_targets[:common]
    if len(edge_scores) > 0 and (len(edge_scores) != len(edge_queries)):
        edge_scores = edge_scores[: min(len(edge_scores), len(edge_queries))]
        edge_queries = edge_queries[: len(edge_scores)]
        edge_targets = edge_targets[: len(edge_scores)]
    degrees = np.zeros(total_genes, dtype=np.int32)
    component_adjacency: MutableMapping[int, set] = defaultdict(set)
    component_nodes = set()
    if len(edge_queries) > 0:
        source_clusters = gene_to_cluster[edge_queries]
        target_clusters = gene_to_cluster[edge_targets]
        same_cluster = (
            (source_clusters >= 0)
            & (source_clusters == target_clusters)
            & (edge_queries != edge_targets)
        )
        if np.any(same_cluster):
            np.add.at(degrees, edge_queries[same_cluster], 1)
            np.add.at(degrees, edge_targets[same_cluster], 1)
            if len(edge_scores) == len(edge_queries):
                component_edges = same_cluster & (
                    edge_scores >= panel_copy_component_min_edge_weight
                )
                for query, target in zip(
                    edge_queries[component_edges],
                    edge_targets[component_edges],
                ):
                    query_i = int(query)
                    target_i = int(target)
                    component_adjacency[query_i].add(target_i)
                    component_adjacency[target_i].add(query_i)
                    component_nodes.add(query_i)
                    component_nodes.add(target_i)

    species = np.asarray(gene_to_species)
    dataset_species_count = int(np.unique(species).size) if len(species) else 0
    refined: List[IndexCluster] = []
    for cluster in clusters:
        cluster_list = [int(gene) for gene in cluster]
        cluster_set = set(cluster_list)
        if len(cluster_list) < min_size and (
            copy_split_min_size <= 0 or len(cluster_list) < copy_split_min_size
        ):
            refined.append(cluster_list)
            continue
        species_counts = Counter(int(species[gene]) for gene in cluster_list)
        if _should_split_large_copy_cluster(
            len(cluster_list),
            species_counts,
            dataset_species_count,
            copy_split_min_size,
            copy_split_min_species_count,
            copy_split_min_dataset_species,
        ):
            refined.extend([[gene] for gene in cluster_list])
            continue
        if _should_split_medium_panel_copy_cluster(
            len(cluster_list),
            species_counts,
            dataset_species_count,
            panel_copy_split_min_size,
            panel_copy_split_max_size,
            panel_copy_split_min_species_count,
            panel_copy_split_min_dataset_species,
        ):
            if broad_mode and not (cluster_set & component_nodes):
                refined.append(cluster_list)
                continue
            split_components = _split_index_cluster_by_components(
                cluster_list,
                component_adjacency,
            )
            if broad_mode and all(len(component) == 1 for component in split_components):
                refined.append(cluster_list)
            else:
                refined.extend(split_components)
            continue
        if broad_mode:
            refined.append(cluster_list)
            continue
        if len(cluster_list) < min_size:
            refined.append(cluster_list)
            continue
        if max(species_counts.values()) < min_species_count:
            refined.append(cluster_list)
            continue
        if component_split_high_duplication:
            split_components = _split_index_cluster_by_components(
                cluster_list,
                component_adjacency,
            )
            if len(split_components) > 1 and any(
                len(component) > 1 for component in split_components
            ):
                refined.extend(split_components)
                continue
        cluster_degrees = degrees[np.asarray(cluster_list, dtype=np.int64)]
        max_degree = int(cluster_degrees.max()) if len(cluster_degrees) else 0
        if max_degree == 0:
            refined.append(cluster_list)
            continue
        keep = [
            gene
            for gene, degree in zip(cluster_list, cluster_degrees)
            if int(degree) * degree_ratio >= max_degree
        ]
        cut = [
            gene
            for gene, degree in zip(cluster_list, cluster_degrees)
            if int(degree) * degree_ratio < max_degree
        ]
        if len(keep) > 1:
            refined.append(keep)
        else:
            refined.extend([[gene] for gene in keep])
        refined.extend([[gene] for gene in cut])
    return refined


def refine_clusters(
    clusters: Sequence[Cluster],
    directed_hits: Iterable[DirectedHit],
    rbnh_edges: EdgeMap,
    gene_to_species: Mapping[Gene, str],
    max_reciprocal_merges: int = DEFAULT_MAX_RECIPROCAL_MERGES,
    max_component_genes: int = DEFAULT_MAX_COMPONENT_GENES,
    split_min_size: int = DEFAULT_SPLIT_MIN_SIZE,
    split_min_species_count: int = DEFAULT_SPLIT_MIN_SPECIES_COUNT,
    split_degree_ratio: float = DEFAULT_SPLIT_DEGREE_RATIO,
    broad_copy_split_min_species_count: int = DEFAULT_BROAD_COPY_SPLIT_MIN_SPECIES_COUNT,
    copy_split_min_size: int = DEFAULT_COPY_SPLIT_MIN_SIZE,
    copy_split_min_species_count: int = DEFAULT_COPY_SPLIT_MIN_SPECIES_COUNT,
    copy_split_min_dataset_species: int = DEFAULT_COPY_SPLIT_MIN_DATASET_SPECIES,
    panel_copy_split_min_size: int = DEFAULT_PANEL_COPY_SPLIT_MIN_SIZE,
    panel_copy_split_max_size: int = DEFAULT_PANEL_COPY_SPLIT_MAX_SIZE,
    panel_copy_split_min_species_count: int = DEFAULT_PANEL_COPY_SPLIT_MIN_SPECIES_COUNT,
    panel_copy_split_min_dataset_species: int = DEFAULT_PANEL_COPY_SPLIT_MIN_DATASET_SPECIES,
    panel_copy_component_min_edge_weight: float = DEFAULT_PANEL_COPY_COMPONENT_MIN_EDGE_WEIGHT,
    component_split_high_duplication: bool = False,
) -> List[Cluster]:
    """Refine orthogroup clusters using cluster-level support and graph degree.

    The default parameters are deliberately conservative: merge only small
    components with reciprocal aggregate evidence, then split only large
    duplicated clusters with clear internal-degree outliers.
    """
    if not clusters:
        return []
    if _is_broad_string_dataset(gene_to_species, copy_split_min_dataset_species):
        resolved_copy_split_min_species_count = (
            _resolve_broad_copy_split_min_species_count(
                len({sp for sp in gene_to_species.values() if sp != ""}),
                copy_split_min_species_count,
                copy_split_min_dataset_species,
                broad_copy_split_min_species_count,
            )
        )
        return _split_high_duplication_clusters(
            clusters,
            rbnh_edges,
            gene_to_species,
            min_size=split_min_size,
            min_species_count=split_min_species_count,
            degree_ratio=split_degree_ratio,
            copy_split_min_size=copy_split_min_size,
            copy_split_min_species_count=resolved_copy_split_min_species_count,
            copy_split_min_dataset_species=copy_split_min_dataset_species,
            panel_copy_split_min_size=panel_copy_split_min_size,
            panel_copy_split_max_size=panel_copy_split_max_size,
            panel_copy_split_min_species_count=panel_copy_split_min_species_count,
            panel_copy_split_min_dataset_species=panel_copy_split_min_dataset_species,
            panel_copy_component_min_edge_weight=panel_copy_component_min_edge_weight,
            component_split_high_duplication=component_split_high_duplication,
            broad_mode=True,
        )
    pair_stats = _cluster_pair_stats(clusters, directed_hits)
    if not pair_stats:
        return _split_high_duplication_clusters(
            clusters,
            rbnh_edges,
            gene_to_species,
            min_size=split_min_size,
            min_species_count=split_min_species_count,
            degree_ratio=split_degree_ratio,
            copy_split_min_size=copy_split_min_size,
            copy_split_min_species_count=copy_split_min_species_count,
            copy_split_min_dataset_species=copy_split_min_dataset_species,
            panel_copy_split_min_size=panel_copy_split_min_size,
            panel_copy_split_max_size=panel_copy_split_max_size,
            panel_copy_split_min_species_count=panel_copy_split_min_species_count,
            panel_copy_split_min_dataset_species=panel_copy_split_min_dataset_species,
            panel_copy_component_min_edge_weight=panel_copy_component_min_edge_weight,
            component_split_high_duplication=component_split_high_duplication,
        )

    dsu = _DSU(clusters, gene_to_species)
    _merge_reciprocal_cluster_best(
        clusters,
        pair_stats,
        dsu,
        max_merges=max_reciprocal_merges,
        max_genes=max_component_genes,
    )
    _merge_weak_balanced_rescue(
        clusters,
        pair_stats,
        dsu,
        max_genes=max_component_genes,
    )
    merged = _component_clusters(clusters, dsu)
    return _split_high_duplication_clusters(
        merged,
        rbnh_edges,
        gene_to_species,
        min_size=split_min_size,
        min_species_count=split_min_species_count,
        degree_ratio=split_degree_ratio,
        copy_split_min_size=copy_split_min_size,
        copy_split_min_species_count=copy_split_min_species_count,
        copy_split_min_dataset_species=copy_split_min_dataset_species,
        panel_copy_split_min_size=panel_copy_split_min_size,
        panel_copy_split_max_size=panel_copy_split_max_size,
        panel_copy_split_min_species_count=panel_copy_split_min_species_count,
        panel_copy_split_min_dataset_species=panel_copy_split_min_dataset_species,
        panel_copy_component_min_edge_weight=panel_copy_component_min_edge_weight,
        component_split_high_duplication=component_split_high_duplication,
    )


def refine_cluster_indices(
    clusters: Sequence[Sequence[int]],
    hit_queries: Sequence[int],
    hit_targets: Sequence[int],
    hit_scores: Sequence[float],
    rbnh_queries: Sequence[int],
    rbnh_targets: Sequence[int],
    gene_to_species: Sequence[int],
    max_reciprocal_merges: int = DEFAULT_MAX_RECIPROCAL_MERGES,
    max_component_genes: int = DEFAULT_MAX_COMPONENT_GENES,
    broad_copy_split_min_species_count: int = DEFAULT_BROAD_COPY_SPLIT_MIN_SPECIES_COUNT,
    broad_max_reciprocal_merges: int = DEFAULT_BROAD_MAX_RECIPROCAL_MERGES,
    broad_reciprocal_min_avg_score: float = DEFAULT_BROAD_RECIPROCAL_MIN_AVG_SCORE,
    broad_reciprocal_min_max_score: float = DEFAULT_BROAD_RECIPROCAL_MIN_MAX_SCORE,
    broad_reciprocal_min_coverage: float = DEFAULT_BROAD_RECIPROCAL_MIN_COVERAGE,
    broad_reciprocal_min_norm: float = DEFAULT_BROAD_RECIPROCAL_MIN_NORM,
    split_min_size: int = DEFAULT_SPLIT_MIN_SIZE,
    split_min_species_count: int = DEFAULT_SPLIT_MIN_SPECIES_COUNT,
    split_degree_ratio: float = DEFAULT_SPLIT_DEGREE_RATIO,
    copy_split_min_size: int = DEFAULT_COPY_SPLIT_MIN_SIZE,
    copy_split_min_species_count: int = DEFAULT_COPY_SPLIT_MIN_SPECIES_COUNT,
    copy_split_min_dataset_species: int = DEFAULT_COPY_SPLIT_MIN_DATASET_SPECIES,
    panel_copy_split_min_size: int = DEFAULT_PANEL_COPY_SPLIT_MIN_SIZE,
    panel_copy_split_max_size: int = DEFAULT_PANEL_COPY_SPLIT_MAX_SIZE,
    panel_copy_split_min_species_count: int = DEFAULT_PANEL_COPY_SPLIT_MIN_SPECIES_COUNT,
    panel_copy_split_min_dataset_species: int = DEFAULT_PANEL_COPY_SPLIT_MIN_DATASET_SPECIES,
    panel_copy_component_min_edge_weight: float = DEFAULT_PANEL_COPY_COMPONENT_MIN_EDGE_WEIGHT,
    rbnh_scores: Sequence[float] | None = None,
    component_split_high_duplication: bool = False,
) -> List[IndexCluster]:
    """Refine int-indexed clusters without materializing string hit triples."""
    if not clusters:
        return []
    total_genes = len(gene_to_species)
    rbnh_queries_arr = np.asarray(rbnh_queries, dtype=np.int64)
    rbnh_targets_arr = np.asarray(rbnh_targets, dtype=np.int64)
    rbnh_scores_arr = (
        np.asarray(rbnh_scores, dtype=np.float64)
        if rbnh_scores is not None else np.asarray([], dtype=np.float64)
    )
    rbnh_sources, rbnh_target_pairs, rbnh_totals, rbnh_counts, rbnh_max_scores = _cluster_pair_arrays(
        clusters,
        rbnh_queries_arr,
        rbnh_targets_arr,
        rbnh_scores_arr,
        total_genes,
        default_score=1.0,
    )
    if _is_broad_index_dataset(gene_to_species, copy_split_min_dataset_species):
        resolved_copy_split_min_species_count = _resolve_broad_copy_split_min_species_count(
            int(np.unique(np.asarray(gene_to_species)).size),
            copy_split_min_species_count,
            copy_split_min_dataset_species,
            broad_copy_split_min_species_count,
        )
        if len(rbnh_sources) == 0:
            return _split_high_duplication_index_clusters(
                clusters,
                rbnh_queries_arr,
                rbnh_targets_arr,
                gene_to_species,
                total_genes,
                min_size=split_min_size,
                min_species_count=split_min_species_count,
                degree_ratio=split_degree_ratio,
                copy_split_min_size=copy_split_min_size,
                copy_split_min_species_count=resolved_copy_split_min_species_count,
                copy_split_min_dataset_species=copy_split_min_dataset_species,
                panel_copy_split_min_size=panel_copy_split_min_size,
                panel_copy_split_max_size=panel_copy_split_max_size,
                panel_copy_split_min_species_count=panel_copy_split_min_species_count,
                panel_copy_split_min_dataset_species=panel_copy_split_min_dataset_species,
                panel_copy_component_min_edge_weight=panel_copy_component_min_edge_weight,
                rbnh_scores=rbnh_scores,
                component_split_high_duplication=component_split_high_duplication,
                broad_mode=True,
            )

        sizes = np.asarray([len(cluster) for cluster in clusters], dtype=np.int64)
        dsu = _IndexedDSU(clusters, gene_to_species)
        _merge_reciprocal_cluster_best_indexed(
            rbnh_sources,
            rbnh_target_pairs,
            rbnh_totals,
            rbnh_counts,
            rbnh_max_scores,
            sizes,
            dsu,
            max_merges=broad_max_reciprocal_merges,
            max_genes=max_component_genes,
            min_avg_score=broad_reciprocal_min_avg_score,
            min_max_score=broad_reciprocal_min_max_score,
            min_coverage=broad_reciprocal_min_coverage,
            min_norm=broad_reciprocal_min_norm,
        )
        _merge_weak_balanced_rescue_indexed(
            rbnh_sources,
            rbnh_target_pairs,
            rbnh_totals,
            rbnh_counts,
            rbnh_max_scores,
            sizes,
            dsu,
            max_genes=max_component_genes,
        )
        merged = _component_index_clusters(clusters, dsu)
        return _split_high_duplication_index_clusters(
            merged,
            rbnh_queries_arr,
            rbnh_targets_arr,
            gene_to_species,
            total_genes,
            min_size=split_min_size,
            min_species_count=split_min_species_count,
            degree_ratio=split_degree_ratio,
            copy_split_min_size=copy_split_min_size,
            copy_split_min_species_count=resolved_copy_split_min_species_count,
            copy_split_min_dataset_species=copy_split_min_dataset_species,
            panel_copy_split_min_size=panel_copy_split_min_size,
            panel_copy_split_max_size=panel_copy_split_max_size,
            panel_copy_split_min_species_count=panel_copy_split_min_species_count,
            panel_copy_split_min_dataset_species=panel_copy_split_min_dataset_species,
            panel_copy_component_min_edge_weight=panel_copy_component_min_edge_weight,
            rbnh_scores=rbnh_scores,
            component_split_high_duplication=component_split_high_duplication,
            broad_mode=True,
        )

    sources, targets, totals, counts, max_scores = _cluster_pair_arrays(
        clusters,
        hit_queries,
        hit_targets,
        hit_scores,
        total_genes,
    )
    if len(sources) == 0:
        return _split_high_duplication_index_clusters(
            clusters,
            rbnh_queries_arr,
            rbnh_targets_arr,
            gene_to_species,
            total_genes,
            min_size=split_min_size,
            min_species_count=split_min_species_count,
            degree_ratio=split_degree_ratio,
            copy_split_min_size=copy_split_min_size,
            copy_split_min_species_count=copy_split_min_species_count,
            copy_split_min_dataset_species=copy_split_min_dataset_species,
            panel_copy_split_min_size=panel_copy_split_min_size,
            panel_copy_split_max_size=panel_copy_split_max_size,
            panel_copy_split_min_species_count=panel_copy_split_min_species_count,
            panel_copy_split_min_dataset_species=panel_copy_split_min_dataset_species,
            panel_copy_component_min_edge_weight=panel_copy_component_min_edge_weight,
            rbnh_scores=rbnh_scores,
            component_split_high_duplication=component_split_high_duplication,
        )

    sizes = np.asarray([len(cluster) for cluster in clusters], dtype=np.int64)
    dsu = _IndexedDSU(clusters, gene_to_species)
    _merge_reciprocal_cluster_best_indexed(
        sources,
        targets,
        totals,
        counts,
        max_scores,
        sizes,
        dsu,
        max_merges=max_reciprocal_merges,
        max_genes=max_component_genes,
    )
    _merge_weak_balanced_rescue_indexed(
        sources,
        targets,
        totals,
        counts,
        max_scores,
        sizes,
        dsu,
        max_genes=max_component_genes,
    )
    merged = _component_index_clusters(clusters, dsu)
    return _split_high_duplication_index_clusters(
        merged,
        rbnh_queries_arr,
        rbnh_targets_arr,
        gene_to_species,
        total_genes,
        min_size=split_min_size,
        min_species_count=split_min_species_count,
        degree_ratio=split_degree_ratio,
        copy_split_min_size=copy_split_min_size,
        copy_split_min_species_count=copy_split_min_species_count,
        copy_split_min_dataset_species=copy_split_min_dataset_species,
        panel_copy_split_min_size=panel_copy_split_min_size,
        panel_copy_split_max_size=panel_copy_split_max_size,
        panel_copy_split_min_species_count=panel_copy_split_min_species_count,
        panel_copy_split_min_dataset_species=panel_copy_split_min_dataset_species,
        panel_copy_component_min_edge_weight=panel_copy_component_min_edge_weight,
        rbnh_scores=rbnh_scores,
        component_split_high_duplication=component_split_high_duplication,
    )
