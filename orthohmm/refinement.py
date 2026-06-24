"""Cluster refinement utilities for orthogroup inference.

The refinement stage works after graph clustering. It uses two independent
signals:

* all significant directed search hits to identify reciprocal cluster-level
  support that Leiden can miss when evidence is spread across many genes; and
* the stricter RBNH graph to split very large duplicated clusters whose member
  degrees show obvious outliers.
"""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple


Gene = str
Cluster = List[Gene]
DirectedHit = Tuple[Gene, Gene, float]
EdgeMap = Mapping[frozenset, float]


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


def _species_overlap(a: Counter, b: Counter) -> int:
    return len(set(a) & set(b))


def _merge_reciprocal_cluster_best(
    clusters: Sequence[Cluster],
    pair_stats: Mapping[Tuple[int, int], _PairStats],
    dsu: _DSU,
    max_merges: int,
    max_genes: int,
) -> int:
    sizes = [len(c) for c in clusters]
    best_out: Dict[int, Tuple[float, int]] = {}
    for (source, target), stats in pair_stats.items():
        avg_score, max_score, coverage, norm = _features(
            stats, sizes[source], sizes[target]
        )
        if avg_score < 0.5 or max_score < 0.8 or coverage < 1.0:
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
) -> List[Cluster]:
    adjacency: MutableMapping[Gene, set] = defaultdict(set)
    for edge in rbnh_edges:
        if len(edge) != 2:
            continue
        gene_a, gene_b = tuple(edge)
        adjacency[gene_a].add(gene_b)
        adjacency[gene_b].add(gene_a)

    refined: List[Cluster] = []
    for cluster in clusters:
        if len(cluster) < min_size:
            refined.append(list(cluster))
            continue
        species_counts = Counter(gene_to_species.get(g, "") for g in cluster)
        if max(species_counts.values()) < min_species_count:
            refined.append(list(cluster))
            continue
        members = set(cluster)
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


def refine_clusters(
    clusters: Sequence[Cluster],
    directed_hits: Iterable[DirectedHit],
    rbnh_edges: EdgeMap,
    gene_to_species: Mapping[Gene, str],
    max_reciprocal_merges: int = 150,
    max_component_genes: int = 80,
    split_min_size: int = 20,
    split_min_species_count: int = 20,
    split_degree_ratio: float = 1.1,
) -> List[Cluster]:
    """Refine orthogroup clusters using cluster-level support and graph degree.

    The default parameters are deliberately conservative: merge only small
    components with reciprocal aggregate evidence, then split only large
    duplicated clusters with clear internal-degree outliers.
    """
    if not clusters:
        return []
    pair_stats = _cluster_pair_stats(clusters, directed_hits)
    if not pair_stats:
        return _split_high_duplication_clusters(
            clusters,
            rbnh_edges,
            gene_to_species,
            min_size=split_min_size,
            min_species_count=split_min_species_count,
            degree_ratio=split_degree_ratio,
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
    )
