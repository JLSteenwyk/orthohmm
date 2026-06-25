from itertools import combinations

import numpy as np

from orthohmm.refinement import refine_cluster_indices, refine_clusters


def _cluster(prefix, counts):
    genes = []
    species = {}
    for species_idx, count in enumerate(counts):
        for gene_idx in range(count):
            gene = f"{prefix}_s{species_idx}_{gene_idx}"
            genes.append(gene)
            species[gene] = f"s{species_idx}"
    return genes, species


def _hits(source, target, scores):
    return [
        (source[i % len(source)], target[i % len(target)], score)
        for i, score in enumerate(scores)
    ]


def test_weak_balanced_rescue_merges_matching_components():
    cluster_a, species_a = _cluster("a", [3, 2, 1, 2, 3, 3, 3, 3, 2])
    cluster_b, species_b = _cluster("b", [2, 4, 2, 2, 2, 2, 2, 2, 1])
    other, species_other = _cluster("other", [1])
    gene_to_species = {**species_a, **species_b, **species_other}

    weak_scores = [0.235] * 20 + [0.733]
    directed_hits = (
        _hits(cluster_a, cluster_b, weak_scores)
        + _hits(cluster_b, cluster_a, weak_scores)
    )

    refined = refine_clusters(
        [cluster_a, cluster_b, other],
        directed_hits,
        {},
        gene_to_species,
        max_reciprocal_merges=0,
    )

    refined_sets = {frozenset(cluster) for cluster in refined}
    assert frozenset(cluster_a + cluster_b) in refined_sets
    assert frozenset(other) in refined_sets


def test_reciprocal_cluster_best_merges_strong_mutual_support():
    cluster_a, species_a = _cluster("a", [1, 1, 1])
    cluster_b, species_b = _cluster("b", [1, 1, 1])
    gene_to_species = {**species_a, **species_b}
    directed_hits = (
        _hits(cluster_a, cluster_b, [1.0, 1.1, 1.2])
        + _hits(cluster_b, cluster_a, [1.0, 1.1, 1.2])
    )

    refined = refine_clusters(
        [cluster_a, cluster_b],
        directed_hits,
        {},
        gene_to_species,
        max_reciprocal_merges=1,
    )

    assert [sorted(cluster_a + cluster_b)] == [sorted(c) for c in refined]


def test_degree_split_cuts_low_degree_members_from_large_duplicate_cluster():
    cluster, gene_to_species = _cluster("dup", [22])
    core = cluster[:5]
    rbnh_edges = {frozenset(pair): 1.0 for pair in combinations(core, 2)}

    refined = refine_clusters(
        [cluster],
        [],
        rbnh_edges,
        gene_to_species,
    )

    refined_sets = {frozenset(c) for c in refined}
    assert frozenset(core) in refined_sets
    assert all(frozenset([gene]) in refined_sets for gene in cluster[5:])


def test_indexed_refinement_merges_strong_mutual_support():
    clusters = [[0, 1, 2], [3, 4, 5]]
    gene_to_species = np.asarray([0, 1, 2, 0, 1, 2], dtype=np.int32)
    hit_queries = np.asarray([0, 1, 2, 3, 4, 5], dtype=np.int32)
    hit_targets = np.asarray([3, 4, 5, 0, 1, 2], dtype=np.int32)
    hit_scores = np.asarray([1.0, 1.1, 1.2, 1.0, 1.1, 1.2], dtype=np.float32)

    refined = refine_cluster_indices(
        clusters,
        hit_queries,
        hit_targets,
        hit_scores,
        np.asarray([], dtype=np.int32),
        np.asarray([], dtype=np.int32),
        gene_to_species,
        max_reciprocal_merges=1,
    )

    assert [set(range(6))] == [set(c) for c in refined]


def test_indexed_weak_balanced_rescue_merges_matching_components():
    clusters = [
        list(range(22)),
        list(range(22, 41)),
        [41],
    ]
    species_counts_a = [3, 2, 1, 2, 3, 3, 3, 3, 2]
    species_counts_b = [2, 4, 2, 2, 2, 2, 2, 2, 1]
    gene_to_species = []
    for species_idx, count in enumerate(species_counts_a):
        gene_to_species.extend([species_idx] * count)
    for species_idx, count in enumerate(species_counts_b):
        gene_to_species.extend([species_idx] * count)
    gene_to_species.append(9)

    weak_scores = [0.235] * 20 + [0.733]
    forward_queries = [clusters[0][i % len(clusters[0])] for i in range(21)]
    forward_targets = [clusters[1][i % len(clusters[1])] for i in range(21)]
    reverse_queries = forward_targets
    reverse_targets = forward_queries

    refined = refine_cluster_indices(
        clusters,
        np.asarray(forward_queries + reverse_queries, dtype=np.int32),
        np.asarray(forward_targets + reverse_targets, dtype=np.int32),
        np.asarray(weak_scores + weak_scores, dtype=np.float32),
        np.asarray([], dtype=np.int32),
        np.asarray([], dtype=np.int32),
        np.asarray(gene_to_species, dtype=np.int32),
        max_reciprocal_merges=0,
    )

    refined_sets = {frozenset(c) for c in refined}
    assert frozenset(range(41)) in refined_sets
    assert frozenset([41]) in refined_sets


def test_indexed_refinement_splits_large_low_degree_members():
    cluster = list(range(22))
    gene_to_species = np.zeros(22, dtype=np.int32)
    edge_pairs = list(combinations(cluster[:5], 2))

    refined = refine_cluster_indices(
        [cluster],
        np.asarray([], dtype=np.int32),
        np.asarray([], dtype=np.int32),
        np.asarray([], dtype=np.float32),
        np.asarray([a for a, _b in edge_pairs], dtype=np.int32),
        np.asarray([b for _a, b in edge_pairs], dtype=np.int32),
        gene_to_species,
    )

    refined_sets = {frozenset(c) for c in refined}
    assert frozenset(cluster[:5]) in refined_sets
    assert all(frozenset([gene]) in refined_sets for gene in cluster[5:])


def test_indexed_copy_split_breaks_qfo_scale_high_copy_cluster():
    cluster = list(range(150))
    gene_to_species = [0] * 10
    for species_idx in range(1, 50):
        gene_to_species.extend([species_idx] * 2)
    for species_idx in range(1, 43):
        gene_to_species.append(species_idx)

    refined = refine_cluster_indices(
        [cluster],
        np.asarray([], dtype=np.int32),
        np.asarray([], dtype=np.int32),
        np.asarray([], dtype=np.float32),
        np.asarray([], dtype=np.int32),
        np.asarray([], dtype=np.int32),
        np.asarray(gene_to_species, dtype=np.int32),
    )

    assert len(refined) == 150
    assert all(len(c) == 1 for c in refined)
    assert {c[0] for c in refined} == set(cluster)


def test_indexed_copy_split_is_disabled_for_small_species_panels():
    cluster = list(range(150))
    gene_to_species = np.repeat(np.arange(15, dtype=np.int32), 10)

    refined = refine_cluster_indices(
        [cluster],
        np.asarray([], dtype=np.int32),
        np.asarray([], dtype=np.int32),
        np.asarray([], dtype=np.float32),
        np.asarray([], dtype=np.int32),
        np.asarray([], dtype=np.int32),
        gene_to_species,
    )

    assert [set(cluster)] == [set(c) for c in refined]
