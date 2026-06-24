from itertools import combinations

from orthohmm.refinement import refine_clusters


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
