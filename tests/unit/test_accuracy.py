import numpy as np
import pytest

from orthohmm.accuracy import (
    build_rbnh_edges,
    build_singleton_assignment_edges,
    combine_edges,
    deduplicate_undirected_edges,
    resolve_accuracy_profile,
)


def _edge_map(edges):
    return {
        (int(source), int(target)): float(weight)
        for source, target, weight in zip(
            edges.sources, edges.targets, edges.weights
        )
    }


def test_accuracy_profiles_keep_standard_and_high_search_distinct():
    standard = resolve_accuracy_profile("standard")
    high = resolve_accuracy_profile("high_sensitivity")

    assert (standard.kmer_k, standard.max_candidates_per_query) == (5, 0)
    assert standard.multipass_graph is False
    assert (high.kmer_k, high.max_candidates_per_query) == (4, 100)
    assert high.multipass_graph is True


def test_unknown_accuracy_profile_is_rejected():
    with pytest.raises(ValueError, match="Unknown accuracy profile"):
        resolve_accuracy_profile("maximum")


def test_rbnh_uses_weakest_reciprocal_threshold_from_either_endpoint():
    names = ["a1", "a2", "b1", "b2", "c1"]
    species = np.array([0, 0, 1, 1, 2], dtype=np.int32)
    queries = np.array([0, 0, 2, 1, 3, 0, 4, 0, 3], dtype=np.int32)
    targets = np.array([2, 3, 0, 3, 1, 4, 0, 1, 0], dtype=np.int32)
    scores = np.array([5, 3, 4, 6, 6, 2, 2, 9, 1], dtype=np.float64)

    edges = build_rbnh_edges(names, species, queries, targets, scores)

    assert _edge_map(edges) == {
        (0, 1): 9.0,
        (0, 2): 5.0,
        (0, 3): 3.0,
        (0, 4): 2.0,
        (1, 3): 6.0,
    }


def test_rbnh_preserves_first_best_target_on_exact_tie():
    names = ["a", "b1", "b2"]
    species = np.array([0, 1, 1], dtype=np.int32)
    queries = np.array([0, 0, 1, 2], dtype=np.int32)
    targets = np.array([1, 2, 0, 0], dtype=np.int32)
    scores = np.array([5, 5, 5, 5], dtype=np.float64)

    edges = build_rbnh_edges(names, species, queries, targets, scores)

    assert _edge_map(edges) == {(0, 1): 5.0, (0, 2): 5.0}


def test_singleton_assignment_accepts_cluster_zero_and_all_matching_hits():
    names = ["a", "b", "singleton", "c", "d"]
    clusters = [[0, 1], [2], [3, 4]]
    queries = np.array([2, 2, 2], dtype=np.int32)
    targets = np.array([0, 1, 3], dtype=np.int32)
    scores = np.array([5, 4, 3], dtype=np.float64)

    edges = build_singleton_assignment_edges(
        names, clusters, queries, targets, scores
    )

    assert _edge_map(edges) == {(0, 2): 5.0, (1, 2): 4.0}


def test_singleton_assignment_selects_one_best_cluster():
    names = ["a", "b", "singleton", "c", "d"]
    clusters = [[0, 1], [2], [3, 4]]
    queries = np.array([2, 2, 2], dtype=np.int32)
    targets = np.array([0, 3, 4], dtype=np.int32)
    scores = np.array([3, 5, 2], dtype=np.float64)

    edges = build_singleton_assignment_edges(
        names, clusters, queries, targets, scores
    )

    assert _edge_map(edges) == {(2, 3): 5.0, (2, 4): 2.0}


def test_edge_combination_keeps_maximum_duplicate_weight():
    names = ["a", "b", "c"]
    first = deduplicate_undirected_edges(names, [0, 1], [1, 2], [2, 3])
    second = deduplicate_undirected_edges(names, [1, 0], [0, 2], [4, 1])

    merged = combine_edges(first, second)

    assert _edge_map(merged) == {(0, 1): 4.0, (0, 2): 1.0, (1, 2): 3.0}
