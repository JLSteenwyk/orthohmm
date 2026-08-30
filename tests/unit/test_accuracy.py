import numpy as np
import pytest

from orthohmm.accuracy import (
    build_rbnh_edges,
    build_singleton_assignment_edges,
    combine_edges,
    deduplicate_undirected_edges,
    load_accuracy_checkpoint,
    resolve_accuracy_profile,
    write_accuracy_checkpoint,
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
    assert standard.profile_expansion is False
    assert standard.leiden_seed == 0
    assert (high.kmer_k, high.max_candidates_per_query) == (4, 100)
    assert high.multipass_graph is True
    assert high.profile_expansion is True
    assert high.leiden_seed == 4


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


def test_rbnh_excludes_self_hits_before_estimating_thresholds():
    names = ["a1", "a2"]
    species = np.array([0, 0], dtype=np.int32)
    queries = np.array([0, 1, 0, 1], dtype=np.int32)
    targets = np.array([0, 1, 1, 0], dtype=np.int32)
    scores = np.array([10, 10, 5, 5], dtype=np.float64)

    edges = build_rbnh_edges(names, species, queries, targets, scores)

    assert _edge_map(edges) == {(0, 1): 5.0}


def test_rbnh_threshold_factor_can_expand_edges_without_changing_default():
    names = ["a", "b", "c", "d", "e"]
    species = np.array([0, 1, 2, 1, 0], dtype=np.int32)
    queries = np.array([0, 1, 2, 3, 2, 4, 0, 2], dtype=np.int32)
    targets = np.array([1, 0, 3, 2, 4, 2, 2, 0], dtype=np.int32)
    scores = np.array([10, 10, 10, 10, 10, 10, 4, 4], dtype=np.float64)

    default_edges = build_rbnh_edges(
        names, species, queries, targets, scores
    )
    expanded_edges = build_rbnh_edges(
        names,
        species,
        queries,
        targets,
        scores,
        threshold_factor=0.4,
    )

    assert _edge_map(default_edges) == {
        (0, 1): 10.0,
        (2, 3): 10.0,
        (2, 4): 10.0,
    }
    assert _edge_map(expanded_edges) == {
        (0, 1): 10.0,
        (0, 2): 4.0,
        (2, 3): 10.0,
        (2, 4): 10.0,
    }


@pytest.mark.parametrize("factor", [0.0, -1.0, np.inf, np.nan])
def test_rbnh_rejects_invalid_threshold_factor(factor):
    with pytest.raises(ValueError, match="finite and positive"):
        build_rbnh_edges([], [], [], [], [], threshold_factor=factor)


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


def test_accuracy_checkpoint_round_trip(tmp_path):
    checkpoint = write_accuracy_checkpoint(
        str(tmp_path),
        ["a", "b", "c"],
        [0, 1, 1],
        [0, 1],
        [1, 2],
        [3.5, 2.0],
    )

    names, species, queries, targets, scores = load_accuracy_checkpoint(
        checkpoint
    )

    assert names == ["a", "b", "c"]
    np.testing.assert_array_equal(species, [0, 1, 1])
    np.testing.assert_array_equal(queries, [0, 1])
    np.testing.assert_array_equal(targets, [1, 2])
    np.testing.assert_array_equal(scores, [3.5, 2.0])


def test_accuracy_checkpoint_rejects_corruption(tmp_path):
    checkpoint = write_accuracy_checkpoint(
        str(tmp_path), ["a", "b"], [0, 1], [0], [1], [2.0]
    )
    with (checkpoint / "hit_scores.npy").open("ab") as handle:
        handle.write(b"corrupt")

    with pytest.raises(ValueError, match="checksum mismatch"):
        load_accuracy_checkpoint(checkpoint)
