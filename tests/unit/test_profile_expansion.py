import numpy as np

from orthohmm.search.profile_expansion import (
    ProfileHits,
    build_cluster_profiles,
    build_strict_profile_edges,
    load_global_sequence_database,
)


def _edge_map(edges):
    return {
        (int(source), int(target)): float(weight)
        for source, target, weight in zip(
            edges.sources, edges.targets, edges.weights
        )
    }


def test_load_global_sequence_database_preserves_gene_table_order(tmp_path):
    (tmp_path / "a.fa").write_text(">a1\nACDE\n>a2\nFGHI\n")
    (tmp_path / "b.fa").write_text(">b1\nKLMN\n")

    database = load_global_sequence_database(
        str(tmp_path), ["a.fa", "b.fa"], ["b1", "a2", "a1"]
    )

    assert database.ids == ["b1", "a2", "a1"]
    assert database.lengths.tolist() == [4, 4, 4]


def test_build_cluster_profiles_skips_clusters_outside_size_range(tmp_path):
    (tmp_path / "a.fa").write_text(
        ">a1\nACDEFG\n>a2\nACDEYG\n>a3\nACDEWG\n>a4\nTTTTTT\n"
    )
    names = ["a1", "a2", "a3", "a4"]
    database = load_global_sequence_database(
        str(tmp_path), ["a.fa"], names
    )

    profiles = build_cluster_profiles(
        [[0, 1, 2], [3]],
        names,
        database,
        "BLOSUM62",
        cpu=1,
    )

    assert list(profiles) == [0]
    assert profiles[0].member_ids == ["a1", "a2", "a3"]


def test_strict_profile_edges_use_best_profile_and_one_sequence_anchor():
    names = ["a", "b", "c", "d", "candidate"]
    clusters = [[0, 1], [2, 3], [4]]
    profile_hits = ProfileHits(
        profile_cluster_ids=np.array([0, 1], dtype=np.int32),
        gene_ids=np.array([4, 4], dtype=np.int32),
        scores=np.array([10.0, 9.0]),
        evalues=np.array([1e-20, 1e-20]),
        candidate_count=2,
    )
    # candidate->a exists, so its lower score is preferred over a->candidate.
    # candidate has no outgoing hit to b, so b->candidate is the fallback and
    # becomes the strongest anchor into selected profile 0.
    hit_queries = np.array([4, 0, 1, 4], dtype=np.int32)
    hit_targets = np.array([0, 4, 4, 2], dtype=np.int32)
    hit_scores = np.array([3.0, 7.0, 5.0, 8.0])

    edges = build_strict_profile_edges(
        names,
        clusters,
        profile_hits,
        {0: 8.0, 1: 8.0},
        hit_queries,
        hit_targets,
        hit_scores,
    )

    assert _edge_map(edges) == {(1, 4): 5.0}


def test_strict_profile_edges_reject_hits_below_weakest_member():
    names = ["a", "b", "candidate"]
    hits = ProfileHits(
        profile_cluster_ids=np.array([0], dtype=np.int32),
        gene_ids=np.array([2], dtype=np.int32),
        scores=np.array([7.9]),
        evalues=np.array([1e-20]),
        candidate_count=1,
    )

    edges = build_strict_profile_edges(
        names,
        [[0, 1], [2]],
        hits,
        {0: 8.0},
        np.array([2], dtype=np.int32),
        np.array([0], dtype=np.int32),
        np.array([5.0]),
    )

    assert len(edges) == 0
