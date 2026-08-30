import numpy as np

from orthohmm.search.profile_expansion import (
    ProfileHits,
    build_cluster_profiles,
    build_direct_profile_fallback_edges,
    build_reciprocal_profile_edges,
    build_strict_profile_edges,
    compute_profile_self_thresholds,
    load_global_sequence_database,
    select_single_copy_profile_ids,
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


def test_build_cluster_profiles_requires_cross_species_support(tmp_path):
    (tmp_path / "a.fa").write_text(
        ">a1\nACDEFG\n>a2\nACDEYG\n>a3\nACDEWG\n"
        ">b1\nKLMNPQ\n>b2\nKLMNAQ\n>b3\nKLMNWQ\n"
    )
    names = ["a1", "a2", "a3", "b1", "b2", "b3"]
    database = load_global_sequence_database(
        str(tmp_path), ["a.fa"], names
    )

    profiles = build_cluster_profiles(
        [[0, 1, 2], [3, 4, 5]],
        names,
        database,
        "BLOSUM62",
        cpu=1,
        gene_to_species=[0, 0, 0, 0, 1, 2],
        min_species_count=3,
    )

    assert list(profiles) == [1]


def test_profile_species_gate_requires_species_mapping(tmp_path):
    (tmp_path / "a.fa").write_text(
        ">a1\nACDEFG\n>a2\nACDEYG\n>a3\nACDEWG\n"
    )
    names = ["a1", "a2", "a3"]
    database = load_global_sequence_database(
        str(tmp_path), ["a.fa"], names
    )

    with np.testing.assert_raises_regex(
        ValueError, "gene_to_species is required"
    ):
        build_cluster_profiles(
            [[0, 1, 2]],
            names,
            database,
            "BLOSUM62",
            cpu=1,
            min_species_count=3,
        )


def test_parallel_profile_build_matches_serial(tmp_path):
    (tmp_path / "a.fa").write_text(
        ">a1\nACDEFG\n>a2\nACDEYG\n>a3\nACDEWG\n"
        ">b1\nKLMNPQ\n>b2\nKLMNAQ\n>b3\nKLMNWQ\n"
    )
    names = ["a1", "a2", "a3", "b1", "b2", "b3"]
    database = load_global_sequence_database(
        str(tmp_path), ["a.fa"], names
    )
    clusters = [[0, 1, 2], [3, 4, 5]]

    serial = build_cluster_profiles(
        clusters, names, database, "BLOSUM62", cpu=1
    )
    parallel = build_cluster_profiles(
        clusters, names, database, "BLOSUM62", cpu=2
    )

    assert sorted(parallel) == sorted(serial)
    for cluster_id in serial:
        assert parallel[cluster_id].member_ids == serial[cluster_id].member_ids
        np.testing.assert_array_equal(
            parallel[cluster_id].match_emissions,
            serial[cluster_id].match_emissions,
        )


def test_weakest_member_jackknife_never_tightens_threshold(tmp_path):
    (tmp_path / "a.fa").write_text(
        ">a1\nACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY\n"
        ">a2\nACDEFGHIKLMNPQRSTVWYACDEYGHIKLMNPQRSTVWY\n"
        ">a3\nACDEYGHIKLMAPQRSTVWYFCDEFGHIKLMNPQKSTVWY\n"
    )
    names = ["a1", "a2", "a3"]
    database = load_global_sequence_database(
        str(tmp_path), ["a.fa"], names
    )
    profiles = build_cluster_profiles(
        [[0, 1, 2]], names, database, "BLOSUM62", cpu=1
    )

    in_sample = compute_profile_self_thresholds(
        profiles, names, database, cpu=1
    )
    calibrated = compute_profile_self_thresholds(
        profiles,
        names,
        database,
        cpu=1,
        matrix_name="BLOSUM62",
        calibrate_weakest_member=True,
    )

    assert calibrated[0] < in_sample[0]


def test_profile_thresholds_can_be_expressed_per_target_residue(tmp_path):
    (tmp_path / "a.fa").write_text(
        ">a1\nACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY\n"
        ">a2\nACDEFGHIKLMNPQRSTVWYACDEYGHIKLMNPQRSTVWY\n"
        ">a3\nACDEYGHIKLMAPQRSTVWYFCDEFGHIKLMNPQKSTVWY\n"
    )
    names = ["a1", "a2", "a3"]
    database = load_global_sequence_database(
        str(tmp_path), ["a.fa"], names
    )
    profiles = build_cluster_profiles(
        [[0, 1, 2]], names, database, "BLOSUM62", cpu=1
    )

    raw = compute_profile_self_thresholds(
        profiles, names, database, cpu=1
    )
    per_residue = compute_profile_self_thresholds(
        profiles,
        names,
        database,
        cpu=1,
        score_per_target_residue=True,
    )

    assert per_residue[0] == raw[0] / 40.0


def test_jackknife_can_be_limited_to_single_copy_profiles(tmp_path):
    (tmp_path / "a.fa").write_text(
        ">a1\nACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY\n"
        ">a2\nACDEFGHIKLMNPQRSTVWYACDEYGHIKLMNPQRSTVWY\n"
        ">a3\nACDEYGHIKLMAPQRSTVWYFCDEFGHIKLMNPQKSTVWY\n"
        ">b1\nKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHI\n"
        ">b2\nKLMNPQRSTVWYACDEYGHIKLMNPQRSTVWYACDEFGHI\n"
        ">b3\nKLMAPQRSTVWYFCDEFGHIKLMNPQKSTVWYACDEFGHI\n"
    )
    names = ["a1", "a2", "a3", "b1", "b2", "b3"]
    clusters = [[0, 1, 2], [3, 4, 5]]
    database = load_global_sequence_database(str(tmp_path), ["a.fa"], names)
    profiles = build_cluster_profiles(
        clusters, names, database, "BLOSUM62", cpu=1
    )
    selected = select_single_copy_profile_ids(
        profiles, clusters, [0, 1, 2, 0, 0, 1]
    )
    strict = compute_profile_self_thresholds(
        profiles, names, database, cpu=1
    )
    calibrated = compute_profile_self_thresholds(
        profiles,
        names,
        database,
        cpu=1,
        matrix_name="BLOSUM62",
        calibrate_weakest_member=True,
        calibration_profile_ids=selected,
    )

    assert selected == {0}
    assert calibrated[0] < strict[0]
    assert calibrated[1] == strict[1]


def test_pair_profile_uses_scaled_bidirectional_holdout_threshold(tmp_path):
    (tmp_path / "a.fa").write_text(
        ">a1\nACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY\n"
        ">a2\nACDEFGHIKLMNPQRSTVWYACDEYGHIKLMNPQRSTVWY\n"
    )
    names = ["a1", "a2"]
    clusters = [[0, 1]]
    database = load_global_sequence_database(str(tmp_path), ["a.fa"], names)
    profiles = build_cluster_profiles(
        clusters,
        names,
        database,
        "BLOSUM62",
        cpu=1,
        min_cluster_size=2,
    )
    selected = select_single_copy_profile_ids(profiles, clusters, [0, 1])
    heldout = compute_profile_self_thresholds(
        profiles,
        names,
        database,
        cpu=1,
        matrix_name="BLOSUM62",
        calibrate_weakest_member=True,
        calibration_profile_ids=selected,
    )
    relaxed = compute_profile_self_thresholds(
        profiles,
        names,
        database,
        cpu=1,
        matrix_name="BLOSUM62",
        calibrate_weakest_member=True,
        calibration_profile_ids=selected,
        pair_profile_threshold_ratio=0.7,
    )

    assert selected == {0}
    assert relaxed[0] == heldout[0] * 0.7


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


def test_strict_profile_edges_support_per_residue_thresholds():
    names = ["a", "b", "candidate"]
    hits = ProfileHits(
        profile_cluster_ids=np.array([0], dtype=np.int32),
        gene_ids=np.array([2], dtype=np.int32),
        scores=np.array([7.9]),
        evalues=np.array([1e-20]),
        candidate_count=1,
        scores_per_target_residue=np.array([3.95]),
    )

    edges = build_strict_profile_edges(
        names,
        [[0, 1], [2]],
        hits,
        {0: 3.8},
        np.array([2], dtype=np.int32),
        np.array([0], dtype=np.int32),
        np.array([5.0]),
        score_per_target_residue=True,
    )

    assert _edge_map(edges) == {(0, 2): 5.0}


def test_direct_profile_fallback_anchors_unconnected_calibrated_winner():
    names = ["a", "b", "candidate"]
    hits = ProfileHits(
        profile_cluster_ids=np.array([0], dtype=np.int32),
        gene_ids=np.array([2], dtype=np.int32),
        scores=np.array([10.0]),
        evalues=np.array([1e-20]),
        candidate_count=1,
    )

    edges = build_direct_profile_fallback_edges(
        names,
        [[0, 1], [2]],
        hits,
        {0: 8.0},
        np.array([], dtype=np.int32),
        np.array([], dtype=np.int32),
        np.array([], dtype=np.float64),
        {0},
    )

    assert _edge_map(edges) == {(0, 2): 1.25}


def test_direct_profile_fallback_skips_existing_sequence_anchor():
    names = ["a", "b", "candidate"]
    hits = ProfileHits(
        profile_cluster_ids=np.array([0], dtype=np.int32),
        gene_ids=np.array([2], dtype=np.int32),
        scores=np.array([10.0]),
        evalues=np.array([1e-20]),
        candidate_count=1,
    )

    edges = build_direct_profile_fallback_edges(
        names,
        [[0, 1], [2]],
        hits,
        {0: 8.0},
        np.array([2], dtype=np.int32),
        np.array([1], dtype=np.int32),
        np.array([0.5]),
        {0},
    )

    assert len(edges) == 0


def test_direct_profile_fallback_requires_allowed_profile():
    names = ["a", "b", "candidate"]
    hits = ProfileHits(
        profile_cluster_ids=np.array([0], dtype=np.int32),
        gene_ids=np.array([2], dtype=np.int32),
        scores=np.array([10.0]),
        evalues=np.array([1e-20]),
        candidate_count=1,
    )

    edges = build_direct_profile_fallback_edges(
        names,
        [[0, 1], [2]],
        hits,
        {0: 8.0},
        np.array([], dtype=np.int32),
        np.array([], dtype=np.int32),
        np.array([], dtype=np.float64),
        set(),
    )

    assert len(edges) == 0


def test_reciprocal_profile_edges_use_bidirectional_calibrated_support():
    names = ["a0", "a1", "b0", "b1"]
    hits = ProfileHits(
        profile_cluster_ids=np.array([1, 1, 0, 0], dtype=np.int32),
        gene_ids=np.array([0, 1, 2, 3], dtype=np.int32),
        scores=np.array([8.0, 7.5, 9.0, 8.0]),
        evalues=np.full(4, 1e-20),
        candidate_count=4,
    )

    edges, pair_count = build_reciprocal_profile_edges(
        names,
        [[0, 1], [2, 3]],
        hits,
        {0: 10.0, 1: 10.0},
        [0, 1, 2, 3],
        threshold_ratio=0.7,
        min_support=2,
    )

    assert pair_count == 1
    assert _edge_map(edges) == {
        (0, 2): 0.8,
        (0, 3): 0.8,
        (1, 2): 0.75,
        (1, 3): 0.75,
    }


def test_reciprocal_profile_edges_reject_species_overlap():
    names = ["a0", "a1", "b0", "b1"]
    hits = ProfileHits(
        profile_cluster_ids=np.array([1, 1, 0, 0], dtype=np.int32),
        gene_ids=np.array([0, 1, 2, 3], dtype=np.int32),
        scores=np.full(4, 8.0),
        evalues=np.full(4, 1e-20),
        candidate_count=4,
    )

    edges, pair_count = build_reciprocal_profile_edges(
        names,
        [[0, 1], [2, 3]],
        hits,
        {0: 10.0, 1: 10.0},
        [0, 1, 2, 1],
    )

    assert pair_count == 0
    assert len(edges) == 0


def test_reciprocal_profile_edges_require_minimum_support_each_way():
    names = ["a0", "a1", "b0", "b1"]
    hits = ProfileHits(
        profile_cluster_ids=np.array([1, 1, 0], dtype=np.int32),
        gene_ids=np.array([0, 1, 2], dtype=np.int32),
        scores=np.full(3, 8.0),
        evalues=np.full(3, 1e-20),
        candidate_count=3,
    )

    edges, pair_count = build_reciprocal_profile_edges(
        names,
        [[0, 1], [2, 3]],
        hits,
        {0: 10.0, 1: 10.0},
        [0, 1, 2, 3],
        min_support=2,
    )

    assert pair_count == 0
    assert len(edges) == 0
