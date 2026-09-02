from pathlib import Path

import pytest

from benchmark_tools.analyze_phylogeny_changes import (
    external_relation_count,
    family_size_category,
    mean_pairwise_identity,
    pair_metrics,
    read_root_hogs,
    relation_sets,
    reference_nonreference_relations,
    split_summary,
)


@pytest.mark.parametrize(
    ("size", "expected"),
    [(2, "small_2_20"), (20, "small_2_20"), (21, "medium_21_50"),
     (50, "medium_21_50"), (51, "large_gt_50")],
)
def test_family_size_category_uses_fixed_boundaries(size, expected):
    assert family_size_category(size) == expected


def test_mean_pairwise_identity_ignores_gap_bearing_positions():
    assert mean_pairwise_identity(["AC-D", "ATGD"]) == pytest.approx(2 / 3)
    assert mean_pairwise_identity(["A--", "-C-"]) == 0.0


def test_relation_sets_distinguish_true_and_cross_refog_pairs():
    refog_by_gene = {"a": {0}, "b": {0}, "c": {1}, "d": {1}}
    within, cross = relation_sets([{"a", "b", "c"}, {"d"}], refog_by_gene)

    assert within == {(0, "a", "b")}
    assert cross == {("a", "c"), ("b", "c")}
    metrics = pair_metrics([{"a", "b"}, {"c", "d"}], within, cross)
    assert metrics["precision"] == pytest.approx(1 / 3)
    assert metrics["recall"] == pytest.approx(1 / 2)


def test_relation_sets_preserve_overlapping_refog_membership():
    memberships = {"a": {0, 1}, "b": {0}, "c": {1}}
    within, cross = relation_sets([{"a", "b", "c"}], memberships)

    assert within == {(0, "a", "b"), (1, "a", "c")}
    assert cross == {("b", "c")}


def test_external_relations_count_refog_contamination():
    clusters = [{"a", "b", "x"}, {"c", "y", "z"}]
    index = {
        gene: cluster_id
        for cluster_id, cluster in enumerate(clusters)
        for gene in cluster
    }

    assert external_relation_count({"a", "b", "c"}, clusters, index) == 4
    assert reference_nonreference_relations(clusters, {"a", "b", "c"}) == 4


def test_read_root_hogs_and_split_summary(tmp_path: Path):
    path = tmp_path / "root_hogs.tsv"
    path.write_text(
        "root_hog\tsource_family\tgenes\n"
        "RootHOG0\tFamily0000000\ta,b\n"
        "RootHOG1\tFamily0000000\tc\n"
        "RootHOG2\tFamily0000001\td\n"
    )
    root_hogs, sources = read_root_hogs(path)

    assert root_hogs == [{"a", "b"}, {"c"}, {"d"}]
    assert sources == ["Family0000000", "Family0000000", "Family0000001"]
    assert split_summary([{"a", "b", "c"}, {"d"}], root_hogs, sources, {"a", "c"}) == {
        "candidate_families": 2,
        "root_hogs": 3,
        "split_source_families": 1,
        "reference_bearing_split_source_families": 1,
        "cross_source_merges": 0,
        "all_candidate_genes_preserved": True,
        "candidate_genes": 4,
        "root_hog_genes": 4,
    }


def test_split_summary_rejects_cross_source_gene():
    with pytest.raises(ValueError, match="another family"):
        split_summary([{"a"}, {"b"}], [{"a", "b"}], ["Family0000000"], {"a"})
