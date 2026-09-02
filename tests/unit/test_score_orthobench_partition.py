import pytest

from benchmark_tools.score_orthobench_partition import score_partition


def test_score_partition_matches_equal_group_weighting_and_uncertainty():
    predictions = [{"a", "b", "x"}, {"c", "d"}]
    references = {"RefOG001.txt": {"a", "b"}, "RefOG002.txt": {"c", "d"}}
    uncertain = {"RefOG001.txt": {"x"}, "RefOG002.txt": set()}

    score = score_partition(predictions, references, uncertain)

    assert score["f_score"] == pytest.approx(100.0)
    assert score["precision"] == pytest.approx(100.0)
    assert score["recall"] == pytest.approx(100.0)
    assert score["exact_refogs"] == 2


def test_score_partition_penalizes_splits_and_extra_members():
    score = score_partition(
        [{"a", "x"}, {"b"}],
        {"RefOG001.txt": {"a", "b"}},
        {},
    )

    assert score["weighted_true_positive"] == 0
    assert score["weighted_false_positive"] == 1
    assert score["weighted_false_negative"] == 1
    assert score["exact_refogs"] == 0
