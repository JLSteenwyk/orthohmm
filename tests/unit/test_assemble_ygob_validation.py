from benchmark_tools.assemble_ygob_validation import enumerated_counts
from benchmark_tools.score_ygob_groups import score_groups


def test_enumeration_with_merges_splits_singletons_and_unscored_genes():
    reference = {"r1": ["a", "b", "c"], "r2": ["d", "e"], "r3": ["f"]}
    predicted = {"p1": ["a", "b", "d", "f", "unscored"], "p2": ["c", "e"]}
    expected = {"tp": 1, "fp": 6, "fn": 3}
    assert enumerated_counts(predicted, reference) == expected
    assert score_groups(predicted, reference, set("abcdef") | {"unscored"})["counts"] == expected


def test_empty_predictions_preserve_all_false_negatives():
    assert enumerated_counts({}, {"r": ["a", "b", "c"]}) == {"tp": 0, "fp": 0, "fn": 3}
