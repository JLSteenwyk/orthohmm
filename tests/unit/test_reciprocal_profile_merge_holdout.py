from benchmark_tools.reciprocal_profile_merge_holdout import evaluate_split


def test_reciprocal_profile_holdout_rejects_close_paralog_merges():
    result = evaluate_split(seed=17, families=4, cpu=1)

    assert result["predicted_pairs"] == 3
    assert result["true_positives"] == 3
    assert result["false_positives"] == 0
    assert result["close_paralog_false_pairs"] == 0
