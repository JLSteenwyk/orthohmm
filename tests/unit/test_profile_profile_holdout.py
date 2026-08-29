from benchmark_tools.profile_profile_holdout import evaluate_split


def test_profile_profile_repair_rejects_close_paralog_pairs():
    result = evaluate_split(seed=1103, families=8, cpu=1)

    assert result["passed"] is True
    assert result["true_positives"] == 8
    assert result["false_positives"] == 0
