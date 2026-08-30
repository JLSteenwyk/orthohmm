from benchmark_tools.candidate_superfamily_holdout import evaluate


def test_bounded_candidate_holdout_recovers_fragmented_families():
    one_pass = evaluate(101, 1)
    bounded = evaluate(101, 4)

    assert one_pass["pair_score"]["recall"] < 1.0
    assert bounded["pair_score"] == {
        "precision": 1.0,
        "recall": 1.0,
        "f_score": 1.0,
        "predicted_pairs": bounded["pair_score"]["expected_pairs"],
        "expected_pairs": bounded["pair_score"]["expected_pairs"],
        "true_positive_pairs": bounded["pair_score"]["expected_pairs"],
    }
