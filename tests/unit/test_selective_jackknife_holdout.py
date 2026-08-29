from benchmark_tools.selective_jackknife_holdout import evaluate_split


def test_single_copy_jackknife_keeps_mixed_profiles_strict():
    result = evaluate_split(seed=17, families=4, cpu=1)

    assert result["passed"] is True
    assert result["selected_profiles"] == 4
    assert result["contaminated_profiles"] == 2
    assert result["selectively_relaxed_contaminated_profiles"] == 0
    assert result["evaluations"]["single_copy_jackknife"]["precision"] == 1.0
    normalized = result["evaluations"][
        "single_copy_jackknife_per_target_residue"
    ]
    assert normalized["precision"] == 1.0
    assert normalized["negative_rejection_rate"] == 1.0
    assert normalized["contaminated_profile_wins"] == 0
    assert normalized["recall"] > result["evaluations"][
        "single_copy_jackknife"
    ]["recall"]
