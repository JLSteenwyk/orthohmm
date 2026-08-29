from benchmark_tools.profile_species_support_holdout import evaluate_species_gate


def test_species_gate_retains_true_and_removes_lineage_specific_profiles():
    result = evaluate_species_gate(seed=17, families=4, cpu=1)

    assert result["passed"] is True
    assert result["true_profiles_retained"] == 4
    assert result["lineage_specific_profiles"] == 2
    assert result["lineage_specific_profiles_retained"] == 0
