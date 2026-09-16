import pytest

from benchmark_tools.prepare_simulation_panel import CONDITIONS, SEEDS, VARIABLE_SEEDS, parameter_sets
from benchmark_tools.run_zombi_seeded import family_lengths
from benchmark_tools.zombi_truth import check_family_lengths


def test_new_seed_panel_preserves_original_parameters_except_seed():
    assert len(VARIABLE_SEEDS) == 10 and not set(VARIABLE_SEEDS) & set(SEEDS)
    assert VARIABLE_SEEDS[0] == 20261101 and VARIABLE_SEEDS[-1] == 20261110
    defaults = dict.fromkeys("TGS", {})
    for condition in CONDITIONS:
        original = parameter_sets(defaults, SEEDS[0], condition)
        updated = parameter_sets(defaults, VARIABLE_SEEDS[0], condition, allowed_seeds=VARIABLE_SEEDS)
        for mode in original:
            assert {k: v for k, v in original[mode].items() if k != "SEED"} == {
                k: v for k, v in updated[mode].items() if k != "SEED"}
    with pytest.raises(ValueError):
        parameter_sets(defaults, VARIABLE_SEEDS[0], "baseline")


def test_exported_lengths_match_frozen_rule_and_family():
    mapping = {"schema_version": 1, "seed": VARIABLE_SEEDS[0], "lengths": family_lengths(VARIABLE_SEEDS[0], ["1", "2"])}
    sequences = {"F1__a": ("species", "A" * mapping["lengths"]["1"])}
    assert check_family_lengths(sequences, mapping)["verified_sequences"] == 1
    with pytest.raises(ValueError, match="length differs"):
        check_family_lengths({"F1__a": ("species", "A")}, mapping)
    with pytest.raises(ValueError, match="no assigned"):
        check_family_lengths({"F3__a": ("species", "A")}, mapping)
    mapping["lengths"]["1"] += 1
    with pytest.raises(ValueError, match="mapping differs"):
        check_family_lengths(sequences, mapping)
