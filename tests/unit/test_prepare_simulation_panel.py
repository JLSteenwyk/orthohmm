import pytest

from benchmark_tools.prepare_simulation_panel import CONDITIONS, SEEDS, parameter_sets


def test_exact_panel_dimensions_and_pairing():
    assert len(SEEDS) == 10
    assert SEEDS[0] == 20261001 and SEEDS[-1] == 20261010
    assert len(CONDITIONS) == 4
    defaults = {mode: {"sentinel": "preserved"} for mode in "TGS"}
    for seed in SEEDS:
        sets = {c: parameter_sets(defaults, seed, c) for c in CONDITIONS}
        for first, second in (("baseline", "divergent"), ("turnover", "divergent_turnover")):
            assert sets[first]["T"] == sets[second]["T"]
            assert sets[first]["G"] == sets[second]["G"]
            difference = {key for key in sets[first]["S"] if sets[first]["S"][key] != sets[second]["S"][key]}
            assert difference == {"SCALING"}
        for parameters in sets.values():
            assert parameters["T"]["TOTAL_LINEAGES"] == "8"
            assert parameters["G"]["INITIAL_GENOME_SIZE"] == "100"
            assert parameters["S"]["SEQUENCE_SIZE"] == "300"
            assert parameters["S"]["AA_MODEL"] == "WAG"
            assert all(p["SEED"] == str(seed) for p in parameters.values())
            assert all(p["sentinel"] == "preserved" for p in parameters.values())
    assert defaults == {mode: {"sentinel": "preserved"} for mode in "TGS"}


def test_turnover_rates_and_no_transfer():
    defaults = dict.fromkeys("TGS", {})
    baseline = parameter_sets(defaults, SEEDS[0], "baseline")["G"]
    turnover = parameter_sets(defaults, SEEDS[0], "turnover")["G"]
    assert (baseline["DUPLICATION"], baseline["LOSS"]) == ("f:2", "f:1")
    assert (turnover["DUPLICATION"], turnover["LOSS"]) == ("f:10", "f:8")
    assert baseline["TRANSFER"] == turnover["TRANSFER"] == "f:0"
    assert baseline["LOSS_EXTENSION"] == "g:1"


@pytest.mark.parametrize("seed,condition", [(1, "baseline"), (20261001, "other")])
def test_out_of_protocol_input_rejected(seed, condition):
    with pytest.raises(ValueError):
        parameter_sets(dict.fromkeys("TGS", {}), seed, condition)
