import copy
import itertools

import pytest

from benchmark_tools.readback_viterbi_center_recurrence import timing_summary


def case():
    arms = ("original_division", "division", "recurrence")
    return {"query": "q", "target": "t", "cap": 100, "candidate_pairs": 3200,
        "median_recurrence_over_division_wall": .5, "median_recurrence_over_original_wall": .25,
        "repeats": [{"repeat": i, "order": order, "modes": {
            "original_division": {"wall_s": 4., "process_cpu_s": 16.},
            "division": {"wall_s": 2., "process_cpu_s": 8.},
            "recurrence": {"wall_s": 1., "process_cpu_s": 4.}}}
            for i, order in enumerate(itertools.permutations(arms))]}


def test_readback_recalculates_paired_wall_and_cpu_ratios_without_mutating_input():
    original = case()
    before = copy.deepcopy(original)
    result = timing_summary(original)
    assert result["median_recurrence_over_division_wall_s"] == .5
    assert result["median_recurrence_over_original_division_wall_s"] == .25
    assert result["median_recurrence_over_original_division_process_cpu_s"] == .25
    assert result["median_recurrence_over_division_process_cpu_s"] == .5
    assert original == before


@pytest.mark.parametrize("mutation", ["missing_repeat", "reordered", "wrong_number", "missing_arm", "wrong_median"])
def test_invalid_schedule_or_producer_summary_is_rejected(mutation):
    item = case()
    if mutation == "missing_repeat":
        item["repeats"].pop()
    elif mutation == "reordered":
        item["repeats"].reverse()
    elif mutation == "wrong_number":
        item["repeats"][0]["repeat"] = 10
    elif mutation == "missing_arm":
        del item["repeats"][0]["modes"]["division"]
    else:
        item["median_recurrence_over_original_wall"] = .3
    with pytest.raises(ValueError):
        timing_summary(item)


@pytest.mark.parametrize("field", ["wall_s", "process_cpu_s"])
@pytest.mark.parametrize("value", [0., -1., float("inf"), float("nan")])
def test_nonpositive_or_nonfinite_measurements_are_rejected(field, value):
    item = case()
    item["repeats"][0]["modes"]["division"][field] = value
    with pytest.raises(ValueError, match="Invalid measured time"):
        timing_summary(item)
