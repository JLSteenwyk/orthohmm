import copy
import itertools

import pytest

from benchmark_tools.readback_viterbi_boundary_clear import ARMS, timing_summary


def case():
    return {"query": "q", "target": "t", "cap": 100, "candidate_pairs": 2000,
        "ratios": {"boundary_over_original_full_wall_s": .25, "boundary_over_original_full_process_cpu_s": .25,
                   "boundary_over_compiled_full_wall_s": .5, "boundary_over_compiled_full_process_cpu_s": .5},
        "repeats": [{"repeat": i, "order": order, "arms": {
            "original_full": {"wall_s": 4., "process_cpu_s": 16.},
            "compiled_full": {"wall_s": 2., "process_cpu_s": 8.},
            "boundary_only": {"wall_s": 1., "process_cpu_s": 4.}}}
            for i, order in enumerate(itertools.permutations(ARMS))]}


def test_raw_paired_ratios_are_recomputed_without_mutating_input():
    item = case()
    original = copy.deepcopy(item)
    assert timing_summary(item)["ratios"] == item["ratios"]
    assert item == original


@pytest.mark.parametrize("mutation", ["missing", "order", "number", "inventory", "ratio"])
def test_incomplete_or_modified_schedule_or_ratios_are_rejected(mutation):
    item = case()
    if mutation == "missing":
        item["repeats"].pop()
    elif mutation == "order":
        item["repeats"].reverse()
    elif mutation == "number":
        item["repeats"][0]["repeat"] = 10
    elif mutation == "inventory":
        del item["repeats"][0]["arms"]["compiled_full"]
    else:
        item["ratios"]["boundary_over_original_full_wall_s"] = .3
    with pytest.raises(ValueError):
        timing_summary(item)


@pytest.mark.parametrize("key", ["wall_s", "process_cpu_s"])
@pytest.mark.parametrize("value", [0., -1., float("nan"), float("inf")])
def test_invalid_measurements_are_not_admitted(key, value):
    item = case()
    item["repeats"][0]["arms"]["boundary_only"][key] = value
    with pytest.raises(ValueError, match="measured"):
        timing_summary(item)
