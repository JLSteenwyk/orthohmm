from copy import deepcopy

import pytest

from benchmark_tools.bootstrap_orthobench import paired_bootstrap
from benchmark_tools.parameter_neighborhood_statistics import BASELINE, VARIANTS, summarize


def records():
    return [{"refog": "RefOG1", "genes": 3, "true_positive": 2, "false_positive": 1, "false_negative": 1},
            {"refog": "RefOG2", "genes": 4, "true_positive": 3, "false_positive": 0, "false_negative": 3}]


def test_complete_panel_fixed_statistics():
    inputs = {name: records() for name in (BASELINE, *VARIANTS)}
    result = summarize(inputs, {})
    assert result["replicates"] == 20000 and result["seed"] == 20260918
    assert result["planned_endpoints"] == 18
    assert result["missing_intervals"] == {}
    assert len(result["comparisons"]) == 6
    for comparison in result["comparisons"].values():
        assert comparison["family_f1_ties"] == 2
        assert all(metric["bonferroni_percentile_ci"] == [0., 0.] for metric in comparison["metrics"].values())


def test_failed_variants_keep_eighteen_endpoint_correction():
    inputs = {BASELINE: records(), VARIANTS[0]: records()}
    inputs[VARIANTS[0]][0].update(true_positive=3, false_negative=0)
    failures = {name: "terminal inference failure" for name in VARIANTS[1:]}
    result = summarize(inputs, failures)
    expected = paired_bootstrap(inputs, BASELINE, replicates=20000, seed=20260918, multiplicity_endpoints=18)
    assert result["comparisons"] == expected["comparisons"]
    assert "18 planned" in result["multiplicity"]
    assert set(result["missing_intervals"]) == set(failures)
    assert set(result["point_estimates_percent"]) == set(inputs)


def test_all_failures_report_baseline_without_fabricated_intervals():
    result = summarize({BASELINE: records()}, {name: "terminal failure" for name in VARIANTS})
    assert result["replicates"] == 0 and result["planned_replicates"] == 20000
    assert result["comparisons"] == {}
    assert len(result["missing_intervals"]) == 6
    assert set(result["point_estimates_percent"]) == {BASELINE}


@pytest.mark.parametrize("problem", ["missing", "extra", "overlap", "reason", "control", "families"])
def test_incomplete_or_inconsistent_panel_rejected(problem):
    inputs = {name: deepcopy(records()) for name in (BASELINE, *VARIANTS)}
    failures = {}
    if problem == "missing":
        inputs.pop(VARIANTS[0])
    elif problem == "extra":
        inputs["extra"] = records()
    elif problem == "overlap":
        failures[VARIANTS[0]] = "failure"
    elif problem == "reason":
        inputs.pop(VARIANTS[0])
        failures[VARIANTS[0]] = " "
    elif problem == "control":
        inputs.pop(BASELINE)
    else:
        inputs[VARIANTS[0]][0]["refog"] = "different"
    with pytest.raises(ValueError):
        summarize(inputs, failures)
