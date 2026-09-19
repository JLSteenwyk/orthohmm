import copy

import pytest

from benchmark_tools.bootstrap_corrected_swiss_comparators import bootstrap
from benchmark_tools.reproduce_corrected_swiss_comparison import verify
from tests.unit.test_bootstrap_corrected_swiss_comparators import fixture


@pytest.fixture(scope="module")
def result():
    counts = fixture()
    for index in (6, 7):
        counts["methods"][index] = {"method": counts["methods"][index]["method"],
                                    "status": "not_admitted", "reason": "pending"}
    result = bootstrap(counts)
    result.update(status="corrected_swiss_comparison_intervals_audited", scientific_inputs_admitted=True,
                  uncertainty_admitted=True, complete_panel=False, estimated_contrasts=6,
                  reconstructed_counts=counts)
    return result


def test_independent_weighted_family_sums_match_all_available_endpoints(result):
    assert verify(result) == 18


@pytest.mark.parametrize("problem", ["seed", "order", "point", "interval", "family", "wins", "missing", "complete", "count"])
def test_incorrect_arithmetic_or_admission_rejected(result, problem):
    report = copy.deepcopy(result)
    if problem == "seed":
        report["seed"] += 1
    elif problem == "order":
        report["comparisons"].reverse()
    elif problem == "point":
        report["comparisons"][0]["metrics"]["F1"]["difference"] += .001
    elif problem == "interval":
        report["comparisons"][0]["metrics"]["PPV"]["bonferroni_percentile_ci"][0] += .001
    elif problem == "family":
        report["comparisons"][0]["family_differences"][0]["F1"] += .001
    elif problem == "wins":
        report["comparisons"][0]["metrics"]["F1"]["family_wins"] += 1
    elif problem == "missing":
        report["comparisons"][5]["metrics"] = {}
    elif problem == "complete":
        report["complete_panel"] = True
    else:
        report["estimated_contrasts"] = 8
    with pytest.raises(ValueError):
        verify(report)
