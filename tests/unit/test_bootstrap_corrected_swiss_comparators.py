import copy

import pytest

from benchmark_tools import bootstrap_corrected_swiss_comparators as module
from benchmark_tools.bootstrap_qfo_swiss_comparators import bootstrap as historical
from tests.unit.test_bootstrap_qfo_swiss_comparators import fixture as original_fixture


def fixture():
    report = original_fixture()
    report.update(status="corrected_comparison_swiss_counts_verified", publication_ready=False,
                  uncertainty_admitted=False, reference={"sha256": module.REFERENCE_SHA})
    for row in report["methods"]:
        row["status"] = "counts_verified"
    return report


def test_exact_match_to_frozen_historical_algorithm_with_same_synthetic_counts():
    expected = historical(original_fixture(), replicates=200, seed=42)
    result = module.bootstrap(fixture(), replicates=200, seed=42)
    assert result["point_estimates"] == expected["point_estimates"]
    for row, old in zip(result["comparisons"], expected["comparisons"]):
        assert row == dict(old, status="estimated")


@pytest.mark.parametrize("missing", [(6, 7), (2,), tuple(range(8))])
def test_missing_methods_are_not_imputed_or_removed_from_multiplicity(missing):
    report = fixture()
    for index in missing:
        report["methods"][index] = {"method": module.METHODS[index], "status": "not_admitted", "reason": "pending"}
    result = module.bootstrap(report, replicates=100)
    assert result["multiplicity_endpoints"] == 24
    assert len(result["comparisons"]) == 8
    for row, (candidate, reference) in zip(result["comparisons"], module.CONTRASTS):
        if candidate in missing or reference in missing:
            assert row["status"] == "not_estimable"
            assert row["metrics"] is None and row["family_differences"] is None
        else:
            assert row["status"] == "estimated"


@pytest.mark.parametrize("problem", ["order", "reference", "overlap", "bool", "negative", "stat", "nan",
                                   "aggregate", "truth", "family_order", "imputed", "missing_reason"])
def test_invalid_counts_fail_closed(problem):
    report = fixture()
    row = report["methods"][1]["families"][0]
    if problem == "order":
        report["methods"].reverse()
    elif problem == "reference":
        report["reference"]["sha256"] = "other"
    elif problem == "overlap":
        row["represented_genes"][0] = "f1g0"
    elif problem in ("bool", "negative"):
        row["counts_without_prior"]["TP"] = True if problem == "bool" else -1
    elif problem in ("stat", "nan"):
        row["statistics_with_prior"]["F1"] = .1 if problem == "stat" else float("nan")
    elif problem == "aggregate":
        report["methods"][1]["aggregate"]["PPV"] += .1
    elif problem == "truth":
        row["counts_without_prior"]["TN"] += 1
    elif problem == "family_order":
        report["methods"][1]["families"].reverse()
    elif problem == "imputed":
        report["methods"][1].update(status="not_admitted", reason="pending")
    else:
        report["methods"][1] = {"method": module.METHODS[1], "status": "not_admitted", "reason": ""}
    with pytest.raises(ValueError):
        module.validated_values(report)


def test_identical_methods_have_zero_intervals():
    report = fixture()
    for row in report["methods"]:
        row["families"] = copy.deepcopy(report["methods"][0]["families"])
        row["aggregate"] = report["methods"][0]["aggregate"].copy()
    result = module.bootstrap(report, replicates=100)
    assert all(metric["difference"] == 0 and metric["bonferroni_percentile_ci"] == [0, 0]
               and metric["family_ties"] == 18 for row in result["comparisons"] for metric in row["metrics"].values())
