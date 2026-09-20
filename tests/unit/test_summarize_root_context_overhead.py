from copy import deepcopy

import pytest

from benchmark_tools import summarize_root_context_overhead as module
from tests.unit.test_prepare_root_context_overhead import PARENT, PROTOCOL
from benchmark_tools.prepare_root_context_overhead import build


def fixture():
    plan = build(PARENT, PROTOCOL)
    rows = [dict((key, task[key]) for key in ("index", "block", "pair", "arm", "method")) for task in plan["runs"]]
    for row in rows:
        row.update(status="validated", scientific_timings_admitted=False, native_wall_s=100.,
                   output_equivalent=True, work_identity={"method": row["method"], "groups": 3})
    return plan, rows


def test_complete_pairs_reproduce_signed_ratios_and_medians():
    plan, rows = fixture()
    expected = [-.03, .02, .01, .04, -.02, -.01, .02, .01, .03]
    for row in rows:
        if row["arm"] == "root_context":
            row["native_wall_s"] *= 1 + expected[row["pair"]]
    result = module.summarize(plan, rows, [])
    assert [p["signed_ratio"] for p in result["pairs"]] == pytest.approx(expected)
    assert [m["median_signed_ratio"] for m in result["methods"]] == pytest.approx([-.01, .03, .01])
    assert result["engineering_budget_passed"] is True
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("status", ["retained_failure", "retained_unrun", "output_mismatch", "failed_or_invalid"])
def test_partial_method_median_not_computed(status):
    plan, rows = fixture()
    rows[0]["status"] = status
    result = module.summarize(plan, rows, [])
    assert result["pairs"][0]["signed_ratio"] is None
    assert result["methods"][0]["valid_pairs"] == 2
    assert result["methods"][0]["median_signed_ratio"] is None
    assert result["methods"][0]["median_budget_passed"] is None
    assert result["engineering_budget_passed"] is None
    assert len(result["pairs"]) == 9


@pytest.mark.parametrize("value", [0, -1, True, float("inf"), float("nan")])
def test_invalid_wall_time_invalidates_pair(value):
    plan, rows = fixture()
    rows[1]["native_wall_s"] = value
    assert module.summarize(plan, rows, [])["pairs"][0]["status"] == "invalid_pair"


@pytest.mark.parametrize("fault", ["prior", "within_pair"])
def test_native_output_mismatch_invalidates_pair(fault):
    plan, rows = fixture()
    if fault == "prior":
        rows[1]["output_equivalent"] = False
    else:
        rows[1]["work_identity"]["groups"] += 1
    result = module.summarize(plan, rows, [])
    assert result["pairs"][0]["status"] == "invalid_pair"
    assert result["engineering_budget_passed"] is None


def test_panel_issue_prevents_overall_budget_conclusion():
    plan, rows = fixture()
    result = module.summarize(plan, rows, [dict(reason="overlapping observation windows")])
    assert result["engineering_budget_passed"] is None
    assert result["all_pairs_and_panel_valid"] is False
    assert all(m["median_signed_ratio"] == 0 for m in result["methods"])


@pytest.mark.parametrize("walls", [(111., 100., 100.), (106., 106., 100.)])
def test_pair_or_median_budget_failures_are_not_success(walls):
    plan, rows = fixture()
    selected = [row for row in rows if row["method"] == module.METHODS[0] and row["arm"] == "root_context"]
    for row, wall in zip(selected, walls):
        row["native_wall_s"] = wall
    assert module.summarize(plan, rows, [])["engineering_budget_passed"] is False


def test_exact_budget_boundaries_do_not_fail_from_float_roundoff():
    plan, rows = fixture()
    selected = [row for row in rows if row["method"] == module.METHODS[0] and row["arm"] == "root_context"]
    for row, wall in zip(selected, [110., 105., 100.]):
        row["native_wall_s"] = wall
    assert module.summarize(plan, rows, [])["engineering_budget_passed"] is True
    assert module.within_budget(.1 + 1e-9, .1) is False


@pytest.mark.parametrize("fault", ["order", "missing", "extra", "budget", "bool_index", "admission"])
def test_design_drift_rejected(fault):
    plan, rows = fixture()
    if fault == "order":
        rows.reverse()
    elif fault == "missing":
        rows.pop()
    elif fault == "extra":
        rows.append(deepcopy(rows[-1]))
    elif fault == "budget":
        plan["engineering_budget"]["every_pair_max"] = .2
    elif fault == "bool_index":
        rows[0]["index"] = False
    else:
        rows[0]["scientific_timings_admitted"] = 0
    with pytest.raises(ValueError):
        module.summarize(plan, rows, [])
