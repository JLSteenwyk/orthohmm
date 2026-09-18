import copy
import json
from pathlib import Path

import pytest

from benchmark_tools.summarize_frontier_overhead import summarize, terminal_scheduler_rows


@pytest.fixture
def data():
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/dgx_frontier_overhead_plan_20260918.json"
    plan = json.loads(path.read_text())
    rows = [dict(index=t["index"], method=t["method"], mode=t["mode"], pair=t["pair"],
                 status="validated", native_wall_s=100 if t["mode"] == "boundary" else 102,
                 work_identity={"test_canonical_partition": "same"},
                 whole_command_screen_passed=True, flagged_intervals=None if t["mode"] == "boundary" else [])
            for t in plan["runs"]]
    return plan, rows


def test_complete_panel_arithmetic_never_admits_timings(data):
    plan, rows = data
    result = summarize(plan, rows)
    assert len(result["pairs"]) == 9
    assert result["complete_numeric_panel"] is True
    assert result["numerical_budget_met"] is True
    assert all(pair["wall_ratio_minus_one"] == pytest.approx(.02) for pair in result["pairs"])
    assert result["raw_evidence_audited"] is False
    assert result["environmental_validity_established"] is False
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("status", ["absent", "failed", "missing_evidence", "invalid_evidence"])
def test_missing_or_failed_pair_prevents_survivor_only_summary(data, status):
    plan, rows = data
    if status == "absent":
        rows.pop(0)
    else:
        rows[0]["status"] = status
    result = summarize(plan, rows)
    assert result["numerical_budget_met"] is None
    assert result["complete_numeric_panel"] is False
    assert result["pairs"][0]["wall_ratio_minus_one"] is None
    assert result["pairs"][0]["unavailable_reasons"]
    assert result["methods"][0]["available_pairs"] == 2
    assert result["methods"][0]["median"] is None
    assert result["methods"][1]["median"] == pytest.approx(.02)


@pytest.mark.parametrize("periodic,expected", [(110, False), (111, False), (95, True)])
def test_numerical_budget_and_negative_differences(data, periodic, expected):
    plan, rows = data
    for row in rows:
        if row["mode"] == "periodic":
            row["native_wall_s"] = periodic
    result = summarize(plan, rows)
    assert result["numerical_budget_met"] is expected
    assert result["methods"][0]["median"] == pytest.approx(periodic / 100 - 1)


def test_one_pair_over_budget_cannot_be_hidden_by_median(data):
    plan, rows = data
    rows[1]["native_wall_s"] = 111
    result = summarize(plan, rows)
    assert result["methods"][0]["median"] == pytest.approx(.02)
    assert result["methods"][0]["numerical_budget_met"] is False


def test_exact_ten_percent_pair_boundary_passes_without_tolerance(data):
    plan, rows = data
    rows[1]["native_wall_s"] = 110
    result = summarize(plan, rows)
    assert result["pairs"][0]["numerical_pair_budget_met"] is True
    assert result["methods"][0]["numerical_budget_met"] is True


def test_exact_five_percent_median_boundary_passes_without_tolerance(data):
    plan, rows = data
    for row in rows:
        if row["mode"] == "periodic":
            row["native_wall_s"] = 105
    assert summarize(plan, rows)["numerical_budget_met"] is True


def test_original_flags_preserved_separately_from_budget(data):
    plan, rows = data
    rows[1]["flagged_intervals"] = [2, 4]
    rows[1]["whole_command_screen_passed"] = False
    result = summarize(plan, rows)
    assert result["numerical_budget_met"] is True
    assert result["pairs"][0]["arms"]["periodic"]["flagged_intervals"] == [2, 4]
    assert result["environmental_validity_established"] is False


@pytest.mark.parametrize("fault", ["work", "duration"])
def test_non_equivalent_or_short_work_retains_descriptive_ratio(data, fault):
    plan, rows = data
    if fault == "work":
        rows[1]["work_identity"] = {"test_canonical_partition": "different"}
    else:
        rows[0]["native_wall_s"] = 50
    result = summarize(plan, rows)
    assert result["pairs"][0]["wall_ratio_minus_one"] is not None
    assert result["equivalent_work_and_duration_met"] is False


@pytest.mark.parametrize("fault", ["duplicate", "index", "bool_index", "mode", "pair", "status",
                                    "nan", "infinity", "zero", "boolean_time", "boundary_flags", "negative_flag"])
def test_malformed_audit_input_rejected(data, fault):
    plan, rows = data
    if fault == "duplicate":
        rows.append(copy.deepcopy(rows[0]))
    elif fault == "index":
        rows[0]["index"] = 18
    elif fault == "bool_index":
        rows[0]["index"] = False
    elif fault in {"mode", "pair", "status"}:
        rows[0][fault] = "wrong"
    elif fault in {"nan", "infinity", "zero", "boolean_time"}:
        rows[0]["native_wall_s"] = {"nan": float("nan"), "infinity": float("inf"), "zero": 0, "boolean_time": True}[fault]
    elif fault == "boundary_flags":
        rows[0]["flagged_intervals"] = []
    else:
        rows[1]["flagged_intervals"] = [-1]
    with pytest.raises(ValueError):
        summarize(plan, rows)


def accounting():
    return [f"21838_{i}|COMPLETED|0:0|00:10:00|20|96Gn|spark-7ff0" for i in range(18)]


def test_scheduler_failures_are_terminal_and_retained():
    rows = accounting()
    rows[4] = rows[4].replace("COMPLETED|0:0", "TIMEOUT|0:15")
    rows[5] = rows[5].replace("COMPLETED|0:0", "CANCELLED by 1000|0:15")
    result = terminal_scheduler_rows("\n".join(rows), 21838)
    assert len(result) == 18
    assert result[4]["state"] == "TIMEOUT"
    assert result[5]["state"] == "CANCELLED by 1000"


@pytest.mark.parametrize("fault", ["missing", "duplicate", "running", "pending", "fields", "exit", "wrong_array"])
def test_scheduler_gate_rejects_partial_or_ambiguous_evidence(fault):
    rows = accounting()
    if fault == "missing":
        rows.pop()
    elif fault == "duplicate":
        rows.append(rows[0])
    elif fault in {"running", "pending"}:
        rows[0] = rows[0].replace("COMPLETED", fault.upper())
    elif fault == "fields":
        rows[0] += "|extra"
    elif fault == "exit":
        rows[0] = rows[0].replace("0:0", "invalid")
    else:
        rows[0] = rows[0].replace("21838", "21839")
    with pytest.raises(ValueError):
        terminal_scheduler_rows("\n".join(rows), 21838)
