import copy
import json
from pathlib import Path

import pytest

from benchmark_tools import bind_native_orthobench_uncertainty as bridge


BASE = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


@pytest.fixture
def retained():
    return json.loads((BASE / "orthobench_factorial_results_20260916.json").read_text())


def scores(retained):
    return {cell: copy.deepcopy(score) for cell, score in retained["scores"].items()
            if cell.startswith("p0_")}


def test_all_available_contrasts_reuse_full_records_and_adjustment(retained):
    result = bridge.project(scores(retained), retained)
    assert len(result) == 12
    matched = [row for row in result if row["status"] == "native_records_matched"]
    assert len(matched) == 4
    original = {(row["on"], row["off"]): row for row in retained["comparisons"]}
    for row in matched:
        assert row["metrics"] == original[row["on"], row["off"]]["metrics"]
        assert sum(row[k] for k in ("family_f1_wins", "family_f1_ties", "family_f1_losses")) == 70
    assert all(row["metrics"] is None for row in result if row["missing_cells"])


@pytest.mark.parametrize("field", ["true_positive", "false_positive", "false_negative", "genes", "splits", "exact"])
def test_scalar_agreement_does_not_allow_changed_family_records(retained, field):
    current = scores(retained)
    row = current["p0_c1_r1"]["refog_records"][0]
    row[field] = not row[field] if field == "exact" else row[field] + 1
    with pytest.raises(ValueError, match="family records differ"):
        bridge.project(current, retained)


def test_record_order_is_not_an_analysis_difference(retained):
    current = scores(retained)
    for score in current.values():
        score["refog_records"].reverse()
    assert bridge.project(current, retained) == bridge.project(scores(retained), retained)


@pytest.mark.parametrize("mutation", ["duplicate", "missing", "unknown", "aggregate", "multiplicity", "contrast_count"])
def test_invalid_binding_refused(retained, mutation):
    current = scores(retained)
    records = current["p0_c0_r0"]["refog_records"]
    if mutation == "duplicate":
        records[-1] = copy.deepcopy(records[0])
    elif mutation == "missing":
        records.pop()
    elif mutation == "unknown":
        current["unknown"] = current.pop("p0_c0_r0")
    elif mutation == "aggregate":
        current["p0_c0_r0"]["f_score"] += .01
    elif mutation == "multiplicity":
        retained["multiplicity_endpoints"] = 12
    else:
        retained["comparisons"].pop()
    with pytest.raises(ValueError):
        bridge.project(current, retained)


def test_unavailable_cells_do_not_inherit_cached_intervals(retained):
    single = {"p0_c0_r0": scores(retained)["p0_c0_r0"]}
    result = bridge.project(single, retained)
    assert all(row["metrics"] is None and row["missing_cells"] for row in result)


def test_no_scores_is_not_a_success(retained):
    with pytest.raises(ValueError, match="absent"):
        bridge.project({}, retained)
