"""Terminal reporting keeps native outcome, accuracy and measured resources distinct."""

import copy
import json
from pathlib import Path

import pytest

from benchmark_tools import export_native_factorial_terminal_resources as current


ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def actual():
    values = {}
    for key, (path, _) in current.SOURCES.items():
        value = json.loads((ROOT / path).read_text())
        if isinstance(key, int):
            value = {k: v for k, v in value.items() if k != "evidence"}
        values[key] = value
    return values


def projected(actual):
    return current.collect(actual["baseline"], actual["qfo"], {k: actual[k] for k in current.REVIEW_SCHEMAS})


def test_complete_terminal_inventory_keeps_failed_measurements_and_failed_accuracy(actual):
    before = copy.deepcopy(actual)
    rows = projected(actual)
    assert actual == before
    assert len(rows) == 13
    assert sum(row["wall_seconds"] is not None for row in rows) == 11
    assert rows[7]["accuracy_admitted"] is True and rows[7]["wall_seconds"] is None
    assert rows[9]["job_id"] == 23891 and rows[9]["wall_seconds"] is None
    assert rows[11]["accuracy_admitted"] is False and rows[11]["wall_seconds"] == 60568.563058126
    assert rows[12]["native_outcome"] == "native_failure_retained"
    assert rows[12]["resource_observation"] == "failed_native_command_not_successful_inference"
    assert rows[12]["accuracy_admitted"] is False and rows[12]["wall_seconds"] == 64526.103655242
    assert rows[8]["wall_seconds"] == 40676.339487507
    assert all(row["new_timing_admission"] is False for row in rows)
    for index in range(7):
        assert all(rows[index][key] == actual["baseline"]["rows"][index][key] for key in current.SCOPES)


@pytest.mark.parametrize("change", [
    "inventory", "duplicate", "cell", "job", "schema", "scope", "admission", "unreviewed",
    "failed_success", "missing_resources", "mismatch", "negative", "noninteger", "missing_score",
])
def test_inconsistent_or_unsafe_projection_refused(actual, change):
    values = copy.deepcopy(actual)
    if change == "inventory": values["baseline"]["rows"].pop()
    elif change == "duplicate": values["qfo"]["rows"].append(values["qfo"]["rows"][0])
    elif change == "cell": values[8]["cell"] = "p1_c1_r1"
    elif change == "job": values[8]["job_id"] = 999
    elif change == "schema": values[11]["schema"] = "native_factorial_terminal_review_v1"
    elif change == "scope": values[8]["resource_scopes"]["wall_seconds"] = "cached_stage"
    elif change == "admission": values[8]["scientific_timings_admitted"] = True
    elif change == "unreviewed": values[8]["primary_resources_replayed"] = False
    elif change == "failed_success": values[12]["native_outputs_validated"] = True
    elif change == "missing_resources": values[7]["resources"] = dict(values[8]["resources"])
    elif change == "mismatch": values[10]["resources"]["wall_seconds"] += 1
    elif change == "negative": values[8]["resources"]["wall_seconds"] = -1
    elif change == "noninteger": values[8]["resources"]["peak_memory_bytes"] = 1.5
    elif change == "missing_score": values["qfo"]["rows"][0]["scores"]["GO"] = None
    with pytest.raises(ValueError): projected(values)


def test_actual_export_is_fresh_and_preserves_original_snapshot(tmp_path, actual):
    output = tmp_path / "export"
    result = current.run(ROOT, output)
    report = json.loads((output / "report.json").read_text())
    assert result["rows"] == 13 and result["measured_rows"] == 11
    assert report["rows"] == projected(actual)
    assert report["baseline_rows_preserved"] == actual["baseline"]["rows"]
    assert report["new_scientific_admission"] is report["scoring_repeated"] is False
    assert "NA" in (output / "resources.tsv").read_text()
    assert "unknown and potentially tool-dependent" in (output / "resources.md").read_text()
    with pytest.raises(FileExistsError): current.run(ROOT, output)
