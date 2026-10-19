from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import report_completed_hmm_resources as report
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.summarize_matched_resources import verbose_time
from tests.unit.test_summarize_matched_resources import sample


def cell_fixture(root, *, local=False):
    root.mkdir(parents=True)
    resources, pins = {}, []
    for phase in (("native",) if local else ("graph", "native")):
        log, clock = root / (phase + ".log"), root / (phase + ".time.txt")
        log.write_text("saved phase output\n")
        clock.write_text(sample())
        execution = {"command": ["python", "replay.py", "--cpu", "16"], "returncode": 0,
            "log": record(log), "time": record(clock), "resources": verbose_time(clock.read_text())}
        resources[phase + "_execution"] = execution
        pins += [record(log), record(clock)]
    if not local:
        resources.update(candidate_resources={"wall_s": .5, "cpu_s": .25, "peak_parent_process_rss_GiB": 8.,
            "memory_scope": "Parent lifetime maximum"}, total_incremental_wall_s=4.)
    return {"dataset": "OrthoBench", "label": "fixture", "status": "independent_content_readback_verified",
            "resources": resources, "checked": deepcopy(pins)}, pins


def test_phase_rows_never_sum_rss_or_mix_arm_and_phase_times(tmp_path):
    cell, pins = cell_fixture(tmp_path / "arm")
    rows, checked = report.cell_rows("weighting:test", cell, pins)
    assert [row["phase"] for row in rows] == ["graph", "native", "candidate_parent", "arm_total"]
    assert rows[0]["elapsed_seconds"] == rows[1]["elapsed_seconds"] == 1.25
    assert rows[0]["cpu_seconds"] == pytest.approx(2.4)
    assert rows[0]["max_process_rss_gib"] == 12345 / 1024 ** 2
    assert rows[0]["parent_lifetime_max_rss_gib"] is None
    assert rows[2]["max_process_rss_gib"] is None and rows[2]["parent_lifetime_max_rss_gib"] == 8.
    assert rows[3]["elapsed_seconds"] == 4. and rows[3]["cpu_seconds"] is None
    assert "Do not sum" in rows[3]["scope"] and "not independently" in rows[2]["evidence"].replace("_", " ")
    assert checked == pins


def test_local_native_only_does_not_fabricate_search_or_nomination_costs(tmp_path):
    cell, pins = cell_fixture(tmp_path / "arm", local=True)
    rows, checked = report.cell_rows("local:test", cell, pins, local=True)
    assert len(rows) == 1 and rows[0]["phase"] == "native"
    assert "excludes nomination" in rows[0]["scope"] and "joint alignment/tree preparation" in rows[0]["scope"]
    assert checked == pins


def test_nested_admitted_cell_inventory_is_valid_even_when_parent_omits_duplicate_pins(tmp_path):
    cell, pins = cell_fixture(tmp_path / "arm")
    assert report.cell_rows("alignment:stopped", cell, [])[1] == pins


@pytest.mark.parametrize("change", ["partial", "scope", "failed", "time", "log", "parsed", "missing_pins", "cpu"])
def test_unadmitted_incomplete_modified_or_unbound_receipts_are_rejected(tmp_path, change):
    cell, pins = cell_fixture(tmp_path / "arm")
    execution = cell["resources"]["native_execution"]
    if change == "partial": cell["status"] = "running"
    elif change == "scope": del cell["resources"]["graph_execution"]
    elif change == "failed": execution["returncode"] = 1
    elif change == "time": Path(execution["time"]["path"]).write_text(sample(elapsed="0:03.00"))
    elif change == "log": Path(execution["log"]["path"]).write_text("changed\n")
    elif change == "parsed": execution["resources"]["elapsed_seconds"] += 1
    elif change == "missing_pins": pins.clear(); cell["checked"] = []
    else: execution["command"][-1] = "32"
    with pytest.raises(ValueError): report.cell_rows("test", cell, pins)


@pytest.mark.parametrize("value", [-1., float("inf"), float("nan"), True, "1"])
@pytest.mark.parametrize("field", ["candidate", "total"])
def test_parent_metadata_requires_finite_nonnegative_numeric_values(tmp_path, field, value):
    cell, pins = cell_fixture(tmp_path / "arm")
    if field == "candidate": cell["resources"]["candidate_resources"]["peak_parent_process_rss_GiB"] = value
    else: cell["resources"]["total_incremental_wall_s"] = value
    with pytest.raises(ValueError): report.cell_rows("test", cell, pins)


def run_fixture(tmp_path, monkeypatch, change=None):
    source = tmp_path / "benchmark_tools/report_completed_hmm_resources.py"
    source.parent.mkdir()
    source.write_text("committed resource reader\n")
    proof_path = tmp_path / report.ADMISSION
    proof_path.parent.mkdir()
    proof = {"integration_readback": {}}
    for i, (kind, status) in enumerate(report.STATUS.items()):
        cell, pins = cell_fixture(tmp_path / kind, local=kind == "local")
        result = {"status": status, "accuracy_result_admitted": True, "reporter_revision": "reader" + str(i),
            "defaults_promoted": False, "goal_completed": False, "checked": pins, "cells": [cell]}
        if kind == "weighting": result["methods"] = [{"key": "weighting:test", "cells": [cell]}]
        if change == "report_status" and i == 0: result["status"] = "running"
        path = tmp_path / (kind + ".json")
        path.write_text(json.dumps(result))
        proof["integration_readback"][kind] = {"summary": {**record(path), "path": path.name},
            "reader_revision": result["reporter_revision"], "job_id": str(24500 + i)}
    blob = json.dumps(proof).encode()
    proof_path.write_bytes(blob)
    if change == "proof": proof_path.write_text("changed proof\n")
    if change == "summary": (tmp_path / "weighting.json").write_text("changed report\n")
    monkeypatch.setattr(report, "__file__", str(source))
    calls = []
    def job(*args):
        calls.append(args)
        if change == "job": raise ValueError("Reader is not successful")
        return {"fixture": list(args)}
    monkeypatch.setattr(report, "successful_job", job)
    def git(command, **kwargs):
        if command[1] == "rev-parse": return "source_revision"
        if command[-1].startswith(report.ADMISSION_REVISION + ":"): return blob
        return b"dirty\n" if change == "source" else source.read_bytes()
    monkeypatch.setattr(report.subprocess, "check_output", git)
    return tmp_path / "output", calls


def test_full_export_preserves_evidence_scopes_na_and_original_report_bindings(tmp_path, monkeypatch):
    output, calls = run_fixture(tmp_path, monkeypatch)
    result = report.run(tmp_path, output)
    assert len(result["rows"]) == 5 and len(result["checked"]) == 10
    assert calls == [("24500", 1, "16G"), ("24501", 1, "16G")]
    assert result["new_biological_inference"] is False and result["goal_completed"] is False
    assert json.loads((output / "summary.json").read_text()) == result
    table = (output / "resources.tsv").read_text()
    assert "\tNA\t" in table and "raw_GNU_time" in table and "producer_metadata" in table
    assert result["resource_table"] == record(output / "resources.tsv")
    with pytest.raises(FileExistsError): report.run(tmp_path, output)


@pytest.mark.parametrize("change", ["source", "proof", "summary", "report_status", "job"])
def test_export_rejects_modified_source_admission_or_unsuccessful_reader(tmp_path, monkeypatch, change):
    output, _ = run_fixture(tmp_path, monkeypatch, change)
    with pytest.raises(ValueError): report.run(tmp_path, output)
    assert not output.exists()
