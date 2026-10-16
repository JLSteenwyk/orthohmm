from copy import deepcopy
import json

import pytest

from benchmark_tools import report_additional_hillclimb as reporter
from benchmark_tools import report_cap400_scores as recovery
from benchmark_tools.metric_frontier import ENDPOINTS
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.recover_cap400_scores import READY
from benchmark_tools.report_hillclimb import endpoint_vector
from benchmark_tools.run_membership_hillclimb import PRODUCTS, write_json
from benchmark_tools.summarize_matched_resources import LABELS

ACCOUNTING = "JobIDRaw|State|ExitCode|AllocCPUS|NodeList\n24355|COMPLETED|0:0|8|bizon\n"


def test_terminal_scoring_only_resources_required(monkeypatch):
    monkeypatch.setattr(recovery.subprocess, "check_output", lambda *args, **kwargs: ACCOUNTING)
    assert recovery.terminal({"job_id": "24355"})["AllocCPUS"] == "8"
    with pytest.raises(ValueError):
        recovery.terminal({"job_id": "24350"})


@pytest.mark.parametrize("text", [ACCOUNTING.replace("COMPLETED", "RUNNING"),
    ACCOUNTING.replace("COMPLETED", "FAILED"), ACCOUNTING.replace("0:0", "1:0"),
    ACCOUNTING.replace("|8|", "|32|"), ACCOUNTING.replace("bizon", "other"),
    ACCOUNTING + "24355|COMPLETED|0:0|8|bizon\n"])
def test_nonterminal_failed_or_wrong_scheduler_refused(monkeypatch, text):
    monkeypatch.setattr(recovery.subprocess, "check_output", lambda *args, **kwargs: text)
    with pytest.raises(ValueError):
        recovery.terminal({"job_id": "24355"})


def panel(tmp_path, monkeypatch, status):
    root = tmp_path / "producer"
    root.mkdir()
    write_json(root / "results.json", {"status": status, "job_id": "24355", "scores": {"GO": 1}})
    table = tmp_path / "old.tsv"
    table.write_text("old table\n")
    other = {"key": "other", "label": "other", "scores": {ENDPOINTS[0]: .5}, "cells": []}
    skipped = {"dataset": "QfO", "status": "not_run_failed_control", "scores": {}}
    method = {"key": "full_search:cap400", "label": "cap400",
        "scores": {ENDPOINTS[0]: .72, ENDPOINTS[-1]: .91}, "cells": [skipped]}
    parent = tmp_path / "parent.json"
    write_json(parent, {"methods": [other, method], "score_table": record(table),
                        "round_screening": {}, "limitations": []})
    failures = [{"dataset": "QfO", "status": "failed", "scores": {}, "producer_job_id": job,
                 "error": error} for job, error in (("24245", "original final Leiden crash"),
                                                    ("24350", "scoring converter import failure"))]
    monkeypatch.setattr(recovery, "load_plan", lambda *args: ({"runtime": {"files": []}}, [], deepcopy(failures)))
    monkeypatch.setattr(reporter, "runtime", lambda: {"files": []})
    monkeypatch.setattr(recovery.subprocess, "check_output", lambda *args, **kwargs: ACCOUNTING)
    return root, parent, other, method, failures


@pytest.mark.parametrize("status", ["scoring_original_cap400_native", "failed_cap400_scoring_only_recovery"])
def test_pending_or_failed_scores_not_admitted_and_both_failures_retained(tmp_path, monkeypatch, status):
    root, parent, other, old, failures = panel(tmp_path, monkeypatch, status)
    monkeypatch.setattr(recovery, "cell", lambda *args: pytest.fail("Partial scores admitted"))
    result = reporter.build(tmp_path, parent, tmp_path / "out", cap400_scoring_root=root)
    assert result["methods"][0] == other
    method = result["methods"][1]
    assert method["scores"] == old["scores"]
    assert all(failed in method["cells"] for failed in failures)
    assert old["cells"][0] in method["cells"]
    assert "NA" in (tmp_path / "out/scores.tsv").read_text()
    assert result["defaults_promoted"] is result["goal_completed"] is False


def test_completed_scoring_updates_only_qfo_with_both_failures_preserved(tmp_path, monkeypatch):
    root, parent, other, old, failures = panel(tmp_path, monkeypatch, recovery.SUCCESS)
    scores = {key: .81 for key in ENDPOINTS if key.startswith("QfO:")}
    good = {"dataset": "QfO", "status": "independent_content_readback_verified", "scores": scores,
        "products": {name: {"path": name, "bytes": 1, "sha256": "same"} for name in PRODUCTS}}
    monkeypatch.setattr(recovery, "cell", lambda *args: deepcopy(good))
    result = reporter.build(tmp_path, parent, tmp_path / "out", cap400_scoring_root=root)
    assert result["methods"][0] == other
    method = result["methods"][1]
    assert method["scores"] == {**old["scores"], **scores}
    assert all(failed in method["cells"] for failed in failures)
    assert (tmp_path / "out/snapshots/cap400_scoring_reference.json").exists()


@pytest.fixture
def completed(tmp_path, monkeypatch):
    root = tmp_path / "producer"
    root.mkdir()
    artifact = tmp_path / "native"
    artifact.write_text("native fixture\n")
    pin = record(artifact)
    raw = {key.split(": ", 1)[1]: .81 for key in ENDPOINTS if key.startswith("QfO:")}
    control = {"scores": {key: .8 for key in raw}}
    original_path = tmp_path / "fresh_failed.json"
    write_json(original_path, {"execution": {"failed_fresh": True}, "control": control})
    original_pin = record(original_path)
    downstream_path = tmp_path / "downstream_failed.json"
    write_json(downstream_path, {"graph_execution": {"cached_graph": True},
        "execution": {"native_succeeded_job_failed": True}, "control": control})
    graph_report = tmp_path / "graph.json"
    write_json(graph_report, {"execution": {"final_graph_recovery": True}})
    graph_admission = tmp_path / "graph_admission.json"
    write_json(graph_admission, {"producer": record(graph_report)})
    parent_path = tmp_path / "parent_plan.json"
    write_json(parent_path, {"admission": record(graph_admission), "control": control})
    converter = tmp_path / "executor/qfo_benchmark/og_to_pairwise.py"
    converter.parent.mkdir(parents=True)
    converter.write_text("# trusted converter fixture\n")
    plan = {"source_root": str(converter.parent.parent), "scoring_id": "fixture",
        "parent_root": str(tmp_path), "parent_plan": record(parent_path),
        "parent_report": record(downstream_path), "admission": pin,
        "scientific_revision": "original"}
    native = {"genes": 984137, "species": 78, "groups": 370680, "pairs": 5836092,
        "products": {name: pin for name in PRODUCTS}, "checked": [pin], "summary": {"native": True}}
    evidence = {"status": READY, "scores": {}, "native_inference_rerun": False,
        "producer": plan["parent_report"], "scheduler": {"State": "FAILED"},
        "scientific_revision": "original", "native": native, "checked": [pin],
        "original_failed_fresh_attempt": {"source_snapshot": original_pin}}
    admission = {key: deepcopy(evidence[key]) for key in ("status", "scores", "native_inference_rerun",
        "producer", "scheduler", "scientific_revision")}
    admission.update(audit=pin, coverage={key: native[key] for key in ("genes", "species", "groups", "pairs")})
    write_json(root / "execution_plan.json", plan)
    report = {"status": recovery.SUCCESS, "job_id": "24355", "postflight_pins_verified": True,
        "native_inference_rerun": False, "defaults_promoted": False, "goal_completed": False,
        "execution_plan": record(root / "execution_plan.json"), "native_products": native["products"],
        "native_summary": native["summary"], "native_audit": pin, "converter": record(converter),
        "coverage": admission["coverage"], "scoring": {"scores": deepcopy(raw)}, "scores": deepcopy(raw),
        "deltas": {key: .81 - .8 for key in raw}}
    write_json(root / "results.json", report)
    (root / "recovery.time.txt").write_text("\n".join(f"{label}: {value}" for label, value in
        zip(LABELS, ("0:30.00", "20", "2", "1000", "0"))) + "\n")
    monkeypatch.setattr(recovery, "read_admission", lambda *args: (deepcopy(admission), deepcopy(evidence)))
    monkeypatch.setattr(recovery, "audit", lambda *args: deepcopy(evidence))
    def raw_reader(repo, dataset, scored, actual, folder, participant):
        assert dataset == "QfO" and actual == native and folder == root
        assert participant == "ohmm_hillclimb_r1_fixture"
        return deepcopy(raw), {}, [pin]
    monkeypatch.setattr(recovery, "readback_scores", raw_reader)
    monkeypatch.setattr(recovery, "control_readback", lambda *args: (
        {"report": original_pin, "scores": endpoint_vector("QfO", control["scores"])}, [original_pin]))
    monkeypatch.setattr(recovery.subprocess, "check_output", lambda *args, **kwargs: ACCOUNTING)
    return root, plan, report


def test_completed_cell_reads_raw_scores_and_keeps_all_cost_scopes(tmp_path, completed):
    cell = recovery.cell(tmp_path, *completed)
    assert len(cell["scores"]) == 6 and cell["coverage"]["genes"] == 984137
    resources = cell["resources"]
    assert resources["failed_fresh_execution"] == {"failed_fresh": True}
    assert resources["cached_final_graph_recovery"] == {"final_graph_recovery": True}
    assert resources["completed_native_stage_with_failed_scoring_import"] == {"native_succeeded_job_failed": True}
    assert resources["scoring_only_recovery"]["gnu_time"]["exit_status"] == 0


@pytest.mark.parametrize("field,value", [("status", "running"), ("postflight_pins_verified", False),
    ("native_inference_rerun", True), ("defaults_promoted", True), ("goal_completed", True),
    ("error", "failed"), ("postflight_error", "changed")])
def test_unverified_execution_refused(tmp_path, completed, field, value):
    root, plan, report = completed
    report[field] = value
    with pytest.raises(ValueError, match="successfully verified"):
        recovery.cell(tmp_path, root, plan, report)


@pytest.mark.parametrize("field", ["native_products", "native_summary", "coverage", "native_audit", "converter"])
def test_changed_native_provenance_refused(tmp_path, completed, field):
    root, plan, report = completed
    report[field] = {}
    with pytest.raises(ValueError, match="authenticated original native"):
        recovery.cell(tmp_path, root, plan, report)


@pytest.mark.parametrize("field", ["scores", "deltas"])
def test_changed_score_or_delta_refused(tmp_path, completed, field):
    root, plan, report = completed
    report[field]["GO"] = 1
    with pytest.raises(ValueError, match="matched-control deltas"):
        recovery.cell(tmp_path, root, plan, report)


def test_incomplete_resource_record_refused(tmp_path, completed):
    root, plan, report = completed
    (root / "recovery.time.txt").write_text("Exit status: 0\n")
    with pytest.raises(ValueError, match="Incomplete GNU-time"):
        recovery.cell(tmp_path, root, plan, report)


@pytest.mark.parametrize("field,value", [("job_id", "other"), ("source_revision", "other"),
    ("scientific_revision", "other")])
def test_report_different_from_actual_launch_refused_before_content_readback(tmp_path, monkeypatch, field, value):
    root = tmp_path / "producer"
    root.mkdir()
    source = tmp_path / "executor"
    plan = {"repo": str(tmp_path), "output": str(root), "source_root": str(source),
        "source_revision": "executor", "scientific_revision": "original",
        "native_inference_rerun": False, "scoring_id": "fixture"}
    path = root / "execution_plan.json"
    write_json(path, plan)
    launch = {"job_id": "24355", "execution_plan": record(path), "source_root": str(source),
        "source_revision": "executor", "scientific_revision": "original"}
    launch_path = tmp_path / recovery.SUBMISSION
    launch_path.parent.mkdir(parents=True)
    write_json(launch_path, launch)
    report = {"job_id": "24355", "execution_plan": record(path), "source_revision": "executor",
              "scientific_revision": "original"}
    report[field] = value
    monkeypatch.setattr(recovery, "short_id", lambda *args: "fixture")
    monkeypatch.setattr(recovery, "read_admission", lambda *args: pytest.fail("Wrong launch admitted"))
    with pytest.raises(ValueError):
        recovery.load_plan(tmp_path, root, report)


def test_failed_scoring_resource_record_refused(tmp_path, completed):
    root, plan, report = completed
    path = root / "recovery.time.txt"
    path.write_text(path.read_text().replace("Exit status: 0", "Exit status: 1"))
    with pytest.raises(ValueError, match="Failed stage"):
        recovery.cell(tmp_path, root, plan, report)
