from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import report_full_qfo_recovery as reporting
from benchmark_tools import report_hillclimb as combined
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))
    return record(path)


def fixture_reports(tmp_path):
    root = tmp_path / "recovery"
    plan = {"output": str(root), "source_revision": "worker"}
    producer = {"source_revision": "worker", "execution_plan": write(root / "execution_plan.json", plan),
                "status": "running", "job_id": "24204", "arms": [
                    {"label": "control100", "status": "independent_content_readback_verified", "scores": {}},
                    {"label": "cap400", "status": "running", "scores": {"GO": .99}}]}
    write(root / "results.json", producer)
    snapshots = tmp_path / "snapshots"
    snapshots.mkdir()
    return root, snapshots, producer, plan


def admitted():
    scores = {key: .8 for key in combined.METRIC_NAMES["QfO"].values()}
    return {"status": "independent_content_readback_verified", "scores": combined.endpoint_vector("QfO", scores),
            "products": {}}


def test_missing_recovery_does_not_create_measured_zero(tmp_path):
    cells, captures = reporting.collect_recovered_cells(tmp_path, tmp_path, tmp_path)
    assert [c["method_key"] for c in cells] == ["full_search:control100", "full_search:cap400"]
    assert all(c["scores"] == {} and "missing" in c["status"] for c in cells)
    assert captures == []


def test_collector_admits_only_readback_and_keeps_immutable_snapshots(tmp_path, monkeypatch):
    root, snapshots, _, _ = fixture_reports(tmp_path)
    old = {path: path.read_bytes() for path in root.glob("*.json")}
    monkeypatch.setattr(reporting, "audit_recovered_cell", lambda *a: admitted())
    cells, captures = reporting.collect_recovered_cells(tmp_path, root, snapshots)
    assert cells[0]["scores"] == admitted()["scores"]
    assert cells[1]["scores"] == {} and cells[1]["claimed_scores"] == {"GO": .99}
    assert len(captures) == 2 and all(c["attempt_role"] == "scoring_recovery" for c in cells)
    assert all(path.read_bytes() == value for path, value in old.items())


@pytest.mark.parametrize("mutation", ["revision", "output", "plan_hash", "readback"])
def test_collector_retains_bad_admission_without_scores(tmp_path, monkeypatch, mutation):
    root, snapshots, producer, plan = fixture_reports(tmp_path)
    if mutation in ("revision", "output"):
        plan["source_revision" if mutation == "revision" else "output"] = "wrong"
        producer["execution_plan"] = write(root / "execution_plan.json", plan)
        write(root / "results.json", producer)
    elif mutation == "plan_hash":
        write(root / "execution_plan.json", {})
    def fail(*a):
        raise ValueError("Native readback mismatch")
    monkeypatch.setattr(reporting, "audit_recovered_cell", fail)
    cells, _ = reporting.collect_recovered_cells(tmp_path, root, snapshots)
    assert cells[0]["status"] == "readback_failed" and cells[0]["scores"] == {} and cells[0]["error"]


def test_variant_admission_is_gated_by_independent_control(tmp_path, monkeypatch):
    root, snapshots, producer, _ = fixture_reports(tmp_path)
    producer["arms"][0]["status"] = "failed"
    producer["arms"][1]["status"] = "independent_content_readback_verified"
    write(root / "results.json", producer)
    monkeypatch.setattr(reporting, "audit_recovered_cell", lambda *a: pytest.fail("Control failed; no variant audit"))
    cells, _ = reporting.collect_recovered_cells(tmp_path, root, snapshots)
    assert cells[1]["status"] == "readback_failed" and "matched-control" in cells[1]["error"]


def test_terminal_recovery_gap_and_original_failure_are_not_hidden(tmp_path):
    root, snapshots, producer, _ = fixture_reports(tmp_path)
    producer.update(status="failed_with_explicit_gaps", arms=[], error="Parent still RUNNING")
    write(root / "results.json", producer)
    cells, _ = reporting.collect_recovered_cells(tmp_path, root, snapshots)
    assert all(c["status"] == "not_run_terminal_recovery_gap" and c["error"] == "Parent still RUNNING" for c in cells)


def fixture_audit(tmp_path, monkeypatch):
    recovery = reporting.recovery
    monkeypatch.setattr(recovery, "short_id", lambda repo, output, index: f"short{index}")
    source, output = tmp_path / "native_source", tmp_path / "recovery"
    cases = []
    for label, cap in zip(recovery.LABELS, (100, 400)):
        folder = tmp_path / "original" / label
        cases.append({"label": label, "dataset": "QfO", "output": str(folder), "fasta": "fasta",
            "inputs": [{"path": f"species{i}.fa"} for i in range(78)], "search_max_candidates_per_query": cap,
            "scoring_id": "original_long_id", "command": recovery.native_command(source, Path("fasta"), folder,
                                                                                     ["mafft", "fasttree"], cap, 32)})
    parent = {"source_revision": recovery.REVISION, "cpu": 32, "source_root": str(source), "cases": cases}
    parent_plan = write(tmp_path / "parent_plan.json", parent)
    rows = [{"label": "control100", "status": "failed", "exit_code": 0,
        "native": {"status": "fresh_native_content_verified"}, "scores": {}, "cpu_backend_probe": {"passed": True},
        "inputs_sources_runtime_reverified": True, "error_type": "ValueError", "error": recovery.PATH_ERROR},
        {"label": "cap400", "status": "not_run_failed_control", "scores": {}, "inputs_sources_runtime_reverified": True}]
    original = {"status": "finished_with_explicit_gaps", "job_id": "24201", "source_revision": recovery.REVISION,
                "execution_plan": parent_plan, "arms": rows}
    original_identity = write(tmp_path / "original/results.json", original)
    snapshot = write(output / "parent_results_snapshot.json", original)
    plan = {"repo": str(tmp_path), "output": str(output), "parent_job": "24201", "parent_plan": parent_plan,
            "parent_results": original_identity["path"], "source_revision": "worker", "cases": [
                recovery.planned_case(tmp_path, output, case, i) for i, case in enumerate(cases)]}
    accounting = "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n24201|FAILED|1:0|01:00:00|bizon|32\n"
    resolution = {"status": "terminal_parent_native_outputs_resolved", "parent": original_identity, "snapshot": snapshot,
        "scheduler": recovery.terminal_scheduler(accounting, "24201"), "checked": [original_identity, snapshot],
        "actions": [{"label": "control100", "action": "score_existing"}, {"label": "cap400", "action": "run_missing"}]}
    producer = {"execution_plan": write(output / "execution_plan.json", plan),
                "resolved_parent_inputs": write(output / "resolved_parent_inputs.json", resolution)}
    native = {"genes": 984137, "species": 78, "groups": 1, "pairs": 1, "products": {}, "checked": [],
              "summary": {}, "metrics": {"timing": "original fresh native run"}}
    scores = {key: .8 for key in combined.METRIC_NAMES["QfO"].values()}
    row = {"label": "control100", "status": "independent_content_readback_verified", "action": "score_existing",
        "inputs_sources_runtime_reverified": True, "native_case": cases[0], "native": native,
        "scoring_output": plan["cases"][0]["scoring_output"], "scoring_id": "short0", "scores": scores,
        "scoring_recovery_wall_s": 20}
    monkeypatch.setattr(recovery, "verify_inputs", lambda *a: None)
    monkeypatch.setattr(recovery, "verified_native", lambda *a: native)
    monkeypatch.setattr(reporting.subprocess, "check_output", lambda *a, **k: accounting)
    calls = []
    def readback(repo, dataset, scored, native, folder, participant):
        calls.append((folder, participant))
        return scores, {"precision": .8}, []
    monkeypatch.setattr(reporting, "readback_scores", readback)
    return row, plan, producer, resolution, native, calls


def test_admission_uses_original_native_and_new_scoring_namespace(tmp_path, monkeypatch):
    row, plan, producer, _, _, calls = fixture_audit(tmp_path, monkeypatch)
    result = reporting.audit_recovered_cell(tmp_path, row, plan, producer)
    assert result["scores"] == admitted()["scores"]
    assert calls == [(Path(row["scoring_output"]), "ohmm_hillclimb_r1_short0")]
    assert result["native_source_revision"] == reporting.recovery.REVISION
    assert result["resources"]["scoring_recovery_wall_s"] == 20
    assert "total_full_pipeline_and_scoring_wall_s" not in result["resources"]


@pytest.mark.parametrize("mutation", ["postflight", "receipt_path", "receipt_status", "parent_path", "parent_changed",
    "scheduler", "actions", "native_case", "scoring_folder", "scoring_id", "products", "planned_case", "repository"])
def test_admission_rejects_wrong_provenance_or_handoff(tmp_path, monkeypatch, mutation):
    row, plan, producer, resolution, native, _ = fixture_audit(tmp_path, monkeypatch)
    if mutation == "postflight":
        row["inputs_sources_runtime_reverified"] = False
    elif mutation == "receipt_path":
        producer["resolved_parent_inputs"] = write(tmp_path / "foreign.json", resolution)
    elif mutation in ("receipt_status", "parent_path", "scheduler", "actions"):
        if mutation == "receipt_status":
            resolution["status"] = "unresolved"
        elif mutation == "parent_path":
            resolution["parent"] = write(tmp_path / "foreign.json", {})
        elif mutation == "scheduler":
            resolution["scheduler"]["State"] = "RUNNING"
        else:
            resolution["actions"][0]["action"] = "run_missing"
        producer["resolved_parent_inputs"] = write(Path(producer["resolved_parent_inputs"]["path"]), resolution)
    elif mutation == "parent_changed":
        write(Path(resolution["parent"]["path"]), {})
    elif mutation == "native_case":
        row["native_case"] = deepcopy(plan["cases"][0]["fresh_case"])
    elif mutation == "scoring_folder":
        row["scoring_output"] = row["native_case"]["output"]
    elif mutation == "scoring_id":
        row["scoring_id"] = "foreign"
    elif mutation == "products":
        row["native"] = {**native, "products": {"foreign": {}}}
    elif mutation == "planned_case":
        plan["cases"][0]["fresh_case"]["command"].append("--unconstrained-membership")
    else:
        plan["repo"] = str(tmp_path / "other")
    with pytest.raises(ValueError):
        reporting.audit_recovered_cell(tmp_path, row, plan, producer)


def test_combined_report_keeps_failed_original_alongside_recovery(tmp_path, monkeypatch):
    registry = tmp_path / "registry.json"
    write(registry, {"rows": [], "limitations": []})
    original = tmp_path / "full/QfO/results.json"
    identity = write(tmp_path / "original_plan.json", {"source_revision": "original"})
    write(original, {"source_revision": "original", "execution_plan": identity, "arms": [
        {"label": "control100", "status": "failed", "scores": {}, "error": reporting.recovery.PATH_ERROR},
        {"label": "cap400", "status": "not_run_failed_control", "scores": {}}]})
    recovered = {"method_key": "full_search:control100", "dataset": "QfO", "attempt_role": "scoring_recovery",
                 **admitted()}
    monkeypatch.setattr(reporting, "collect_recovered_cells", lambda *a: ([recovered], []))
    result = combined.build_report(tmp_path, registry, tmp_path / "missing_membership", tmp_path / "missing_supplement",
                                  original.parent.parent, tmp_path / "comparison", full_qfo_recovery_root=tmp_path / "recovery")
    method = next(m for m in result["methods"] if m["key"] == "full_search:control100")
    cells = [c for c in method["cells"] if c["dataset"] == "QfO"]
    assert len(cells) == 2 and cells[0]["status"] == "failed" and cells[0]["error"] == reporting.recovery.PATH_ERROR
    assert cells[1]["attempt_role"] == "scoring_recovery" and method["scores"] == admitted()["scores"]
