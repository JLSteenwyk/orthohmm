from copy import deepcopy
import csv
import json

import pytest

from benchmark_tools import report_additional_hillclimb as reporter
from benchmark_tools import report_protected_union_pilot as pilot
from benchmark_tools.metric_frontier import ENDPOINTS
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.report_hillclimb import endpoint_vector
from benchmark_tools.run_membership_hillclimb import PRODUCTS, write_json
from benchmark_tools.summarize_matched_resources import LABELS

ACCOUNTING = "JobIDRaw|State|ExitCode|AllocCPUS|NodeList\n24363|COMPLETED|0:0|16|bizon\n"


@pytest.mark.parametrize("text", [ACCOUNTING.replace("COMPLETED", "RUNNING"),
    ACCOUNTING.replace("COMPLETED", "FAILED"), ACCOUNTING.replace("0:0", "1:0"),
    ACCOUNTING.replace("|16|", "|8|"), ACCOUNTING.replace("bizon", "other"),
    ACCOUNTING + "24363|COMPLETED|0:0|16|bizon\n", ""])
def test_terminal_success_with_exact_resources_is_required(monkeypatch, text):
    monkeypatch.setattr(pilot.subprocess, "check_output", lambda *args, **kwargs: text)
    with pytest.raises(ValueError):
        pilot.terminal({"job_id": "24363"})


def test_wrong_job_refused_and_terminal_success_accepted(monkeypatch):
    monkeypatch.setattr(pilot.subprocess, "check_output", lambda *args, **kwargs: ACCOUNTING)
    with pytest.raises(ValueError):
        pilot.terminal({"job_id": "24362"})
    assert pilot.terminal({"job_id": "24363"})["AllocCPUS"] == "16"


def panel(tmp_path, monkeypatch, status):
    root = tmp_path / "producer"
    root.mkdir()
    write_json(root / "results.json", {"job_id": "24363", "status": status,
        "scores": {"OrthoBench F1": 1}, **({"error": "original failure"} if status == "failed" else {})})
    table = tmp_path / "old.tsv"
    table.write_text("old table\n")
    control = {"scores": {"OrthoBench F1": .74}, "details": {"precision": .8, "recall": .69},
        "coverage": {"genes": 251378, "species": 12}, "products": {}, "source_revision": "control"}
    canonical = {"key": "canonical_collection:canonical", "label": "canonical", "scores": {
        **{key: .8 for key in ENDPOINTS}, **endpoint_vector("OrthoBench", control["scores"])},
        "cells": [{"dataset": "OrthoBench", "status": "independent_content_readback_verified",
                   **control, "scores": endpoint_vector("OrthoBench", control["scores"])}]}
    others = [{"key": "other" + str(i), "label": "other" + str(i),
               "scores": {ENDPOINTS[0]: .5}, "cells": []} for i in range(26)]
    methods = [canonical, *others]
    parent = tmp_path / "parent.json"
    write_json(parent, {"methods": methods, "score_table": record(table), "round_screening": {}, "limitations": []})
    monkeypatch.setattr(pilot, "load_plan", lambda *args: ({"runtime": {}}, [], []))
    monkeypatch.setattr(reporter, "runtime", lambda: {})
    monkeypatch.setattr(pilot.subprocess, "check_output", lambda *args, **kwargs: ACCOUNTING)
    good = {"dataset": "OrthoBench", "status": "independent_content_readback_verified",
        "scores": endpoint_vector("OrthoBench", {"OrthoBench F1": .75}), "matched_control": control,
        "products": {}, "gate": {"passes_development_gate": True}}
    return root, parent, methods, good


@pytest.mark.parametrize("status", ["fresh_full_orthobench_native_running", "native_verified_scoring_pending", "failed"])
def test_pending_and_failed_never_supply_scores_or_modify_existing_27(tmp_path, monkeypatch, status):
    root, parent, old, _ = panel(tmp_path, monkeypatch, status)
    monkeypatch.setattr(pilot, "cell", lambda *args: pytest.fail("Partial scores admitted"))
    result = reporter.build(tmp_path, parent, tmp_path / "out", protected_pilot_root=root)
    assert result["methods"][:27] == old and len(result["methods"]) == 28
    method = result["methods"][-1]
    assert method["key"] == pilot.KEY and method["scores"] == {}
    assert method["complete_eight_endpoints"] is False
    if status == "failed":
        assert method["cells"][0]["error"] == "original failure"
    rows = list(csv.DictReader((tmp_path / "out/scores.tsv").open(), delimiter="\t"))
    assert all(rows[-1][key] == "NA" for key in ENDPOINTS)
    assert result["defaults_promoted"] is result["goal_completed"] is False


@pytest.mark.parametrize("passes", [True, False])
def test_completed_adds_only_ob_and_does_not_rewrite_control_or_claim_full_vector(tmp_path, monkeypatch, passes):
    root, parent, old, good = panel(tmp_path, monkeypatch, pilot.SUCCESS)
    good["gate"]["passes_development_gate"] = passes
    monkeypatch.setattr(pilot, "cell", lambda *args: deepcopy(good))
    result = reporter.build(tmp_path, parent, tmp_path / "out", protected_pilot_root=root)
    assert result["methods"][:27] == old
    method = result["methods"][-1]
    assert method["scores"] == good["scores"] and method["complete_eight_endpoints"] is False
    assert method["cells"][0]["gate"]["passes_development_gate"] is passes
    assert (tmp_path / "out/snapshots/protected_pilot_reference.json").exists()
    screen = result["round_screening"]["protected_union"]["comparisons"][pilot.KEY]
    assert screen["complete_eight_endpoint_contrast"] is screen["passes_point_estimate_screen"] is False


def test_terminal_failed_success_claim_refused_even_before_admission(tmp_path, monkeypatch):
    root, _, _, _ = panel(tmp_path, monkeypatch, pilot.SUCCESS)
    monkeypatch.setattr(pilot.subprocess, "check_output", lambda *args, **kwargs: ACCOUNTING.replace("COMPLETED", "FAILED"))
    captures = tmp_path / "snapshots"
    captures.mkdir()
    with pytest.raises(ValueError):
        pilot.capture(root / "results.json", captures)


def test_failed_attempt_survives_later_verified_pilot(tmp_path, monkeypatch):
    root, parent, old, good = panel(tmp_path, monkeypatch, "failed")
    failed = reporter.build(tmp_path, parent, tmp_path / "failed", protected_pilot_root=root)
    attempt = deepcopy(failed["methods"][-1]["cells"][0])
    write_json(root / "results.json", {"job_id": "24363", "status": pilot.SUCCESS})
    monkeypatch.setattr(pilot, "cell", lambda *args: deepcopy(good))
    result = reporter.build(tmp_path, tmp_path / "failed/summary.json", tmp_path / "verified", protected_pilot_root=root)
    assert result["methods"][:27] == old
    assert attempt in result["methods"][-1]["cells"]
    assert result["methods"][-1]["scores"] == good["scores"]


@pytest.mark.parametrize("change", ["missing", "duplicate", "score", "details", "coverage", "products", "source_revision"])
def test_existing_canonical_record_must_match_without_replacement(tmp_path, monkeypatch, change):
    _, _, methods, good = panel(tmp_path, monkeypatch, pilot.SUCCESS)
    if change == "missing":
        methods.pop(0)
    elif change == "duplicate":
        methods.append(deepcopy(methods[0]))
    elif change == "score":
        methods[0]["scores"][ENDPOINTS[0]] = .1
    else:
        methods[0]["cells"][0][change] = "changed"
    with pytest.raises(ValueError, match="Existing canonical control"):
        pilot.require_existing_control(methods, good)


@pytest.fixture
def completed(tmp_path, monkeypatch):
    root = tmp_path / "producer"
    folder = root / "native"
    folder.mkdir(parents=True)
    artifact = tmp_path / "artifact"
    artifact.write_text("fixture\n")
    pin = record(artifact)
    raw = {"OrthoBench F1": .75}
    details = {"precision": .81, "recall": .70, "family_details": {"precision": 81, "recall": 70}}
    control = {"scores": {"OrthoBench F1": .74}, "details": {"precision": .80, "recall": .69}, "checked": [pin]}
    plan = {"source_revision": pilot.REVISION, "case": {"output": str(folder), "command": ["fresh", "command"]},
        "control": control}
    write_json(root / "execution_plan.json", plan)
    native = {"genes": 251378, "species": 12, "groups": 60000, "pairs": 1000000,
        "prediction": str(artifact), "products": {name: pin for name in PRODUCTS}, "summary": {"native": True}}
    routes = {"actual_initial_and_profile_cpu_kernels": True}
    (folder / "native.time.txt").write_text("\n".join(f"{label}: {value}" for label, value in
        zip(LABELS, ("0:30.00", "20", "2", "1000", "0"))) + "\n")
    (folder / "native.log").write_text("completed native command\n")
    execution = {"command": plan["case"]["command"], "returncode": 0,
        "time": record(folder / "native.time.txt"), "log": record(folder / "native.log"),
        "resources": pilot.verbose_time((folder / "native.time.txt").read_text())}
    report = {"status": pilot.SUCCESS, "job_id": "24363", "source_revision": pilot.REVISION,
        "postflight_pins_verified": True, "defaults_promoted": False, "goal_completed": False,
        "execution_plan": record(root / "execution_plan.json"), "execution": execution,
        "coverage": {k: native[k] for k in ("genes", "species", "groups", "pairs")},
        "native_products": native["products"], "native_summary": native["summary"], "route_evidence": routes,
        "control": deepcopy(control), "scores": deepcopy(raw), "details": deepcopy(details),
        "scoring": {"scores": deepcopy(raw), "details": deepcopy(details["family_details"])},
        "gate": pilot.pilot_gate(raw, details, control), "total_pilot_wall_s": 40,
        "host_before": {"competing_cores": 40}, "host_after": {"competing_cores": 50}}
    write_json(root / "results.json", report)
    monkeypatch.setattr(pilot, "validate_variant", lambda *args: (deepcopy(native), routes, [pin]))
    def reader(repo, dataset, scored, actual, path, participant):
        assert dataset == "OrthoBench" and actual == native and path == folder
        assert participant == "unused_for_orthobench"
        return deepcopy(raw), deepcopy(details), [pin]
    monkeypatch.setattr(pilot, "readback_scores", reader)
    monkeypatch.setattr(pilot.subprocess, "check_output", lambda *args, **kwargs: ACCOUNTING)
    return root, plan, report, native


def test_completed_cell_checks_native_routes_and_independent_pr_gate(tmp_path, completed):
    root, plan, report, _ = completed
    cell = pilot.cell(tmp_path, root, plan, report)
    assert cell["scores"] == {ENDPOINTS[0]: .75}
    assert cell["gate"]["passes_development_gate"] is True
    assert cell["source_revision"] == pilot.REVISION
    assert cell["resources"]["native"] == report["execution"]
    assert cell["resources"]["host_before"] == {"competing_cores": 40}


def test_native_validation_failure_never_falls_back_to_reported_scores(tmp_path, completed, monkeypatch):
    root, plan, report, _ = completed
    def reject(*args):
        raise ValueError("Native validation failed")
    monkeypatch.setattr(pilot, "validate_variant", reject)
    monkeypatch.setattr(pilot, "readback_scores", lambda *args: pytest.fail("Invalid native scored"))
    with pytest.raises(ValueError, match="Native validation failed"):
        pilot.cell(tmp_path, root, plan, report)


@pytest.mark.parametrize("axis", ["precision", "recall"])
def test_independent_gate_rejects_pr_decline_despite_f1_gain(tmp_path, completed, monkeypatch, axis):
    root, plan, report, _ = completed
    report["details"][axis] = plan["control"]["details"][axis] - .006
    report["scoring"]["details"] = report["details"]["family_details"]
    monkeypatch.setattr(pilot, "readback_scores", lambda *args: (report["scores"], report["details"], []))
    # The producer's unchanged positive gate is not trusted.
    with pytest.raises(ValueError, match="precommitted gate"):
        pilot.cell(tmp_path, root, plan, report)
    report["gate"] = pilot.pilot_gate(report["scores"], report["details"], plan["control"])
    admitted = pilot.cell(tmp_path, root, plan, report)
    assert admitted["gate"]["passes_development_gate"] is False


@pytest.mark.parametrize("field,value", [("status", "running"), ("postflight_pins_verified", False),
    ("defaults_promoted", True), ("goal_completed", True), ("source_revision", "other"),
    ("error", "failed"), ("postflight_error", "changed")])
def test_unverified_success_claim_refused(tmp_path, completed, field, value):
    root, plan, report, _ = completed
    report[field] = value
    with pytest.raises(ValueError, match="successfully verified"):
        pilot.cell(tmp_path, root, plan, report)


@pytest.mark.parametrize("field", ["native_products", "native_summary", "coverage", "route_evidence", "control"])
def test_native_or_control_mismatch_refused(tmp_path, completed, field):
    root, plan, report, _ = completed
    report[field] = {}
    with pytest.raises(ValueError, match="complete native outputs"):
        pilot.cell(tmp_path, root, plan, report)


@pytest.mark.parametrize("field", ["scores", "details", "gate", "scoring_scores", "scoring_details"])
def test_score_details_or_gate_cannot_be_changed(tmp_path, completed, field):
    root, plan, report, _ = completed
    if field.startswith("scoring_"):
        report["scoring"][field.split("_", 1)[1]] = {}
    else:
        report[field] = {}
    with pytest.raises(ValueError, match="precommitted gate"):
        pilot.cell(tmp_path, root, plan, report)


@pytest.mark.parametrize("field", ["returncode", "command", "time", "log", "resources"])
def test_actual_command_and_resource_pins_cannot_be_changed(tmp_path, completed, field):
    root, plan, report, _ = completed
    report["execution"][field] = "changed"
    with pytest.raises(ValueError, match="native command"):
        pilot.cell(tmp_path, root, plan, report)


@pytest.mark.parametrize("wall", [None, True, 0, -1, float("nan"), float("inf")])
def test_invalid_wall_time_refused(tmp_path, completed, wall):
    root, plan, report, _ = completed
    report["total_pilot_wall_s"] = wall
    with pytest.raises(ValueError, match="wall time"):
        pilot.cell(tmp_path, root, plan, report)


@pytest.mark.parametrize("field", ["genes", "species"])
def test_incomplete_full_proteome_universe_refused(tmp_path, completed, field):
    root, plan, report, native = completed
    native[field] = 1
    with pytest.raises(ValueError, match="complete native outputs"):
        pilot.cell(tmp_path, root, plan, report)


@pytest.fixture
def planned(tmp_path, monkeypatch):
    root, source = tmp_path / "producer", tmp_path / "executor"
    root.mkdir()
    (source / "benchmark_tools").mkdir(parents=True)
    for name in (pilot.PROTOCOL, pilot.SLURM):
        (source / "benchmark_tools" / name).write_text("fixture\n")
    artifact = tmp_path / "artifact"
    artifact.write_text("pinned\n")
    pin = record(artifact)
    plan = {"repo": str(tmp_path), "output": str(root), "source_root": str(source),
        "source_revision": pilot.REVISION, "defaults_promoted": False, "goal_completed": False,
        "inputs": [pin], "sources": [pin, *[record(source / "benchmark_tools" / name)
            for name in (pilot.PROTOCOL, pilot.SLURM)]], "runtime": {"files": [pin]},
        "control": {"scores": "verified"}, "dependency_accounting": ["terminal"], "case": {"full": True}}
    write_json(root / "execution_plan.json", plan)
    launch = {"job_id": "24363", "source_root": str(source), "source_revision": pilot.REVISION,
        "execution_plan": record(root / "execution_plan.json"), "working_tests": pin, "frozen_tests": pin}
    report = {"job_id": "24363", "execution_plan": launch["execution_plan"], "source_revision": pilot.REVISION}
    def receipt(repo, executor, revision, name):
        assert repo == executor == tmp_path and revision == pilot.LAUNCH_REVISION and name == pilot.SUBMISSION
        return launch, pin
    monkeypatch.setattr(pilot, "committed", receipt)
    monkeypatch.setattr(pilot.subprocess, "check_output", lambda *args, **kwargs: pilot.REVISION + "\n")
    monkeypatch.setattr(pilot, "sources", lambda *args: [pin])
    monkeypatch.setattr(pilot, "runtime", lambda: {"files": [pin]})
    monkeypatch.setattr(pilot, "prerequisites", lambda *args: ({}, {}, deepcopy(plan["control"]), [pin], ["terminal"]))
    monkeypatch.setattr(pilot, "scientific_source", lambda *args: [pin])
    monkeypatch.setattr(pilot, "validate_plan", lambda *args: None)
    monkeypatch.setattr(pilot, "verify_input_inventory", lambda *args: None)
    return root, report, launch, plan, artifact


def test_load_plan_authenticates_committed_launch_and_rechecks_science_control(tmp_path, planned):
    root, report, _, expected, _ = planned
    plan, pins, failed = pilot.load_plan(tmp_path, root, report)
    assert plan == expected and pins and failed == []


@pytest.mark.parametrize("field,value", [("job_id", "other"), ("source_revision", "other"), ("execution_plan", {})])
def test_report_identity_mismatch_refused(tmp_path, planned, field, value):
    root, report, _, _, _ = planned
    report[field] = value
    with pytest.raises(ValueError, match="committed launch"):
        pilot.load_plan(tmp_path, root, report)


@pytest.mark.parametrize("field,value", [("repo", "/other"), ("output", "/other"),
    ("source_root", "/other"), ("source_revision", "other"),
    ("defaults_promoted", True), ("goal_completed", True)])
def test_plan_identity_mismatch_refused_even_with_new_plan_pin(tmp_path, planned, field, value):
    root, report, launch, plan, _ = planned
    plan[field] = value
    write_json(root / "execution_plan.json", plan)
    report["execution_plan"] = launch["execution_plan"] = record(root / "execution_plan.json")
    with pytest.raises(ValueError):
        pilot.load_plan(tmp_path, root, report)


@pytest.mark.parametrize("change", ["pin", "runtime", "inventory", "control", "accounting", "science", "plan", "input"])
def test_full_readback_checks_each_provenance_boundary(tmp_path, planned, monkeypatch, change):
    root, report, _, _, artifact = planned
    if change == "pin":
        artifact.write_text("changed\n")
    elif change == "runtime":
        monkeypatch.setattr(pilot, "runtime", lambda: {})
    elif change == "inventory":
        monkeypatch.setattr(pilot, "sources", lambda *args: [])
    elif change in ("control", "accounting"):
        monkeypatch.setattr(pilot, "prerequisites", lambda *args: ({}, {}, {}, [], []))
    else:
        def reject(*args):
            raise ValueError("Actual " + change + " changed")
        monkeypatch.setattr(pilot, {"science": "scientific_source", "plan": "validate_plan",
                                  "input": "verify_input_inventory"}[change], reject)
    with pytest.raises(ValueError):
        pilot.load_plan(tmp_path, root, report)
