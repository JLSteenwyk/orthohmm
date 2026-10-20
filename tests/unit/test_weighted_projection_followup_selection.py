from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import plan_weighted_projection_followup as planner
from benchmark_tools import report_profile_evidence_pilot as base
from benchmark_tools import report_weighted_projection as reader
from benchmark_tools import weighted_projection_experiment as experiment
from benchmark_tools.run_profile_weighting_pilot import advancement
from tests.unit.test_weighted_projection_readback import cells_and_rows


def admitted(passing=(True, True)):
    cells, rows = cells_and_rows()
    controls = {c["label"]: experiment.point(c) for c in cells[:4]}
    for index, passes in enumerate(passing, start=4):
        if not passes:
            cells[index]["scores"][base.ENDPOINTS[0]] = controls["repair_none"]["f1"]
        rows[index]["advancement"] = experiment.matched_gate(
            cells[index]["label"], experiment.point(cells[index]), controls, advancement)
    gates = reader.verify_panel(cells, rows)
    for cell in cells:
        cell.update(source_revision=planner.SCIENCE, coverage={"genes": 251378, "species": 12})
    report = {"status": "weighted_projection_independently_read_back", "accuracy_result_admitted": True,
        "source_revision": planner.SCIENCE, "reporter_revision": planner.SCIENCE,
        "cells": cells, "matched_gates": gates,
        "advancing_labels": [key for key, gate in gates.items() if gate["eligible"]],
        "defaults_promoted": False, "goal_completed": False}
    return report, {"arms": rows}


@pytest.mark.parametrize("passing", [(False, False), (True, False), (False, True), (True, True)])
def test_exact_matched_closure_for_zero_one_or_two_passing_combinations(passing):
    report, producer = admitted(passing)
    before = deepcopy(report)
    chosen = planner.selection(report, producer)
    count = sum(passing)
    assert len(chosen["execution_arms"]) == (0 if count == 0 else 4 if count == 1 else 6)
    candidates = [label for label, passes in zip(experiment.LABELS[4:], passing) if passes]
    assert chosen["advancing_candidates"] == candidates
    assert chosen["stopped_candidates"] == [label for label in experiment.LABELS[4:] if label not in candidates]
    needed = set(candidates)
    for label in candidates: needed.update(experiment.configuration(label)["matched_controls"])
    assert chosen["execution_arms"] == [label for label in experiment.LABELS if label in needed]
    assert chosen["arm_parameters"] == [experiment.configuration(label) for label in chosen["execution_arms"]]
    for label in chosen["execution_arms"]:
        assert chosen["roles"][label] == experiment.configuration(label)["role"]
        assert chosen["control_dependencies"][label] == (
            [l for l in chosen["execution_arms"] if l not in candidates] if label in candidates else [])
    assert report == before


@pytest.mark.parametrize("kind", ["running", "failed", "admission", "source", "reader", "promoted", "complete",
    "cell_source", "universe", "missing", "order", "hidden_winner", "revived_failure", "gate", "role", "inactive"])
def test_partial_stale_or_relabelled_panels_cannot_select_transfers(kind):
    report, producer = admitted((True, False))
    if kind in ("running", "failed"): report["status"] = kind
    elif kind == "admission": report["accuracy_result_admitted"] = False
    elif kind == "source": report["source_revision"] = "old"
    elif kind == "reader": report["reporter_revision"] = "new"
    elif kind == "promoted": report["defaults_promoted"] = True
    elif kind == "complete": report["goal_completed"] = True
    elif kind == "cell_source": report["cells"][4]["source_revision"] = "old"
    elif kind == "universe": report["cells"][4]["coverage"]["genes"] -= 1
    elif kind == "missing": report["cells"].pop()
    elif kind == "order": producer["arms"].reverse()
    elif kind == "hidden_winner": report["advancing_labels"] = []
    elif kind == "revived_failure": report["advancing_labels"].append(experiment.LABELS[5])
    elif kind == "gate": report["matched_gates"][experiment.LABELS[4]]["eligible"] = False
    elif kind == "role": report["cells"][0]["activation"]["role"] = "candidate"
    else: report["cells"][4]["activation"]["joint_profiles"] = 0
    with pytest.raises(ValueError): planner.selection(report, producer)


def test_no_applicability_fallback_or_failed_branch_revival_is_added():
    report, producer = admitted((False, False))
    chosen = planner.selection(report, producer)
    assert not chosen["execution_arms"] and not chosen["roles"] and not chosen["control_dependencies"]
    assert chosen["stopped_candidates"] == list(experiment.LABELS[4:])
    assert "applicability_candidates" not in chosen


def launch_fixture(tmp_path, monkeypatch, kind=None):
    binding = {"revision": "f" * 40, "entry_key": "weighted_projection_pilot"}
    science = {"runner_commit": planner.SCIENCE, "scientific_source": "science",
        "inference": {"job_id": "24538"}, "independent_reader": {
            "runner_commit": planner.SCIENCE, "scientific_reader_source": "reader"}}
    entry = {"runner_commit": planner.SCIENCE, "scientific_reader_source": "reader", "launch_binding": deepcopy(binding),
        "job": {"submitted": True, "job_id": "24539", "dependency": "afterany:24538",
                "cpus": 1, "memory_allocation": "16G", "wall_limit": "02:00:00"}}
    if kind == "binding": entry["launch_binding"]["entry_key"] = "different"
    elif kind == "source": science["runner_commit"] = "old"
    elif kind == "reader": entry["runner_commit"] = "old"
    elif kind == "reader_path": entry["scientific_reader_source"] = "other"
    elif kind == "same_source": science["scientific_source"] = "reader"
    elif kind == "job": science["inference"]["job_id"] = "old"
    elif kind == "dependency": entry["job"]["dependency"] = "afterany:old"
    elif kind == "unsubmitted": entry["job"]["submitted"] = False
    elif kind == "cpus": entry["job"]["cpus"] = 16
    elif kind == "memory": entry["job"]["memory_allocation"] = "128G"
    elif kind == "wall": entry["job"]["wall_limit"] = "1-00:00:00"
    blob = json.dumps({"weighted_projection_postrun_readback": entry}).encode()
    monkeypatch.setattr(planner.subprocess, "check_output", lambda *a, **kw: blob)
    monkeypatch.setattr(reader, "launch", lambda *a: (science, {"path": "/plan"}, binding))
    return science, entry


@pytest.mark.parametrize("kind", [None, "binding", "source", "reader", "reader_path", "same_source", "job",
                                 "dependency", "unsubmitted", "cpus", "memory", "wall"])
def test_authentication_uses_the_exact_original_submitted_science_and_separate_reader(tmp_path, monkeypatch, kind):
    science, entry = launch_fixture(tmp_path, monkeypatch, kind)
    if kind:
        with pytest.raises(ValueError): planner.launch(tmp_path)
    else:
        actual, observed, pin, provenance = planner.launch(tmp_path)
        assert actual == entry and observed == science and pin == {"path": "/plan"}
        assert provenance["reader_submission"]["revision"] == planner.SUBMISSION


def test_live_or_failed_reader_cannot_be_replaced_by_current_partial_producer(tmp_path, monkeypatch):
    monkeypatch.setattr(planner, "launch", lambda _: ({}, {}, {}, {}))
    calls = []
    def read(*args):
        calls.append(True)
        raise ValueError("Original reader is not terminal successful")
    monkeypatch.setattr(planner.campaign, "read_admission", read)
    with pytest.raises(ValueError, match="terminal successful"):
        planner.authenticate(tmp_path)
    assert calls == [True] and not list(tmp_path.iterdir())


def test_previous_preparation_and_unfrozen_main_checkout_are_rejected(tmp_path):
    source = Path(planner.__file__).resolve().parents[1]
    with pytest.raises(FileExistsError, match="Preserve every"):
        planner.prepare(source, source, tmp_path)
    with pytest.raises(ValueError, match="frozen committed"):
        planner.prepare(source, source, tmp_path / "new")


@pytest.mark.parametrize("kind", [None, "plan", "binding", "source", "inventory", "runtime"])
def test_full_admission_rechecks_producer_plan_frozen_source_and_runtime(tmp_path, monkeypatch, kind):
    report, producer = admitted()
    source = tmp_path / "science"
    (source / "benchmark_tools").mkdir(parents=True)
    slurm = source / "benchmark_tools/run_profile_weighting_pilot.slurm"
    slurm.write_text("fixture\n")
    from benchmark_tools.prepare_ob_candidate_neighborhood import record
    plan = {"source_root": str(source), "source_revision": planner.SCIENCE, "runtime": {"files": []},
        "sources": [record(slurm)], "arms": list(experiment.LABELS), "cpu": 16, "experiment": "weighted_projection"}
    plan_path = tmp_path / "plan.json"
    plan_path.write_text(json.dumps(plan))
    pin = record(plan_path)
    for row in producer["arms"]: row.update(status="independent_content_readback_verified", dataset="OrthoBench")
    producer.update(status="finished_all_arms_verified", execution_plan=pin, source_revision=planner.SCIENCE,
        experiment="weighted_projection", job_id="24538", defaults_promoted=False, goal_completed=False)
    producer_path = tmp_path / "producer.json"
    producer_path.write_text(json.dumps(producer))
    binding = {"science": {"revision": "f" * 40}}
    report.update(launch_ledger=binding["science"], execution_plan=pin, producer=record(producer_path))
    if kind == "plan": report["execution_plan"] = {}
    elif kind == "binding": report["launch_ledger"] = {}
    science = {"scientific_source": "science", "inference": {"job_id": "24538"}}
    monkeypatch.setattr(planner, "launch", lambda _: ({}, science, pin, binding))
    monkeypatch.setattr(planner.campaign, "read_admission", lambda *a: (report, [], {"reader": {"path": "/reader"}}))
    monkeypatch.setattr(planner, "successful_job", lambda *a: {"State": "COMPLETED"})
    monkeypatch.setattr(planner.subprocess, "check_output", lambda *a, **kw: "other" if kind == "source" else planner.SCIENCE)
    monkeypatch.setattr(planner, "sources", lambda _: [{"changed": True}] if kind == "inventory" else [])
    monkeypatch.setattr(planner, "runtime", lambda: {"changed": True} if kind == "runtime" else plan["runtime"])
    if kind:
        with pytest.raises(ValueError): planner.authenticate(tmp_path)
    else:
        result, original, chosen, pins, provenance = planner.authenticate(tmp_path)
        assert result == report and original == plan and chosen["advancing_candidates"] == list(experiment.LABELS[4:])
        assert pins and provenance["launches"] == binding


@pytest.mark.parametrize("passing", [(False, False), (True, False), (False, True), (True, True)])
def test_frozen_planner_writes_only_selection_not_scientific_case_or_accuracy_claim(tmp_path, monkeypatch, passing):
    from benchmark_tools.prepare_ob_candidate_neighborhood import record
    driver = Path(planner.__file__).resolve().parents[1]
    report, producer = admitted(passing)
    chosen = planner.selection(report, producer)
    plan_path = tmp_path / "original_plan.json"
    plan_path.write_text("{}\n")
    report["execution_plan"] = record(plan_path)
    pilot = {"source_root": "/frozen/science", "runtime": {"files": []}, "tool_path": "/tools"}
    provenance = {"reader": {"reader": {"path": "/reader/summary.json"}}}
    monkeypatch.setattr(planner, "authenticate", lambda _: (report, pilot, chosen, [], provenance))
    monkeypatch.setattr(planner, "sources", lambda _: [])
    def git(command, **kwargs):
        if command[1] == "rev-parse": return "b" * 40
        relative = command[2].split(":", 1)[1]
        return (driver / relative).read_bytes()
    monkeypatch.setattr(planner.subprocess, "check_output", git)
    output = tmp_path / "new"
    result = planner.prepare(tmp_path, driver, output)
    assert result["selection"] == chosen
    assert result["datasets"] == (list(planner.DATASETS) if any(passing) else [])
    assert result["status"] == ("weighted_projection_followup_selection_prepared_execution_driver_required"
        if any(passing) else "no_weighted_projection_followup_required")
    assert not result["goal_completed"] and not result["defaults_promoted"]
    assert "cases" not in result and "accuracy_result_admitted" not in result
    assert sorted(p.name for p in output.iterdir()) == ["selection.json"]
    assert json.loads((output / "selection.json").read_text()) == result
    with pytest.raises(FileExistsError): planner.prepare(tmp_path, driver, output)


@pytest.mark.parametrize("kind", ["prerequisite", "protocol", "source"])
def test_failed_authentication_or_changed_planner_creates_no_selection(tmp_path, monkeypatch, kind):
    driver = Path(planner.__file__).resolve().parents[1]
    report, producer = admitted()
    pilot = {"source_root": "/science", "runtime": {}, "tool_path": "/tools"}
    if kind == "prerequisite":
        def fail(_): raise ValueError("Original reader still running or failed")
        monkeypatch.setattr(planner, "authenticate", fail)
    else:
        monkeypatch.setattr(planner, "authenticate", lambda _: (report, pilot, planner.selection(report, producer), [], {}))
    monkeypatch.setattr(planner, "sources", lambda _: [])
    def git(command, **kwargs):
        if command[1] == "rev-parse": return "b" * 40
        revision, relative = command[2].split(":", 1)
        if kind == "protocol" or kind == "source" and revision == "b" * 40: return b"changed"
        return (driver / relative).read_bytes()
    monkeypatch.setattr(planner.subprocess, "check_output", git)
    output = tmp_path / "new"
    with pytest.raises(ValueError): planner.prepare(tmp_path, driver, output)
    assert not output.exists()
