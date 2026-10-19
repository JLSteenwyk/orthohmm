from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import plan_projection_followup as planner
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_projection_recovery_readback import cells_and_rows


def evidence(passing=None, inactive=()):
    cells, rows = cells_and_rows()
    if passing is not None:
        for i, label in ((1, planner.CANDIDATES[0]), (2, planner.CANDIDATES[1]), (4, planner.CANDIDATES[2])):
            if label not in passing:
                cells[i]["scores"] = deepcopy(cells[0]["scores"])
                cells[i]["details"] = deepcopy(cells[0]["details"])
        cells[3]["scores"] = deepcopy(cells[1]["scores"])
        cells[3]["details"] = deepcopy(cells[1]["details"])
    for i, cell in enumerate(cells):
        cell["source_revision"] = planner.SCIENCE
        if cell["label"] in inactive:
            cell["activation"].update(completely_measured_queries=0, applicability="no_complete_measurements")
    control = planner.recovery_reader.point(cells[0])
    projection = planner.recovery_reader.point(cells[1])
    for i in (1, 2, 4):
        rows[i]["advancement"] = (planner.campaign.advancement(control, planner.recovery_reader.point(cells[i]))
                                    if i == 1 else planner.pilot_reader.selective_gate(control, cells[i]))
        if i != 1:
            rows[i].update(matched_projection_advancement=planner.campaign.advancement(projection, planner.recovery_reader.point(cells[i])),
                           accuracy_interpretation=cells[i]["activation"]["applicability"])
    contrast = planner.recovery_reader.verify_panel(cells, rows)
    report = {"status": "projection_recovery_independently_read_back", "accuracy_result_admitted": True,
        "source_revision": planner.SCIENCE, "reporter_revision": planner.READER, "cells": cells,
        "backend_contrast": contrast, "advancing_labels": [cells[i]["label"] for i in (1, 2, 4) if rows[i]["advancement"]["eligible"]],
        "defaults_promoted": False, "goal_completed": False}
    return report, {"arms": rows}


def test_every_passing_candidate_is_retained_with_control_roles():
    report, producer = evidence()
    chosen = planner.selection(report, producer)
    assert chosen["advancing_candidates"] == list(planner.CANDIDATES)
    assert chosen["execution_arms"] == list(planner.recovery_reader.ARMS)
    assert chosen["roles"]["projection_hmmer_score_only"] == "measurement_control"
    assert "projection_hmmer_score_only" not in chosen["advancing_candidates"]
    assert chosen["control_dependencies"]["projection_hmmer_reselect"] == [
        "legacy_off", "projection_off", "projection_builtin_reselect", "projection_hmmer_score_only"]


@pytest.mark.parametrize("passing,arms", [
    (["projection_off"], ["legacy_off", "projection_off"]),
    (["projection_builtin_reselect"], ["legacy_off", "projection_off", "projection_builtin_reselect"]),
    (["projection_hmmer_reselect"], list(planner.recovery_reader.ARMS)),
    ([], []),
])
def test_minimum_matched_control_closure_does_not_promote_failed_candidates(passing, arms):
    report, producer = evidence(passing)
    chosen = planner.selection(report, producer)
    assert chosen["advancing_candidates"] == passing and chosen["execution_arms"] == arms
    assert chosen["stopped_candidates"] == [label for label in planner.CANDIDATES if label not in passing]
    if passing == ["projection_hmmer_reselect"]:
        assert chosen["roles"]["projection_off"] == "matched_construction_control"
        assert chosen["roles"]["projection_builtin_reselect"] == "matched_decision_control"
        assert not chosen["original_gates"]["projection_off"]["eligible"]
        assert not chosen["original_gates"]["projection_builtin_reselect"]["eligible"]


def test_original_inactivity_fallback_is_distinct_from_a_pilot_pass():
    report, producer = evidence([], inactive=["projection_hmmer_reselect"])
    chosen = planner.selection(report, producer)
    assert chosen["advancing_candidates"] == []
    assert chosen["applicability_candidates"] == ["projection_hmmer_reselect"]
    assert chosen["roles"]["projection_hmmer_reselect"] == "applicability_followup"
    assert chosen["execution_arms"] == list(planner.recovery_reader.ARMS)


@pytest.mark.parametrize("kind", ["partial", "unadmitted", "source", "reader", "promoted", "closed", "universe", "cell_source",
    "missing", "order", "digest", "measurement", "gate", "matched_gate", "backend", "chosen_best_only"])
def test_partial_changed_or_cherry_picked_evidence_cannot_select_transfer(kind):
    report, producer = evidence()
    if kind == "partial": report["status"] = "running"
    elif kind == "unadmitted": report["accuracy_result_admitted"] = False
    elif kind == "source": report["source_revision"] = "other"
    elif kind == "reader": report["reporter_revision"] = "other"
    elif kind == "promoted": report["defaults_promoted"] = True
    elif kind == "closed": report["goal_completed"] = True
    elif kind == "universe": report["cells"][2]["coverage"]["genes"] -= 1
    elif kind == "cell_source": report["cells"][2]["source_revision"] = "other"
    elif kind == "missing": report["cells"].pop()
    elif kind == "order": report["cells"].reverse()
    elif kind == "digest": report["cells"][4]["activation"]["candidate_digest"] = "other"
    elif kind == "measurement": report["cells"][3]["rooted_species_tree"]["sha256"] = "other"
    elif kind == "gate": producer["arms"][4]["advancement"]["eligible"] = False
    elif kind == "matched_gate": producer["arms"][4]["matched_projection_advancement"]["eligible"] = True
    elif kind == "backend": report["backend_contrast"]["f1"] = 1.
    else: report["advancing_labels"] = ["projection_off"]
    with pytest.raises(ValueError): planner.selection(report, producer)


def authenticated_fixture(tmp_path, monkeypatch):
    report, producer = evidence()
    science = tmp_path / "science"
    (science / "benchmark_tools").mkdir(parents=True)
    slurm = science / "benchmark_tools/run_profile_weighting_pilot.slurm"
    slurm.write_text("science\n")
    pilot = {"source_root": str(science), "source_revision": planner.SCIENCE,
        "sources": [record(slurm)], "runtime": {"files": []}, "tool_path": "/tools", "cpu": 16,
        "experiment": "projection_recovery", "arms": list(planner.recovery_reader.ARMS)}
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(pilot))
    producer.update(job_id="24514", status="finished_all_arms_verified", source_revision=planner.SCIENCE,
        experiment="projection_recovery", defaults_promoted=False, goal_completed=False, execution_plan=record(path))
    for row in producer["arms"]: row.update(dataset="OrthoBench", status="independent_content_readback_verified")
    product = tmp_path / "producer.json"
    product.write_text(json.dumps(producer))
    report.update(execution_plan=record(path), producer=record(product), launch_ledger={"revision": "bound"})
    report_path = tmp_path / "report.json"
    report_path.write_text(json.dumps(report))
    entry = {"runner_commit": planner.READER, "launch_binding": report["launch_ledger"]}
    science_entry = {"runner_commit": planner.SCIENCE, "scientific_source": "science", "inference": {"job_id": "24514"}}
    monkeypatch.setattr(planner.recovery_campaign, "launch", lambda _: ({"recovery": entry}, science_entry, {"fixture": "immutable"}))
    monkeypatch.setattr(planner.campaign, "read_admission", lambda *a: (report, [record(report_path)], {"reader": record(report_path)}))
    monkeypatch.setattr(planner, "successful_job", lambda *a: {"fixture": list(a)})
    monkeypatch.setattr(planner, "sources", lambda _: [])
    monkeypatch.setattr(planner, "runtime", lambda: {"files": []})
    monkeypatch.setattr(planner.subprocess, "check_output", lambda *a, **k: planner.SCIENCE)
    return report, producer, pilot, product


def test_authentication_binds_complete_reader_parent_plan_and_frozen_runtime(tmp_path, monkeypatch):
    authenticated_fixture(tmp_path, monkeypatch)
    report, pilot, chosen, pins, provenance = planner.authenticate(tmp_path)
    assert chosen["advancing_candidates"] == list(planner.CANDIDATES) and pins
    assert provenance["science_accounting"]["fixture"] == ["24514", 16, "128G"]
    assert pilot["source_revision"] == report["source_revision"] == planner.SCIENCE
    for pin in pins: planner.check(pin)


@pytest.mark.parametrize("kind", ["pending_reader", "binding", "producer", "runtime", "source_inventory"])
def test_no_active_job_or_changed_input_can_pass_authentication(tmp_path, monkeypatch, kind):
    report, producer, pilot, product = authenticated_fixture(tmp_path, monkeypatch)
    if kind == "pending_reader":
        monkeypatch.setattr(planner.campaign, "read_admission", lambda *a: (_ for _ in ()).throw(ValueError("Reader pending")))
    elif kind == "binding": report["launch_ledger"] = {"revision": "other"}
    elif kind == "producer":
        producer["status"] = "running"
        product.write_text(json.dumps(producer))
        report["producer"] = record(product)
    elif kind == "runtime": monkeypatch.setattr(planner, "runtime", lambda: {"changed": True})
    else: monkeypatch.setattr(planner, "sources", lambda _: [{"changed": True}])
    with pytest.raises(ValueError): planner.authenticate(tmp_path)


@pytest.mark.parametrize("eligible", [True, False])
@pytest.mark.parametrize("change", [None, "source", "protocol"])
def test_preparation_is_inputs_only_and_no_pass_means_no_dataset_work(tmp_path, monkeypatch, eligible, change):
    repo, driver, output = tmp_path / "repo", tmp_path / "driver", tmp_path / "prepared"
    repo.mkdir()
    folder = driver / "benchmark_tools"
    folder.mkdir(parents=True)
    original_file = Path(planner.__file__)
    original_protocol = json.loads(original_file.with_name(planner.PROTOCOL).read_text())
    if change == "protocol": original_protocol["resources"]["cpus_per_run"] = 32
    helper = folder / "plan_projection_followup.py"
    helper.write_text("frozen planner\n")
    (folder / planner.PROTOCOL).write_text(json.dumps(original_protocol))
    report, producer = evidence(None if eligible else [])
    chosen = planner.selection(report, producer)
    pin = record(helper)
    pilot = {"source_root": "/frozen/science", "runtime": {"files": []}, "tool_path": "/tools"}
    report["execution_plan"] = pin
    monkeypatch.setattr(planner, "__file__", str(helper))
    monkeypatch.setattr(planner, "sources", lambda _: [record(helper)])
    monkeypatch.setattr(planner, "authenticate", lambda _: (report, pilot, chosen, [], {"reader": {"reader": pin}}))
    calls = []

    def selected(api, repo, dataset):
        calls.append(dataset)
        checkpoint = tmp_path / dataset / "checkpoint"
        checkpoint.mkdir(parents=True)
        (checkpoint / "manifest.json").write_text("fixed\n")
        return {"dataset": dataset, "checkpoint": str(checkpoint), "case": {"baseline_output": str(tmp_path / dataset / "native")}}, []

    monkeypatch.setattr(planner, "selected_case", selected)

    def git(command, **kwargs):
        if command[1] == "rev-parse": return "driver_revision"
        return b"wrong source\n" if change == "source" else (driver / command[-1].split(":", 1)[1]).read_bytes()

    monkeypatch.setattr(planner.subprocess, "check_output", git)
    if change:
        with pytest.raises(ValueError): planner.prepare(repo, driver, output)
        assert calls == [] and not output.exists()
        return
    result = planner.prepare(repo, driver, output)
    assert calls == (list(planner.DATASETS) if eligible else [])
    assert result["status"] == ("projection_transfer_inputs_prepared_execution_driver_required" if eligible else "no_projection_followup_required")
    assert not result["defaults_promoted"] and not result["goal_completed"]
    assert json.loads((output / "prepared_inputs.json").read_text()) == result
    with pytest.raises(FileExistsError): planner.prepare(repo, driver, output)


def test_primary_or_existing_preparation_is_not_a_frozen_launch(tmp_path):
    with pytest.raises(FileExistsError): planner.prepare(tmp_path, tmp_path, tmp_path)
    with pytest.raises(ValueError, match="frozen"):
        planner.prepare(tmp_path, tmp_path, tmp_path / "out")
