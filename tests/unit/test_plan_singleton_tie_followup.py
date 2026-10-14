from copy import deepcopy

import pytest

from benchmark_tools import plan_singleton_tie_followup as planner
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_membership_hillclimb import write_json


def admitted(variant=.7365896292967041):
    metric = "OrthoBench: OrthoBench F1"
    control = {"singleton_tie_policy": "first_hit", "source_revision": planner.PILOT_REVISION,
               "status": "independent_content_readback_verified", "coverage": {"genes": 251378, "species": 12},
               "scores": {metric: .7410607351873405}}
    cell = {**deepcopy(control), "singleton_tie_policy": "defer_ambiguous", "dataset": "OrthoBench",
            "scores": {metric: variant}, "matched_control": control}
    method = {"key": "singleton_tie:defer_ambiguous", "scores": cell["scores"], "cells": [cell]}
    return {"methods": [method]}


def test_actual_pilot_regression_narrowly_passes_development_gate():
    gate = planner.pilot_gate(admitted())
    assert gate["eligible"] is True
    assert gate["delta"] == pytest.approx(-.004471105890636373)
    assert "not equivalence or improvement" in gate["scope"]


def test_larger_decline_rejects_followup():
    assert planner.pilot_gate(admitted(.735))["eligible"] is False


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -.1, 1.1])
def test_invalid_f1_never_passes_gate(value):
    summary = admitted(value)
    with pytest.raises(ValueError, match="finite|scores disagree"):
        planner.pilot_gate(summary)


@pytest.mark.parametrize("change", [{"singleton_tie_policy": "first_hit"}, {"source_revision": "other"},
                                    {"coverage": {"genes": 251377, "species": 12}}])
def test_changed_variant_identity_rejected(change):
    summary = admitted()
    summary["methods"][0]["cells"][0].update(change)
    with pytest.raises(ValueError, match="Pilot policy"):
        planner.pilot_gate(summary)


@pytest.mark.parametrize("change", [{"singleton_tie_policy": "defer_ambiguous"}, {"source_revision": "other"},
                                    {"status": "failed"}, {"coverage": {"genes": 251378, "species": 11}}])
def test_changed_or_failed_matched_control_rejected(change):
    summary = admitted()
    summary["methods"][0]["cells"][0]["matched_control"].update(change)
    with pytest.raises(ValueError, match="Pilot policy"):
        planner.pilot_gate(summary)


@pytest.mark.parametrize("status", ["failed", "running", "not_run_failed_control_gate"])
def test_unverified_cells_cannot_trigger_followup(status):
    summary = admitted()
    summary["methods"][0]["cells"][0]["status"] = status
    with pytest.raises(ValueError, match="complete dataset cell"):
        planner.pilot_gate(summary)


def test_duplicate_admitted_cells_rejected():
    summary = admitted()
    summary["methods"][0]["cells"] *= 2
    with pytest.raises(ValueError, match="complete dataset cell"):
        planner.pilot_gate(summary)


def test_method_and_admitted_cell_scores_must_agree():
    summary = admitted()
    summary["methods"][0]["scores"] = {"OrthoBench: OrthoBench F1": 1.0}
    with pytest.raises(ValueError, match="scores disagree"):
        planner.pilot_gate(summary)


@pytest.mark.parametrize("accounting", ["24250|RUNNING|0:0\n", "24250|FAILED|1:0\n", ""])
def test_preparation_waits_for_successful_terminal_pilot(tmp_path, monkeypatch, accounting):
    summary = admitted()
    table = tmp_path / "scores.tsv"
    table.write_text("fixture\n")
    summary["score_table"] = record(table)
    readback = tmp_path / "summary.json"
    write_json(readback, summary)
    monkeypatch.setattr(planner.subprocess, "check_output", lambda *args, **kwargs: accounting)
    with pytest.raises(ValueError, match="successful terminal"):
        planner.prepare(tmp_path, readback, tmp_path / "output")
    assert not (tmp_path / "output").exists()


def test_rejected_pilot_does_not_start_processes(tmp_path, monkeypatch):
    summary = admitted(.735)
    table = tmp_path / "scores.tsv"
    table.write_text("fixture\n")
    summary["score_table"] = record(table)
    readback = tmp_path / "summary.json"
    write_json(readback, summary)
    monkeypatch.setattr(planner.subprocess, "check_output", lambda *args, **kwargs: pytest.fail("Process started"))
    with pytest.raises(ValueError, match="decline screen"):
        planner.prepare(tmp_path, readback, tmp_path / "output")
