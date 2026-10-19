from copy import deepcopy
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

from benchmark_tools import run_projection_followup as followup
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_membership_hillclimb import write_json
from tests.unit.test_projection_followup import authentication_fixture, execution_fixture


CAMPAIGN = "stockholm_annotation"


def stockholm_fixture(tmp_path, monkeypatch, **kwargs):
    api, plan, path, executed, scored = execution_fixture(tmp_path, monkeypatch, **kwargs)
    plan.update(campaign=CAMPAIGN, scientific_revision=followup.configuration(CAMPAIGN)["science"])
    write_json(path, plan)
    followup.validate_plan(plan)
    return api, plan, path, executed, scored


@pytest.mark.parametrize("dataset", ["QfO", "ThreeKingdoms"])
@pytest.mark.parametrize("passing", [None, ["projection_off"], ["projection_builtin_reselect"], ["projection_hmmer_reselect"]])
def test_rebound_execution_keeps_original_control_order_scores_and_scientific_flags(tmp_path, monkeypatch, dataset, passing):
    api, plan, path, executed, scored = stockholm_fixture(tmp_path, monkeypatch, dataset=dataset, passing=passing)
    result = followup.execute_arms(api, plan, path, ["a", "b", "c"], [0, 0, 0], ([], [], []), {})
    assert result["status"] == "all_projection_followup_arms_scored", result["arms"]
    assert result["source_revision"] == followup.configuration(CAMPAIGN)["science"]
    assert result["campaign"] == CAMPAIGN
    assert executed == [(label, stage) for label in plan["arms"] for stage in ("graph", "native")]
    assert len(scored) == len(plan["arms"]) - (2 if "projection_hmmer_score_only" in plan["arms"] else 1)
    for row in result["arms"]:
        assert row["role"] == plan["selection"]["roles"][row["label"]]
        graph = json.loads((Path(plan["output"]) / row["label"] / "preflight.json").read_text())["graph"]
        assert "--profile-sequence-weighting" not in graph and "--profile-alignment-guard" not in graph
        policy, mode = api.recovery_configuration(row["label"])
        if policy != "legacy": assert graph[graph.index("--profile-projection-policy") + 1] == policy
        if mode != "off": assert graph[graph.index("--selective-profile-mode") + 1] == mode
    if "projection_hmmer_score_only" in plan["arms"]:
        projection, measurement = result["arms"][1], result["arms"][3]
        assert all(projection[key] == measurement[key] for key in ("scores", "details", "metric_vector"))
        assert measurement["measurement_scoring_reuse"] == record(Path(plan["output"]) / "projection_off/assessment_receipt.json")


@pytest.mark.parametrize("failed,blocked", [("legacy_off", "projection_off"),
    ("projection_off", "projection_hmmer_score_only"), ("projection_builtin_reselect", "projection_hmmer_reselect"),
    ("projection_hmmer_score_only", "projection_hmmer_reselect")])
def test_rebound_source_never_bypasses_a_failed_control(tmp_path, monkeypatch, failed, blocked):
    api, plan, path, executed, _ = stockholm_fixture(tmp_path, monkeypatch, failed=(failed, "graph"))
    result = followup.execute_arms(api, plan, path, ["a", "b", "c"], [0, 0, 0], ([], [], []), {})
    assert result["status"] == "finished_with_explicit_gaps"
    assert next(row for row in result["arms"] if row["label"] == blocked)["status"] == "not_run_failed_required_control"
    assert not any(label == blocked for label, _ in executed)


@pytest.mark.parametrize("change", ["source", "missing_campaign", "unknown_campaign", "budget", "control"])
def test_new_plan_cannot_be_reinterpreted_as_original_or_change_gates(tmp_path, monkeypatch, change):
    _, plan, _, _, _ = stockholm_fixture(tmp_path, monkeypatch)
    if change == "source": plan["scientific_revision"] = followup.SCIENCE
    elif change == "missing_campaign": del plan["campaign"]
    elif change == "unknown_campaign": plan["campaign"] = "stockholm"
    elif change == "budget": plan["cpu"] = 32
    else: plan["selection"]["control_dependencies"]["projection_hmmer_reselect"] = []
    with pytest.raises(ValueError): followup.validate_plan(plan)


def new_authentication_fixture(tmp_path, monkeypatch):
    prepared, audit, calls = authentication_fixture(tmp_path, monkeypatch)
    binding = followup.configuration(CAMPAIGN)
    prepared.update(campaign=CAMPAIGN, scientific_revision=binding["science"], preparer_revision=binding["preparer"],
                    output=str((tmp_path / binding["prepared"]).parent))

    def process(command, **kwargs):
        if command[0] == "git": return binding["preparer"]
        assert kwargs["cwd"] == tmp_path / binding["planner"]
        assert kwargs["env"]["PYTHONPATH"] == str(tmp_path / binding["planner"])
        assert command[-1] == CAMPAIGN and "campaign_name=campaign" in command[2]
        assert "format_recovery" in command[2] and "species_tree_inference" in command[2]
        compile(command[2], "rebound_planner_authentication", "exec")
        return json.dumps(audit)

    monkeypatch.setattr(followup.subprocess, "check_output", process)
    return prepared, audit, calls


def test_reauthentication_uses_exact_new_frozen_planner_job_and_full_inventory(tmp_path, monkeypatch):
    prepared, audit, calls = new_authentication_fixture(tmp_path, monkeypatch)
    result = followup.authenticate_prepared(tmp_path, prepared)
    expected = {key: value for key, value in audit.items() if key != "prepared_evidence_digest"}
    assert result == {**expected, "campaign": CAMPAIGN}
    assert calls == [("24521", 1, "16G")]


@pytest.mark.parametrize("change", ["old_source", "old_preparer", "old_output", "missing_campaign", "unknown_campaign",
    "pending", "partial", "inventory", "case", "protocol", "selection"])
def test_reauthentication_rejects_stale_preparation_and_pending_or_changed_evidence(tmp_path, monkeypatch, change):
    prepared, _, _ = new_authentication_fixture(tmp_path, monkeypatch)
    if change == "old_source": prepared["scientific_revision"] = followup.SCIENCE
    elif change == "old_preparer": prepared["preparer_revision"] = followup.PREPARER
    elif change == "old_output": prepared["output"] = str((tmp_path / followup.PREPARED).parent)
    elif change == "missing_campaign": del prepared["campaign"]
    elif change == "unknown_campaign": prepared["campaign"] = "stockholm"
    elif change == "pending":
        monkeypatch.setattr(followup, "helper", lambda _: SimpleNamespace(successful_job=lambda *a: (_ for _ in ()).throw(ValueError("Preparation pending"))))
    elif change == "partial": prepared["status"] = "running"
    elif change == "inventory": prepared["inputs"] = []
    elif change == "case": prepared["cases"].pop()
    elif change == "protocol": prepared["protocol"]["sha256"] = "changed"
    else: prepared["selection"]["execution_arms"].pop()
    with pytest.raises(ValueError): followup.authenticate_prepared(tmp_path, prepared)


@pytest.mark.parametrize("change", [None, "uncommitted_driver", "old_receipt"])
def test_prepare_writes_only_explicit_new_bound_plan(tmp_path, monkeypatch, change):
    api, original, _, _, _ = stockholm_fixture(tmp_path, monkeypatch)
    repo, driver = Path(original["repo"]), tmp_path / "driver"
    binding = followup.configuration(CAMPAIGN)
    (driver / "benchmark_tools").mkdir(parents=True)
    for name in followup.driver_paths(CAMPAIGN): (driver / name).write_text("frozen\n")
    cases = [deepcopy(original["case"]) for _ in followup.UNIVERSES]
    for case, dataset in zip(cases, followup.UNIVERSES): case["dataset"] = dataset
    prepared = {"campaign": CAMPAIGN, "cases": cases, "selection": original["selection"], "inputs": original["inputs"],
                "protocol": record(driver / binding["protocol"])}
    if change == "old_receipt": prepared.pop("campaign")
    prepared_path = repo / binding["prepared"]
    prepared_path.parent.mkdir(parents=True)
    write_json(prepared_path, prepared)
    audit = {key: original[key] for key in ("campaign", "selection", "pilot_readback", "pilot_plan", "scientific_source", "runtime", "tool_path", "scientific_sources", "hmmer_tools")}
    monkeypatch.setattr(followup, "authenticate_prepared", lambda *a: audit)
    monkeypatch.setattr(followup, "science_api", lambda *a: api)
    monkeypatch.setattr(followup, "helper", lambda _: SimpleNamespace(selected_case=lambda api, repo, dataset: (next(c for c in cases if c["dataset"] == dataset), [])))
    monkeypatch.setattr(followup, "__file__", str(driver / "benchmark_tools/run_projection_followup.py"))
    api.require_profile_runtime = lambda *a: {"tiny_profile": "passed"}
    monkeypatch.setattr(followup.subprocess, "check_output", lambda cmd, **k: "new_driver" if cmd[1] == "rev-parse"
                        else b"changed\n" if change == "uncommitted_driver" else (driver / cmd[-1].split(":", 1)[1]).read_bytes())
    output = tmp_path / "new_execution"
    if change:
        with pytest.raises(ValueError): followup.prepare(repo, driver, output, "ThreeKingdoms", campaign_name=CAMPAIGN)
        assert not output.exists()
    else:
        plan = followup.prepare(repo, driver, output, "ThreeKingdoms", campaign_name=CAMPAIGN)
        assert plan["campaign"] == CAMPAIGN and plan["scientific_revision"] == binding["science"]
        assert plan["prepared_inputs"] == record(prepared_path)
        assert [pin["path"] for pin in plan["driver_sources"]] == [str(driver / name) for name in followup.driver_paths(CAMPAIGN)]
        assert json.loads((output / "execution_plan.json").read_text()) == plan
        assert not (output / "results.json").exists()


@pytest.mark.parametrize("change", [None, "old_receipt", "old_path", "old_source"])
def test_standalone_run_reauthenticates_new_campaign_before_any_science(tmp_path, monkeypatch, change):
    api, plan, path, executed, _ = stockholm_fixture(tmp_path, monkeypatch)
    repo, binding = Path(plan["repo"]), followup.configuration(CAMPAIGN)
    prepared_path = repo / (followup.PREPARED if change == "old_path" else binding["prepared"])
    prepared_path.parent.mkdir(parents=True)
    write_json(prepared_path, {"campaign": "projection" if change == "old_receipt" else CAMPAIGN})
    plan["prepared_inputs"] = record(prepared_path)
    plan["inputs"].append(plan["prepared_inputs"])
    write_json(path, plan)
    audit = {key: plan[key] for key in ("campaign", "selection", "pilot_readback", "pilot_plan", "scientific_source", "runtime", "tool_path", "scientific_sources", "hmmer_tools")}
    monkeypatch.setattr(followup, "authenticate_prepared", lambda *a: audit)
    monkeypatch.setattr(followup, "science_api", lambda *a: api)
    monkeypatch.setattr(followup, "helper", lambda _: SimpleNamespace(selected_case=lambda *a: (plan["case"], [])))
    monkeypatch.setattr(followup.subprocess, "check_output", lambda cmd, **k: (
        followup.SCIENCE if change == "old_source" else binding["science"]) if k["cwd"] == Path(plan["scientific_source"]) else plan["driver_revision"])
    if change:
        with pytest.raises(ValueError): followup.run(path)
        assert not executed
    else:
        result = followup.run(path)
        assert result["status"] == "all_projection_followup_arms_scored" and result["campaign"] == CAMPAIGN
        assert result["source_revision"] == binding["science"]


def test_standalone_api_uses_corrected_science_without_importing_driver_science():
    root = Path(followup.__file__).resolve().parents[1]
    source = root / "benchmarks/work/stockholm_annotation_recovery_pilot_v1_20261019"
    assert source.is_dir(), "Retained corrected scientific checkout is required"
    program = """import json,pathlib,runpy,sys
scope=runpy.run_path(sys.argv[1],run_name='standalone_probe')
source=pathlib.Path(sys.argv[2]).resolve()
api=scope['science_api'](source,pathlib.Path(sys.argv[3]).resolve())
names=('commands','validate_activation','recovery_configuration','native_pair_turnover','control_scores','readback_scores','prepare_partition')
paths={name:str(pathlib.Path(sys.modules[getattr(api,name).__module__].__file__).resolve()) for name in names}
assert all(pathlib.Path(path).is_relative_to(source) for path in paths.values())
from orthohmm.search.selective_hmmer import validate_annotated_alignment
assert pathlib.Path(sys.modules[validate_annotated_alignment.__module__].__file__).resolve().is_relative_to(source)
print(json.dumps(paths))
"""
    env = {**os.environ, "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "1"}
    run = subprocess.run([sys.executable, "-c", program, str(Path(followup.__file__)), str(source), str(root)],
                         cwd=source, env=env, text=True, capture_output=True, timeout=60)
    assert run.returncode == 0, run.stdout + run.stderr
    assert len(json.loads(run.stdout)) == 7


def test_explicit_configuration_does_not_change_legacy_defaults():
    assert followup.configuration()["science"] == followup.SCIENCE
    assert followup.configuration()["preparation_job"] == "24517"
    assert followup.configuration(CAMPAIGN)["preparation_job"] == "24521"
    assert followup.driver_paths()[-1] == followup.PROTOCOL
    assert len(followup.driver_paths(CAMPAIGN)) == 6
    with pytest.raises(ValueError): followup.configuration("stockholm")
