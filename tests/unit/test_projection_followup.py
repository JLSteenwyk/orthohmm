from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

from benchmark_tools import run_projection_followup as followup
from benchmark_tools import run_profile_weighting_followup as weighting
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_profile_weighting_pilot import recovery_configuration
from benchmark_tools.run_membership_hillclimb import write_json
from tests.unit.test_profile_weighting_followup import execution_fixture as weighting_fixture
from tests.unit.test_projection_followup_plan import evidence
from benchmark_tools.plan_projection_followup import selection


def execution_fixture(tmp_path, monkeypatch, *, failed=None, mutation=None, dataset="ThreeKingdoms", passing=None):
    plan, path, executed, scored = weighting_fixture(tmp_path, monkeypatch, failed=failed, dataset=dataset)
    api = weighting.science_api(Path(plan["scientific_source"]))
    original_timed = api.timed_stage
    report, producer = evidence(passing)
    chosen = selection(report, producer)
    plan.update(status="prepared_projection_followup_unrun", scientific_revision=followup.SCIENCE,
        selection=chosen, arms=chosen["execution_arms"], prepared_inputs=plan["pilot_admission"],
        defaults_promoted=False, goal_completed=False,
        hmmer_tools={}, scoring_ids={label: followup.scoring_id(Path(plan["output"]), label)
            for label in chosen["execution_arms"] if label not in ("legacy_off", "projection_hmmer_score_only")})
    monkeypatch.setattr(followup, "UNIVERSES", {"QfO": (3, 1), "ThreeKingdoms": (3, 1)})
    api.recovery_configuration = recovery_configuration
    api.validate_activation = lambda graph, mode, *args, **kwargs: {
        "candidate_digest": mode if mutation == "digest" else "fixed_candidates",
        "completely_measured_queries": 3}
    api.native_pair_turnover = lambda *args: {"files": [], "added_pairs": 0, "lost_pairs": 0}

    def timed(command, folder, stage, science, env):
        if stage == "native":
            result = original_timed(command, folder, stage, science, env)
            if mutation == "tree" and folder.name == "projection_hmmer_score_only":
                (folder / "native/orthohmm_phylogeny/species_tree.rooted.nwk").write_text("changed\n")
            return result
        executed.append((folder.name, stage))
        assert env["NUMBA_NUM_THREADS"] == "16" and env["CUDA_VISIBLE_DEVICES"] == ""
        graph = folder / "graph"
        graph.mkdir()
        seed = graph / "orthogroups_profiles_refined.txt"
        seed.write_text("a b c\n" if mutation != "graph" or folder.name != "projection_hmmer_score_only" else "a b\nc\n")
        write_json(folder / "graph_metrics.json", {"stages": [{"label": "profiles", "output": record(seed)}]})
        return {"returncode": int((folder.name, stage) == failed), "command": command}

    api.timed_stage = timed
    original_partition = api.prepare_partition

    def partition(seed, folder, *args):
        result = original_partition(seed, folder, *args)
        if mutation == "candidate" and folder.parent.name == "legacy_off":
            Path(result["candidate_partition"]["path"]).write_text("a b\nc\n")
            result["candidate_partition"] = record(Path(result["candidate_partition"]["path"]))
        return result

    api.prepare_partition = partition
    original_native = api.replay_native

    def native(row, command):
        result = original_native(row, command)
        result["phylogeny"] = str(Path(result["prediction"]).parents[1] / "orthohmm_phylogeny")
        return result

    api.replay_native = native

    def controls(api, plan, row, native, folder):
        if failed == ("legacy_off", "score"): raise ValueError("Controlled score failure")
        row["scores"] = {"score": .7}
        return row["scores"], {}, []

    api.control_scores = controls
    original_readback = api.readback_scores

    def readback(repo, dataset, row, *args):
        actual, details, pins = original_readback(repo, dataset, row, *args)
        if mutation == "assessment" and row["label"] == "projection_hmmer_score_only":
            actual = {"score": .999}
        return actual, details, pins

    api.readback_scores = readback
    write_json(path, plan)
    followup.validate_plan(plan)
    return api, plan, path, executed, scored


@pytest.mark.parametrize("dataset", ["QfO", "ThreeKingdoms"])
def test_ordered_transfer_reuses_measurement_assessment_without_resampling(tmp_path, monkeypatch, dataset):
    api, plan, path, executed, scored = execution_fixture(tmp_path, monkeypatch, dataset=dataset)
    result = followup.execute_arms(api, plan, path, ["a", "b", "c"], [0, 0, 0], ([], [], []), {})
    assert result["status"] == "all_projection_followup_arms_scored", result["arms"]
    assert executed == [(label, stage) for label in followup.ARMS for stage in ("graph", "native")]
    assert [row["role"] for row in result["arms"]] == [plan["selection"]["roles"][label] for label in followup.ARMS]
    projection, measurement = result["arms"][1], result["arms"][3]
    assert projection["scores"] == measurement["scores"]
    assert projection["details"] == measurement["details"]
    assert projection["metric_vector"] == measurement["metric_vector"]
    assert measurement["measurement_scoring_reuse"] == record(Path(plan["output"]) / "projection_off/assessment_receipt.json")
    assert len(scored) == 3
    if dataset == "QfO": assert scored == list(plan["scoring_ids"].values())
    for row in result["arms"]:
        command = json.loads((Path(plan["output"]) / row["label"] / "preflight.json").read_text())["graph"]
        assert "--profile-sequence-weighting" not in command
        policy, mode = recovery_configuration(row["label"])
        if policy != "legacy": assert command[command.index("--profile-projection-policy") + 1] == policy
        if mode != "off": assert command[command.index("--selective-profile-mode") + 1] == mode


@pytest.mark.parametrize("failure", [("legacy_off", stage) for stage in ("graph", "native", "score")])
def test_failed_legacy_control_blocks_every_dependent_arm(tmp_path, monkeypatch, failure):
    api, plan, path, executed, _ = execution_fixture(tmp_path, monkeypatch, failed=failure)
    result = followup.execute_arms(api, plan, path, ["a", "b", "c"], [0, 0, 0], ([], [], []), {})
    assert result["status"] == "finished_with_explicit_gaps"
    assert result["arms"][0]["status"] == "failed"
    assert all(row["status"] == "not_run_failed_required_control" for row in result["arms"][1:])
    assert all(label == "legacy_off" for label, stage in executed)


@pytest.mark.parametrize("label,blocked", [
    ("projection_off", ["projection_builtin_reselect", "projection_hmmer_score_only", "projection_hmmer_reselect"]),
    ("projection_builtin_reselect", ["projection_hmmer_reselect"]),
    ("projection_hmmer_score_only", ["projection_hmmer_reselect"]),
])
def test_failed_required_controls_cannot_be_bypassed(tmp_path, monkeypatch, label, blocked):
    api, plan, path, executed, _ = execution_fixture(tmp_path, monkeypatch, failed=(label, "graph"))
    result = followup.execute_arms(api, plan, path, ["a", "b", "c"], [0, 0, 0], ([], [], []), {})
    assert all(next(row for row in result["arms"] if row["label"] == name)["status"] == "not_run_failed_required_control" for name in blocked)
    assert not any(name in blocked for name, stage in executed)
    if label == "projection_builtin_reselect":
        assert result["arms"][3]["status"] == followup.SUCCESS


@pytest.mark.parametrize("mutation,failed_label", [("digest", "projection_hmmer_score_only"),
    ("graph", "projection_hmmer_score_only"), ("tree", "projection_hmmer_score_only"),
    ("assessment", "projection_hmmer_score_only"), ("candidate", "legacy_off")])
def test_measurement_and_canonical_identity_checks_are_required(tmp_path, monkeypatch, mutation, failed_label):
    api, plan, path, _, _ = execution_fixture(tmp_path, monkeypatch, mutation=mutation)
    result = followup.execute_arms(api, plan, path, ["a", "b", "c"], [0, 0, 0], ([], [], []), {})
    failed = next(row for row in result["arms"] if row["label"] == failed_label)
    assert failed["status"] == "failed" and failed["error"]
    assert result["status"] == "finished_with_explicit_gaps"


@pytest.mark.parametrize("passing", [["projection_off"], ["projection_builtin_reselect"], ["projection_hmmer_reselect"]])
def test_only_required_original_arm_subsequence_runs_and_failed_pilot_roles_stay_controls(tmp_path, monkeypatch, passing):
    api, plan, path, executed, _ = execution_fixture(tmp_path, monkeypatch, passing=passing)
    result = followup.execute_arms(api, plan, path, ["a", "b", "c"], [0, 0, 0], ([], [], []), {})
    assert result["status"] == "all_projection_followup_arms_scored"
    assert executed == [(label, stage) for label in plan["arms"] for stage in ("graph", "native")]
    assert result["selection"]["advancing_candidates"] == passing
    if passing == ["projection_hmmer_reselect"]:
        assert result["arms"][1]["role"] == "matched_construction_control"
        assert result["arms"][2]["role"] == "matched_decision_control"


@pytest.mark.parametrize("mutation", ["source", "budget", "order", "missing_control", "dependency", "namespace", "fasta", "checkpoint", "unpinned", "promoted"])
def test_changed_execution_contract_is_rejected_before_inference(tmp_path, monkeypatch, mutation):
    _, plan, _, _, _ = execution_fixture(tmp_path, monkeypatch)
    if mutation == "source": plan["scientific_revision"] = "other"
    elif mutation == "budget": plan["cpu"] = 32
    elif mutation == "order": plan["selection"]["execution_arms"].reverse()
    elif mutation == "missing_control": plan["selection"]["execution_arms"].remove("projection_hmmer_score_only")
    elif mutation == "dependency": plan["selection"]["control_dependencies"]["projection_hmmer_reselect"] = []
    elif mutation == "namespace": plan["scoring_ids"]["projection_off"] = "other"
    elif mutation == "fasta": plan["case"]["case"]["fasta_directory"] = "/other"
    elif mutation == "checkpoint": plan["case"]["checkpoint_sha"] = "other"
    elif mutation == "unpinned": plan["inputs"].remove(plan["prepared_inputs"])
    else: plan["defaults_promoted"] = True
    with pytest.raises(ValueError): followup.validate_plan(plan)


def test_no_followup_receipt_and_primary_preparation_are_not_executable(tmp_path):
    with pytest.raises(ValueError): followup.authenticate_prepared(tmp_path, {"status": "no_projection_followup_required"})
    with pytest.raises(FileExistsError): followup.prepare(tmp_path, tmp_path, tmp_path, "QfO")
    with pytest.raises(ValueError, match="frozen"): followup.prepare(tmp_path, tmp_path, tmp_path / "out", "QfO")


def test_standalone_scientific_api_is_bound_to_original_recovery_checkout():
    root = Path(followup.__file__).resolve().parents[1]
    source = root / "benchmarks/work/projection_recovery_pilot_v1_20261019"
    assert source.is_dir(), "Retained scientific checkout is required for this integration test"
    program = """import json,pathlib,runpy,sys
scope=runpy.run_path(sys.argv[1],run_name='standalone_probe')
source=pathlib.Path(sys.argv[2]).resolve()
api=scope['science_api'](source,pathlib.Path(sys.argv[3]).resolve())
names=('commands','validate_activation','recovery_configuration','native_pair_turnover','control_scores','readback_scores','prepare_partition')
paths={name:str(pathlib.Path(sys.modules[getattr(api,name).__module__].__file__).resolve()) for name in names}
assert all(pathlib.Path(path).is_relative_to(source) for path in paths.values())
print(json.dumps(paths))
"""
    env = {**os.environ, "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "1"}
    run = subprocess.run([sys.executable, "-c", program, str(Path(followup.__file__)), str(source), str(root)],
        cwd=source, env=env, text=True, capture_output=True, timeout=60)
    assert run.returncode == 0, run.stdout + run.stderr
    assert len(json.loads(run.stdout)) == 7


def test_mixed_scientific_namespace_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="another source"):
        followup.science_api(tmp_path / "other", Path(followup.__file__).resolve().parents[1])


def authentication_fixture(tmp_path, monkeypatch):
    prepared = {"status": "projection_transfer_inputs_prepared_execution_driver_required",
        "scientific_revision": followup.SCIENCE, "preparer_revision": followup.PREPARER, "cpu": 16,
        "repo": str(tmp_path), "output": str((tmp_path / followup.PREPARED).parent),
        "defaults_promoted": False, "goal_completed": False, "selection": {"execution_arms": list(followup.ARMS)},
        "pilot_readback": {"fixed": "reader"}, "pilot_plan": {"fixed": "plan"},
        "scientific_source": "/frozen/science", "runtime": {"files": []}, "tool_path": "/tools",
        "inputs": [{"path": "all_original_inputs"}], "cases": [{"dataset": "QfO"}, {"dataset": "ThreeKingdoms"}],
        "preparer_files": [{"path": "all_original_source"}], "protocol": {"sha256": "fixed_protocol"}}
    evidence = {key: prepared[key] for key in ("inputs", "cases", "preparer_files", "protocol")}
    audit = {key: deepcopy(prepared[key]) for key in ("selection", "pilot_readback", "pilot_plan", "scientific_source", "runtime", "tool_path")}
    audit.update(scientific_sources=[], hmmer_tools={}, prepared_evidence_digest=hashlib.sha256(
        json.dumps(evidence, sort_keys=True, separators=(",", ":")).encode()).hexdigest())
    calls = []
    def job(*args):
        calls.append(args)
        return {"completed": True}
    monkeypatch.setattr(followup, "helper", lambda _: SimpleNamespace(successful_job=job))
    def process(command, **kwargs):
        if command[0] == "git": return followup.PREPARER
        assert kwargs["cwd"] == tmp_path / followup.PLANNER
        assert kwargs["env"]["PYTHONPATH"] == str(tmp_path / followup.PLANNER)
        assert kwargs["env"]["CUDA_VISIBLE_DEVICES"] == ""
        assert "original.selected_case" in command[2] and "species_tree_inference" in command[2]
        compile(command[2], "frozen_planner_authentication", "exec")
        return json.dumps(audit)
    monkeypatch.setattr(followup.subprocess, "check_output", process)
    return prepared, audit, calls


def test_authentication_revalidates_original_reader_selection_complete_inventory_and_preparation_job(tmp_path, monkeypatch):
    prepared, audit, calls = authentication_fixture(tmp_path, monkeypatch)
    result = followup.authenticate_prepared(tmp_path, prepared)
    assert result == {key: value for key, value in audit.items() if key != "prepared_evidence_digest"}
    assert calls == [("24517", 1, "16G")]


@pytest.mark.parametrize("mutation", ["partial", "source", "preparer", "budget", "selected", "case", "inputs", "protocol", "reader", "runtime", "job"])
def test_changed_preparation_or_unsuccessful_prerequisite_cannot_launch(tmp_path, monkeypatch, mutation):
    prepared, _, _ = authentication_fixture(tmp_path, monkeypatch)
    if mutation == "partial": prepared["status"] = "running"
    elif mutation == "source": prepared["scientific_revision"] = "other"
    elif mutation == "preparer": prepared["preparer_revision"] = "other"
    elif mutation == "budget": prepared["cpu"] = 32
    elif mutation == "selected": prepared["selection"]["execution_arms"].pop()
    elif mutation == "case": prepared["cases"].pop()
    elif mutation == "inputs": prepared["inputs"] = []
    elif mutation == "protocol": prepared["protocol"]["sha256"] = "other"
    elif mutation == "reader": prepared["pilot_readback"] = {"unadmitted": True}
    elif mutation == "runtime": prepared["runtime"] = {"changed": True}
    else:
        monkeypatch.setattr(followup, "helper", lambda _: SimpleNamespace(successful_job=lambda *a: (_ for _ in ()).throw(ValueError("Preparation pending"))))
    with pytest.raises(ValueError): followup.authenticate_prepared(tmp_path, prepared)


@pytest.mark.parametrize("changed_driver", [False, True])
def test_actual_preparation_writes_only_bound_plan_and_rejects_uncommitted_driver(tmp_path, monkeypatch, changed_driver):
    api, original, _, _, _ = execution_fixture(tmp_path, monkeypatch)
    repo = Path(original["repo"])
    driver = tmp_path / "driver"
    (driver / "benchmark_tools").mkdir(parents=True)
    for filename in ("run_projection_followup.py", "run_projection_followup.slurm", "run_profile_weighting_followup.py", "hmm_projection_followup_20261019.json"):
        (driver / "benchmark_tools" / filename).write_text("frozen\n")
    original["case"]["case"]["fastas"] = [record(tmp_path / "input/sp.faa")]
    cases = [deepcopy(original["case"]) for _ in followup.UNIVERSES]
    for case, dataset in zip(cases, followup.UNIVERSES): case["dataset"] = dataset
    prepared = {"cases": cases, "selection": original["selection"], "inputs": original["inputs"],
        "protocol": record(driver / followup.PROTOCOL)}
    prepared_path = repo / followup.PREPARED
    prepared_path.parent.mkdir(parents=True)
    write_json(prepared_path, prepared)
    audit = {key: original[key] for key in ("selection", "pilot_readback", "pilot_plan", "scientific_source", "runtime", "tool_path", "scientific_sources", "hmmer_tools")}
    monkeypatch.setattr(followup, "authenticate_prepared", lambda *args: audit)
    monkeypatch.setattr(followup, "science_api", lambda *args: api)
    monkeypatch.setattr(followup, "helper", lambda _: SimpleNamespace(selected_case=lambda api, repo, dataset: (next(c for c in cases if c["dataset"] == dataset), [])))
    monkeypatch.setattr(followup, "__file__", str(driver / "benchmark_tools/run_projection_followup.py"))
    api.require_profile_runtime = lambda *args: {"tiny_profile": "passed"}
    def git(command, **kwargs):
        if command[1] == "rev-parse": return "frozen_driver_revision"
        return b"changed\n" if changed_driver else (driver / command[-1].split(":", 1)[1]).read_bytes()
    monkeypatch.setattr(followup.subprocess, "check_output", git)
    output = tmp_path / "new_execution"
    if changed_driver:
        with pytest.raises(ValueError, match="committed"): followup.prepare(repo, driver, output, "ThreeKingdoms")
        assert not output.exists()
    else:
        plan = followup.prepare(repo, driver, output, "ThreeKingdoms")
        assert json.loads((output / "execution_plan.json").read_text()) == plan
        assert not (output / "results.json").exists()
        followup.validate_plan(plan)
        with pytest.raises(FileExistsError): followup.prepare(repo, driver, output, "ThreeKingdoms")


@pytest.mark.parametrize("mutation", [None, "runtime", "case", "source", "input", "cuda"])
def test_standalone_run_validates_source_inputs_and_ownership_before_arms(tmp_path, monkeypatch, mutation):
    api, plan, path, executed, _ = execution_fixture(tmp_path, monkeypatch)
    repo = Path(plan["repo"])
    prepared_path = repo / followup.PREPARED
    prepared_path.parent.mkdir(parents=True)
    write_json(prepared_path, {"original": "prepared_receipt"})
    plan["prepared_inputs"] = record(prepared_path)
    plan["inputs"].append(plan["prepared_inputs"])
    write_json(path, plan)
    audit = {key: plan[key] for key in ("selection", "pilot_readback", "pilot_plan", "scientific_source", "runtime", "tool_path", "scientific_sources", "hmmer_tools")}
    monkeypatch.setattr(followup, "authenticate_prepared", lambda *args: audit)
    monkeypatch.setattr(followup, "science_api", lambda *args: api)
    monkeypatch.setattr(followup, "helper", lambda _: SimpleNamespace(selected_case=lambda *args: ({} if mutation == "case" else plan["case"], [])))
    def git(command, **kwargs):
        if kwargs["cwd"] == Path(plan["scientific_source"]): return "other" if mutation == "source" else followup.SCIENCE
        return plan["driver_revision"]
    monkeypatch.setattr(followup.subprocess, "check_output", git)
    if mutation == "runtime": api.runtime = lambda: {"changed": True}
    elif mutation == "input": (tmp_path / "input/sp.faa").write_text(">changed\nAAAA\n")
    elif mutation == "cuda": monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "0")
    if mutation:
        with pytest.raises(ValueError): followup.run(path)
        assert not executed
    else:
        result = followup.run(path)
        assert result["status"] == "all_projection_followup_arms_scored"
        assert result["ownership"]["genes"] == 3 and result["ownership"]["species"] == 1
        with pytest.raises(ValueError, match="prior attempts"): followup.run(path)
