"""Request composition tests; scheduler and scientific admission are not simulated claims."""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import prepare_native_factorial_request as handoff
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def write(path, value):
    path.write_text(json.dumps(value))
    return record(path)


def held(job):
    values = dict(JobId=str(job), JobName="orthohmm_factorial_cost", JobState="PENDING",
                  Reason="JobHeldUser", Partition="gpu", ReqNodeList="bizon", NumCPUs="64",
                  NumTasks="1", MinMemoryNode="128G", Requeue="0", Restarts="0",
                  TimeLimit="1-02:00:00", Command=str(handoff.SCRIPT), WorkDir=str(handoff.ROOT))
    values["CPUs/Task"] = "64"
    return " ".join(f"{key}={value}" for key, value in values.items())


@pytest.fixture
def setup(tmp_path, monkeypatch):
    monkeypatch.setattr(handoff, "ROOT", tmp_path)
    runs = [{"output_root": str(tmp_path / f"run_{i:02d}")} for i in range(13)]
    plan = {"runs": runs, "helper_sources": [], "evidence": [], "panel_root": str(tmp_path)}
    plan_ref = write(tmp_path / "plan.json", plan)
    policy_ref = write(tmp_path / "policy.json", {"plan_sha256": plan_ref["sha256"], "evidence": []})
    preparation_ref = write(tmp_path / "preparation.json", {"plan": plan_ref, "policy": policy_ref})
    monkeypatch.setattr(handoff, "PREPARATION", Path(preparation_ref["path"]))
    monkeypatch.setattr(handoff, "PREPARATION_SHA256", preparation_ref["sha256"])
    history = []
    for i in range(3):
        history.append(write(tmp_path / f"review_{i}.json", {
            "job_id": 22427+i, "index": i, "cell": f"cell{i}", "plan": plan_ref,
            "dataset": "orthobench", "status": "native_success"}))
    score = {"schema": "native_factorial_orthobench_score_v1", "status": "terminal_native_orthobench_scored",
             "terminal_review": history[-1], "index": 2, "cell": "cell2", "job_id": 22429,
             "plan": plan_ref, "native_outputs_validated": True, "accuracy_evaluated": True}
    score_ref = write(tmp_path / "score.json", score)
    calls = []
    monkeypatch.setattr(handoff, "validate_plan", lambda p: p["runs"])
    monkeypatch.setattr(handoff, "validate_request", lambda *args: calls.append("validate_request"))
    monkeypatch.setattr(handoff, "reviewed_history", lambda *args: calls.append("reviewed_history") or [{"fresh": True}])
    monkeypatch.setattr(handoff, "available_memory", lambda raw: handoff.MEMORY)
    monkeypatch.setattr(handoff, "snapshot", lambda: {"diagnostic": True})
    monkeypatch.setattr(handoff.subprocess, "check_output", lambda *args, **kw: "prepared-source\n")
    def scheduler(command, **kwargs):
        calls.append(command)
        assert command[:3] == ["scontrol", "show", "job"]
        return SimpleNamespace(stdout=held(int(command[3])), stderr="")
    monkeypatch.setattr(handoff.subprocess, "run", scheduler)
    return tmp_path, history, score_ref, calls


def test_prepare_only_writes_request_no_submission_or_release(setup):
    tmp, history, score_ref, calls = setup
    ref = handoff.prepare(22430, history, score_ref, tmp / "request.json")
    request = json.loads(Path(ref["path"]).read_text())
    assert request["index"] == 3 and request["history"] == history
    assert request["prior_history_scheduler_checks"] == [{"fresh": True}]
    assert request["automatic_retry"] is request["timing_success_established"] is False
    assert request["capacity_precheck"]["background_cpu_used_for_eligibility"] is False
    assert "reviewed_history" in calls and "validate_request" in calls
    assert all(call[:3] == ["scontrol", "show", "job"] for call in calls if isinstance(call, list))
    with pytest.raises(ValueError, match="fresh direct"):
        handoff.prepare(22430, history, score_ref, tmp / "request.json")


@pytest.mark.parametrize("key,value", [("JobState", "RUNNING"), ("Reason", "Resources"),
    ("JobId", "999"), ("NumCPUs", "32"), ("NumTasks", "2"), ("MinMemoryNode", "192G"),
    ("Requeue", "1"), ("Restarts", "1"), ("Partition", "cpu"), ("ReqNodeList", "other"),
    ("Command", "/wrong"), ("WorkDir", "/wrong"), ("CPUs/Task", "32")])
def test_wrong_held_job_refused(setup, key, value):
    raw = held(22430)
    fields = dict(item.split("=", 1) for item in raw.split())
    fields[key] = value
    with pytest.raises(ValueError, match="Held job"):
        handoff.held_job(" ".join(f"{k}={v}" for k,v in fields.items()), 22430)


def test_duplicate_controller_key_refused(setup):
    with pytest.raises(ValueError, match="Held job"):
        handoff.held_job(held(22430) + " JobId=22430", 22430)


@pytest.mark.parametrize("key,value", [("status", "queued"), ("index", 1), ("cell", "other"),
    ("job_id", 999), ("native_outputs_validated", False), ("accuracy_evaluated", False),
    ("terminal_review", {})])
def test_wrong_previous_score_refused(setup, key, value):
    tmp, history, score_ref, _ = setup
    score = json.loads(Path(score_ref["path"]).read_text())
    score[key] = value
    score_ref = write(tmp / "changed_score.json", score)
    with pytest.raises(ValueError, match="bound separate score"):
        handoff.prepare(22430, history, score_ref, tmp / "request.json")
    assert not (tmp / "request.json").exists()


def test_unresolved_history_aborts_before_scheduler_or_file(setup, monkeypatch):
    tmp, history, score_ref, calls = setup
    def unresolved(*args):
        raise ValueError("Previous factorial identity unresolved")
    monkeypatch.setattr(handoff, "reviewed_history", unresolved)
    with pytest.raises(ValueError, match="unresolved"):
        handoff.prepare(22430, history, score_ref, tmp / "request.json")
    assert calls == [] and not (tmp / "request.json").exists()


def test_missing_score_and_old_job_refused(setup):
    tmp, history, score_ref, _ = setup
    with pytest.raises(ValueError, match="preceding OrthoBench score"):
        handoff.prepare(22430, history, None, tmp / "request.json")
    with pytest.raises(ValueError, match="new held job"):
        handoff.prepare(22429, history, score_ref, tmp / "request.json")


def test_capacity_is_real_gate_not_background_cpu(setup, monkeypatch):
    tmp, history, score_ref, _ = setup
    monkeypatch.setattr(handoff, "available_memory", lambda raw: handoff.MEMORY - 1)
    with pytest.raises(ValueError, match="Unsafe"):
        handoff.prepare(22430, history, score_ref, tmp / "request.json")


@pytest.mark.parametrize("directory", ["run_03", "sessions/run_03"])
def test_previous_output_never_implicitly_resumed(setup, directory):
    tmp, history, score_ref, _ = setup
    (tmp / directory).mkdir(parents=True)
    with pytest.raises(ValueError, match="no implicit retry"):
        handoff.prepare(22430, history, score_ref, tmp / "request.json")


def test_incomplete_or_complete_plan_has_no_next_identity(setup):
    tmp, history, score_ref, _ = setup
    with pytest.raises(ValueError, match="remaining reviewed-prefix"):
        handoff.prepare(22430, history[:2], score_ref, tmp / "request.json")
    with pytest.raises(ValueError, match="remaining reviewed-prefix"):
        handoff.prepare(22430, history * 4 + history[:1], score_ref, tmp / "request.json")


def test_failed_or_qfo_prior_does_not_borrow_orthobench_score(setup):
    tmp, history, score_ref, _ = setup
    prior = json.loads(Path(history[-1]["path"]).read_text())
    prior["dataset"] = "qfo_corrected"
    history[-1] = write(tmp / "qfo_review.json", prior)
    with pytest.raises(ValueError, match="Do not substitute"):
        handoff.prepare(22430, history, score_ref, tmp / "request.json")
    ref = handoff.prepare(22430, history, None, tmp / "request.json")
    assert json.loads(Path(ref["path"]).read_text())["index"] == 3


def test_explicit_pin_mismatch_refused(setup):
    _, history, _, _ = setup
    with pytest.raises(ValueError, match="checksum"):
        handoff.pinned(Path(history[0]["path"]), "0" * 64)
