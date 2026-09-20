import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import run_root_context_scaling as module
from tests.unit.test_verify_scaling_allocation import raw, PLAN


def write(path, value):
    path.write_text(json.dumps(value))
    return hashlib.sha256(path.read_bytes()).hexdigest()


def permit(tmp_path):
    policy, ready, auth = [tmp_path / (s + ".json") for s in ("policy", "ready", "auth")]
    psha = write(policy, dict(status="approved_frozen_environment_policy", approval_reference="TEST ONLY"))
    rsha = write(ready, dict(status="environment_preflight_passed", job_id=12345,
        policy_sha256=psha, whole_run_observer_ready=True))
    data = dict(schema="root_context_scaling_authorization_v1", execution_authorized=True,
        plan_sha256=module.PLAN_SHA, recipe_sha256="a" * 64, index=0, job_id=12345,
        environment_policy=dict(path=str(policy), sha256=psha),
        environment_preflight=dict(path=str(ready), sha256=rsha))
    return auth, write(auth, data), data


def test_permit_binds_three_documents(tmp_path):
    path, sha, _ = permit(tmp_path)
    assert len(module.authorization(path, sha, "a" * 64, 0, 12345)) == 3


@pytest.mark.parametrize("key,value", [("execution_authorized", False), ("execution_authorized", 1),
    ("index", 1), ("index", False), ("job_id", 12346), ("job_id", 12345.),
    ("plan_sha256", "b" * 64), ("recipe_sha256", "b" * 64), ("schema", "other")])
def test_wrong_permit_rejected(tmp_path, key, value):
    path, _, data = permit(tmp_path)
    data[key] = value
    with pytest.raises(ValueError):
        module.authorization(path, write(path, data), "a" * 64, 0, 12345)


@pytest.mark.parametrize("which,key,value", [
    ("environment_policy", "status", "unresolved_requires_separate_freeze"),
    ("environment_policy", "approval_reference", ""),
    ("environment_preflight", "status", "failed"),
    ("environment_preflight", "job_id", 12346),
    ("environment_preflight", "whole_run_observer_ready", False),
    ("environment_preflight", "whole_run_observer_ready", 1),
    ("environment_preflight", "policy_sha256", "b" * 64)])
def test_incomplete_environment_rejected(tmp_path, which, key, value):
    path, _, data = permit(tmp_path)
    target = Path(data[which]["path"])
    changed = json.loads(target.read_text())
    changed[key] = value
    data[which]["sha256"] = write(target, changed)
    with pytest.raises(ValueError):
        module.authorization(path, write(path, data), "a" * 64, 0, 12345)


def test_changed_policy_bytes_rejected(tmp_path):
    path, sha, data = permit(tmp_path)
    Path(data["environment_policy"]["path"]).write_text("{}")
    with pytest.raises(ValueError):
        module.authorization(path, sha, "a" * 64, 0, 12345)


def fixture(tmp_path, monkeypatch):
    auth, sha, data = permit(tmp_path)
    recipe = tmp_path / "recipe.json"
    recipe.write_text("{}")
    task = dict(index=0, run=dict(measurement_directory=str(tmp_path / "native/measurement")))
    monkeypatch.setattr(module, "OUTPUT_ROOT", tmp_path / "output")
    monkeypatch.setattr(module, "select", lambda *args: ({}, task))
    monkeypatch.setattr(module, "preflight", lambda plan: 12345)
    monkeypatch.setattr(module.sys, "pycache_prefix", str(tmp_path / "cache"))
    monkeypatch.setattr(module.subprocess, "run", lambda *args, **kwargs:
        SimpleNamespace(returncode=0, stdout=raw(), stderr=""))
    calls = []

    def measure(*args):
        calls.append(args)
        (tmp_path / "native").mkdir()
        write(tmp_path / "native/verification.json", dict(status="command_exited_zero"))
        return dict(status="command_exited_zero")

    monkeypatch.setattr(module, "measure_task", measure)
    args = (PLAN, 0, recipe, "a" * 64, auth, sha)
    return args, calls, data


def test_one_task_and_no_repeat(tmp_path, monkeypatch):
    args, calls, _ = fixture(tmp_path, monkeypatch)
    result = module.execute(*args)
    assert len(calls) == 1
    assert calls[0][-1] == 12345
    assert result["status"] == "measurement_returned_pending_audit"
    for key in ("scientific_timings_admitted", "native_outputs_validated", "scheduler_terminal_verified",
                "environmental_validity_established", "next_submission_authorized", "automatic_retry"):
        assert result[key] is False
    with pytest.raises(FileExistsError):
        module.execute(*args)
    assert len(calls) == 1


@pytest.mark.parametrize("returncode,stdout", [(1, "error"), (0, raw(NumCPUs="32")),
    (0, raw(JobState="PENDING")), (0, raw(JobState="COMPLETED"))])
def test_controller_failure_prevents_inference(tmp_path, monkeypatch, returncode, stdout):
    args, calls, _ = fixture(tmp_path, monkeypatch)
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k:
        SimpleNamespace(returncode=returncode, stdout=stdout, stderr=""))
    with pytest.raises(ValueError):
        module.execute(*args)
    assert not calls
    assert json.loads((tmp_path / "output/sessions/task_00/result.json").read_text())["status"] == "executor_failed"


@pytest.mark.parametrize("exception", [RuntimeError, KeyboardInterrupt])
def test_measurement_failure_saved_without_retry(tmp_path, monkeypatch, exception):
    args, calls, _ = fixture(tmp_path, monkeypatch)

    def fail(*a):
        calls.append(a)
        raise exception("test failure")

    monkeypatch.setattr(module, "measure_task", fail)
    with pytest.raises(exception):
        module.execute(*args)
    assert len(calls) == 1
    saved = json.loads((tmp_path / "output/sessions/task_00/result.json").read_text())
    assert saved["error_type"] == exception.__name__
    assert saved["automatic_retry"] is False


def test_unauthorized_does_not_create_receipts(tmp_path, monkeypatch):
    args, calls, data = fixture(tmp_path, monkeypatch)
    data["execution_authorized"] = False
    args = (*args[:-1], write(args[-2], data))
    with pytest.raises(ValueError):
        module.execute(*args)
    assert not calls and not (tmp_path / "output").exists()


def test_current_checkout_is_not_a_deployed_recipe():
    with pytest.raises(ValueError, match="deployed"):
        module.select(PLAN, 0, module.RECIPE_PATH, "a" * 64)


def test_frozen_plan_does_not_itself_authorize():
    plan, _, _ = module.load_task(PLAN, 0)
    assert plan["execution_authorized"] is False
    assert plan["environment_policy"]["service_change_authorized"] is False


def test_select_checks_deployed_source_bytes(tmp_path, monkeypatch):
    root = tmp_path / "recipe"
    folder = root / "benchmark_tools"
    folder.mkdir(parents=True)
    source, script, plan = [folder / name for name in ("runner.py", "run.sh", "plan.json")]
    source.write_text("# constructed source inventory\n")
    script.write_text("#!/bin/bash\n")
    plan.write_bytes(PLAN.read_bytes())
    recipe = tmp_path / "recipe.json"
    sha = write(recipe, dict(roots=[str(root)], records=[dict(module.record(p), kind="file")
        for p in (source, script, plan)]))
    monkeypatch.setattr(module, "__file__", str(source))
    monkeypatch.setattr(module, "RECIPE_ROOT", root)
    monkeypatch.setattr(module, "RECIPE_PATH", recipe)
    monkeypatch.setattr(module, "SUBMISSION_SCRIPT", script)
    assert module.select(plan, 0, recipe, sha)[1]["index"] == 0
    source.write_text("# changed source\n")
    with pytest.raises(ValueError, match="not pinned"):
        module.select(plan, 0, recipe, sha)


def test_changed_permit_after_native_preserves_failure(tmp_path, monkeypatch):
    args, calls, data = fixture(tmp_path, monkeypatch)
    original = module.measure_task

    def changed(*a):
        result = original(*a)
        Path(data["environment_preflight"]["path"]).write_text("{}")
        return result

    monkeypatch.setattr(module, "measure_task", changed)
    with pytest.raises(ValueError):
        module.execute(*args)
    assert len(calls) == 1
    saved = json.loads((tmp_path / "output/sessions/task_00/result.json").read_text())
    assert saved["status"] == "executor_failed"
    assert saved["wrapper_status"] == "command_exited_zero"
    assert saved["verification"]["sha256"]


@pytest.mark.parametrize("status,exitcode", [("command_exited_zero", 0),
    ("command_exited_nonzero", 1), ("runtime_changed_or_unverifiable", 1)])
def test_cli_exit_is_not_an_admission(tmp_path, monkeypatch, status, exitcode):
    monkeypatch.setattr(module, "execute", lambda *args: dict(wrapper_status=status))
    assert module.main(["--plan", str(PLAN), "--recipe", str(tmp_path / "recipe.json"),
        "--index", "0", "--recipe-sha", "a" * 64, "--authorization", str(tmp_path / "auth.json"),
        "--authorization-sha", "b" * 64]) == exitcode
