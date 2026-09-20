from copy import deepcopy
import gzip
import json
from pathlib import Path
import shutil

import pytest

from benchmark_tools import audit_root_context_overhead as module
from tests.unit.test_verify_root_context_overhead_provenance import fixture as task_fixture, RESULTS


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


def panel_fixture(tmp_path, failed=None):
    context, _, _, _, scheduler, job = task_fixture(tmp_path)
    directory = tmp_path / module.OUTPUT
    launch = dict(job_id=job, plan_sha256=module.PLAN_SHA, recipe_sha256=context["recipe_sha256"],
        protocol_sha256=module.PROTOCOL_SHA, executable=context["plan"]["launcher_python"],
        source_directory=str(module.RECIPE_ROOT), scientific_timings_admitted=False,
        affinity=list(range(20)), uname=["Linux", "spark-7ff0", "test", "test", "aarch64"])
    save(directory / "launch.json", launch)
    rows = []
    for index, task in enumerate(context["plan"]["runs"]):
        row = {key: task[key] for key in module.IDENTITY}
        row.update(task=deepcopy(task), job_id=job, scientific_timings_admitted=False, status="measurement_completed")
        if index == failed:
            row.update(status="failed", error_type="ValueError", error="test failure")
        elif failed is not None and index > failed:
            row.update(status="not_run_after_failure", failed_index=failed)
        rows.append(row)
        save(directory / f"task_{index:02d}.json", row)
        save(directory / f"progress_{index:02d}.json", dict(tasks=deepcopy(rows), job_id=job,
            failed_index=failed if failed is not None and index >= failed else None))
    panel = dict(status="root_context_overhead_completed" if failed is None else "root_context_overhead_stopped",
        tasks=rows, job_id=job, plan_sha256=module.PLAN_SHA, recipe_sha256=context["recipe_sha256"],
        scientific_timings_admitted=False, native_outputs_validated=False, publication_ready=False)
    save(directory / "result.json", panel)
    if failed is not None:
        scheduler = scheduler.replace("JobState=COMPLETED", "JobState=FAILED").replace("ExitCode=0:0", "ExitCode=1:0")
    return directory, context, job, scheduler


@pytest.mark.parametrize("failed", [None, 0, 7, 17])
def test_fixed_panel_outcomes_bind(tmp_path, failed):
    directory, context, job, raw = panel_fixture(tmp_path, failed)
    allocation = module.scheduler_identity(raw, job, require_completed=False)
    assert len(module.bind_panel(directory, context, job, allocation)["tasks"]) == 18


@pytest.mark.parametrize("fault", ["checkpoint", "receipt", "pair", "arm", "missing", "extra", "unrun", "launch", "scheduler"])
def test_panel_drift_rejected(tmp_path, fault):
    directory, context, job, raw = panel_fixture(tmp_path, 7)
    allocation = module.scheduler_identity(raw, job, require_completed=False)
    if fault == "checkpoint":
        path = directory / "progress_07.json"
        value = module.read(path)
        value["failed_index"] = None
        save(path, value)
    elif fault in ("receipt", "pair", "arm"):
        path = directory / "result.json" if fault != "receipt" else directory / "task_00.json"
        value = module.read(path)
        row = value if fault == "receipt" else value["tasks"][0]
        row[{"receipt": "job_id", "pair": "pair", "arm": "arm"}[fault]] = "wrong"
        save(path, value)
    elif fault == "missing":
        (directory / "task_17.json").unlink()
    elif fault in ("extra", "unrun"):
        (directory / ("run_18" if fault == "extra" else "run_08")).mkdir()
    elif fault == "launch":
        path = directory / "launch.json"
        value = module.read(path)
        value["affinity"][0] = True
        save(path, value)
    else:
        allocation["ExitCode"] = "0:0"
    with pytest.raises((ValueError, OSError)):
        module.bind_panel(directory, context, job, allocation)


@pytest.mark.parametrize("index", [0, 1, 4, 5])
@pytest.mark.parametrize("equivalent", [True, False])
def test_real_binding_dispatches_correct_replay_and_native_validation(tmp_path, monkeypatch, index, equivalent):
    context, _, original, receipt, scheduler, job = task_fixture(tmp_path, index)
    task = context["plan"]["runs"][index]
    directory = tmp_path / Path(task["run"]["measurement_directory"]).parent.relative_to(module.ROOT)
    shutil.copytree(original, directory)
    points = [dict(host=[dict(started_monotonic_ns=1), dict(finished_monotonic_ns=19)], ticks=100,
        lineage=dict(boot_before="boot", identities_before={"/": [1, 2], "/system.slice": [1, 3]}),
        root_context=dict(identities_before={"/": [1, 2]}, host_after=dict(finished_ns=20)))]
    lineage = dict(measured=dict(points=points), native_wall_s=15., memory=dict(finished_ns=25),
        screening=dict(observation_window=dict(native_pressure={"retained": True})),
        original_flagged_intervals=[0], narrow_flagged_intervals=[0], evidence=[])
    calls = []

    def replay(path, actual_job, command, **kwargs):
        calls.append("root" if task["arm"] == "root_context" else "lineage")
        assert actual_job == job and kwargs == dict(expected_timeout_s=900)
        assert command == module.read(directory / "preparation.json")["measured_argv"]
        if task["arm"] == "root_context":
            return dict(lineage=lineage, context={"root": True}, evidence=[])
        return lineage

    monkeypatch.setattr(module, "replay_root" if task["arm"] == "root_context" else "replay_lineage", replay)
    monkeypatch.setattr(module, "replay_lineage" if task["arm"] == "root_context" else "replay_root",
                        lambda *a, **k: pytest.fail("wrong collector replay"))

    def validate(run, adapted, roots):
        calls.append("validate")
        assert roots == {module.ROOT: tmp_path}
        assert run == module.read(directory / "preparation.json")["run"]
        return dict(input_genes=73266, checked_files=[], gnu_time_companion=dict(
            source=module.record(directory / "verification.json"), accounting={}))

    monkeypatch.setattr(module, "validate", validate)
    monkeypatch.setattr(module, "fingerprint", lambda *a: dict(identity={"groups": 1}, evidence=[]))
    prior = dict(runs=[dict(index=task["native_parent_index"], status="validated", method=task["method"],
                            work_identity={"groups": 1 if equivalent else 2})])
    result = module.successful_task(tmp_path, context, task, receipt, scheduler, job, prior)
    assert calls == (["root", "validate"] if task["arm"] == "root_context" else ["lineage", "validate"])
    assert result["status"] == ("validated" if equivalent else "output_mismatch")
    assert result["observation_finish_ns"] == 25
    assert result["narrow_flagged_intervals"] == [0]
    assert result["root_context"] == ({"root": True} if task["arm"] == "root_context" else None)
    assert result["pressure_observation_window"] == {"retained": True}


@pytest.mark.parametrize("case", ["complete", "failure", "incomplete", "raw_invalid", "overlap", "scope", "root_scope", "session_invalid"])
def test_whole_audit_retains_all_pairs_and_nulls_invalid_conclusions(tmp_path, monkeypatch, case):
    directory, context, job, raw = panel_fixture(tmp_path, 7 if case == "failure" else None)
    scheduler_path = tmp_path / "scheduler.txt"
    scheduler_path.write_text(raw)
    prior_path = tmp_path / "prior.gz"
    prior_path.write_bytes(gzip.compress(json.dumps(dict(all_tasks_validated=True, issues=[], runs=[{}, {}, {}], inventory=[])).encode()))
    monkeypatch.setattr(module, "PRIOR_AUDIT_SHA", module.record(prior_path)["sha256"])
    for name in ("dgx_root_context_overhead_plan_20260919.json", "dgx_root_context_native_plan_20260919.json",
                 "ROOT_CONTEXT_OVERHEAD_PROTOCOL_20260919.md"):
        target = str(module.RECIPE_ROOT / "benchmark_tools/results" / name)
        context["recipe"]["records"] = [r for r in context["recipe"]["records"] if r["path"] != target]
        context["recipe"]["records"].append(dict(module.record(RESULTS / name), path=target, kind="file"))
    monkeypatch.setattr(module, "load_context", lambda *a: context)
    monkeypatch.setattr(module, "recipe_evidence", lambda *a: [])
    def session(directory, scheduler, recipe, sha, actual_job, timezone):
        assert directory == tmp_path / "session" and scheduler == scheduler_path
        assert actual_job == job and timezone == "America/New_York"
        if case == "session_invalid":
            raise ValueError("Waiting receipt incomplete")
        return dict(status="overhead_bounded_session_verified", evidence=[])
    monkeypatch.setattr(module, "audit_session", session)
    called = []

    def success(archive, context, task, *args):
        index = task["index"]
        called.append(index)
        if case == "raw_invalid" and index == 1:
            raise ValueError("raw replay fails")
        row = {key: task[key] for key in module.IDENTITY}
        row.update(status="validated", output_equivalent=True, work_identity={"groups": task["method"]},
            native_wall_s=100., scientific_timings_admitted=False,
            common_scope_identity="changed" if case == "scope" and index == 1 else "same",
            root_scope_identity="changed" if case == "root_scope" and index == 1 else "same",
            observation_start_ns=0 if case == "overlap" and index == 1 else index*10,
            observation_finish_ns=index*10+9)
        return row

    monkeypatch.setattr(module, "successful_task", success)
    if case == "incomplete":
        (directory / "result.json").unlink()
    result = module.audit(tmp_path, RESULTS, tmp_path / "recipe.json", context["recipe_sha256"], scheduler_path, job, prior_path,
                          tmp_path / "session", "America/New_York")
    assert len(result["runs"]) == 18 and len(result["comparison"]["pairs"]) == 9
    assert result["scientific_timings_admitted"] is False
    assert result["comparison"]["engineering_budget_passed"] is (True if case == "complete" else None)
    if case == "failure":
        assert called == list(range(7))
        assert result["runs"][7]["status"] == "retained_failure"
        assert all(r["status"] == "retained_unrun" for r in result["runs"][8:])
    elif case == "incomplete":
        assert called == []
    elif case == "raw_invalid":
        assert result["runs"][1]["status"] == "failed_or_invalid"
    elif case == "session_invalid":
        assert result["validated_tasks"] == 18
        assert result["issues"][0]["stage"] == "waiting_session"
        assert result["waiting_session"]["status"] == "failed_or_invalid"


def test_nonterminal_scheduler_rejected_before_archive_read(tmp_path, monkeypatch):
    _, _, job, raw = panel_fixture(tmp_path)
    path = tmp_path / "scheduler.txt"
    path.write_text(raw.replace("JobState=COMPLETED", "JobState=RUNNING"))
    monkeypatch.setattr(module, "inventory", lambda *a: pytest.fail("premature archive read"))
    with pytest.raises(ValueError, match="terminal"):
        module.audit(tmp_path, RESULTS, tmp_path / "missing", "sha", path, job, tmp_path / "missing",
                     tmp_path / "session", "America/New_York")
