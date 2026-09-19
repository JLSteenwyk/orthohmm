from copy import deepcopy
import gzip
import json
from pathlib import Path

import pytest

from benchmark_tools import audit_root_context_native as module
from tests.unit.test_verify_root_context_native_provenance import fixture as task_fixture
from tests.unit.test_verify_root_context_native_provenance import RESULTS


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


def panel_fixture(tmp_path, failed=None):
    args = task_fixture(tmp_path)
    context, job, scheduler = args[0], args[-1], args[-2]
    directory = tmp_path / module.OUTPUT
    launch = dict(job_id=job, plan_sha256=module.PLAN_SHA, recipe_sha256=context["recipe_sha256"],
        protocol_sha256=module.PROTOCOL_SHA, executable=context["plan"]["launcher_python"],
        source_directory=str(module.RECIPE_ROOT), scientific_timings_admitted=False,
        affinity=list(range(20)), uname=["Linux", "spark-7ff0", "test", "test", "aarch64"])
    save(directory / "launch.json", launch)
    rows = []
    for index, task in enumerate(context["plan"]["runs"]):
        row = dict(index=index, task=deepcopy(task), job_id=job, scientific_timings_admitted=False,
                   status="measurement_completed")
        if failed is not None and index == failed:
            row.update(status="failed", error_type="ValueError", error="test native failure")
        elif failed is not None and index > failed:
            row.update(status="not_run_after_failure", failed_index=failed)
        rows.append(row)
        save(directory / f"task_{index:02d}.json", row)
        save(directory / f"progress_{index:02d}.json", dict(tasks=deepcopy(rows), job_id=job,
             failed_index=failed if failed is not None and index >= failed else None))
    panel = dict(status="root_context_native_completed" if failed is None else "root_context_native_stopped",
        job_id=job, plan_sha256=module.PLAN_SHA, recipe_sha256=context["recipe_sha256"], tasks=rows,
        scientific_timings_admitted=False, native_outputs_validated=False, publication_ready=False)
    save(directory / "result.json", panel)
    if failed is not None:
        scheduler = scheduler.replace("JobState=COMPLETED", "JobState=FAILED").replace("ExitCode=0:0", "ExitCode=1:0")
    return directory, context, job, scheduler


@pytest.mark.parametrize("failed", [None, 0, 1, 2])
def test_panel_preserves_all_outcomes(tmp_path, failed):
    args = panel_fixture(tmp_path, failed)
    result = module.bind_panel(*args)
    assert len(result["tasks"]) == 3
    assert result["publication_ready"] is False


@pytest.mark.parametrize("fault", ["checkpoint", "receipt", "launch", "affinity", "host", "admission",
    "task_order", "task_count", "unrun_artifacts", "after_failure", "scheduler", "missing", "extra"])
def test_panel_tampering(tmp_path, fault):
    args = panel_fixture(tmp_path, 1)
    directory = args[0]
    name = "result.json"
    if fault in ("launch", "affinity", "host"):
        name = "launch.json"
    elif fault == "checkpoint":
        name = "progress_01.json"
    elif fault == "receipt":
        name = "task_00.json"
    value = module.read(directory / name)
    if fault == "checkpoint":
        value["failed_index"] = None
    elif fault == "receipt":
        value["job_id"] += 1
    elif fault == "launch":
        value["recipe_sha256"] = "0" * 64
    elif fault == "affinity":
        value["affinity"][0] = True
    elif fault == "host":
        value["uname"][1] = "bizon"
    elif fault == "admission":
        value["scientific_timings_admitted"] = 0
    elif fault == "task_order":
        value["tasks"].reverse()
    elif fault == "task_count":
        value["tasks"].pop()
    elif fault == "after_failure":
        value["tasks"][2]["status"] = "measurement_completed"
    elif fault == "unrun_artifacts":
        (directory / "run_02").mkdir()
    elif fault == "extra":
        (directory / "run_03").mkdir()
    elif fault == "scheduler":
        args = (*args[:-1], args[-1].replace("ExitCode=1:0", "ExitCode=0:0"))
    save(directory / name, value)
    if fault == "missing":
        (directory / "progress_02.json").unlink()
    with pytest.raises((ValueError, OSError)):
        module.bind_panel(*args)


def test_terminal_gate_precedes_archive_reads(tmp_path, monkeypatch):
    args = panel_fixture(tmp_path)
    scheduler = tmp_path / "scheduler.txt"
    scheduler.write_text(args[-1].replace("JobState=COMPLETED", "JobState=RUNNING"))
    monkeypatch.setattr(module, "inventory", lambda *a: pytest.fail("archive read before terminal"))
    with pytest.raises(ValueError, match="terminal"):
        module.audit(tmp_path, tmp_path, tmp_path / "missing", "sha", scheduler, args[2], tmp_path / "missing")


@pytest.mark.parametrize("equivalent", [True, False])
def test_successful_task_replays_and_checks_native_products(tmp_path, monkeypatch, equivalent):
    args = task_fixture(tmp_path, 2)
    context, index, prep, verification, receipt, measured, scheduler, job = args
    task = context["plan"]["runs"][index]
    directory = tmp_path / Path(task["run"]["measurement_directory"]).parent.relative_to(module.ROOT)
    save(directory / "preparation.json", prep)
    save(directory / "measurement/lineage_report.json", measured)
    (directory / "verification.json").write_bytes(verification.read_bytes())
    points = [dict(host=[dict(started_monotonic_ns=1)], root_context=dict(boot_before="boot", ticks=100,
        identities_before={"root": "root"}, host_after=dict(finished_ns=9)))]
    lineage = dict(measured=dict(points=points), native_wall_s=15., memory={}, screening={},
                   original_flagged_intervals=[1], narrow_flagged_intervals=[1])
    calls = []

    def replay(path, actual_job, command, **kwargs):
        calls.append("replay")
        assert actual_job == job and command == prep["measured_argv"]
        assert kwargs == dict(expected_timeout_s=900)
        return dict(lineage=lineage, context={"retained": True}, evidence=[])

    def validate(run, adapted, roots):
        calls.append("validate")
        assert run == prep["run"] and adapted["command"] == prep["measured_argv"]
        assert roots == {module.ROOT: tmp_path}
        return dict(checked_files=[], input_genes=73266,
                    gnu_time_companion=dict(source=module.record(verification), accounting={}))

    monkeypatch.setattr(module, "replay", replay)
    monkeypatch.setattr(module, "validate", validate)
    monkeypatch.setattr(module, "fingerprint", lambda *a: dict(identity={"groups": 1}, evidence=[]))
    prior = dict(runs=[dict(index=2, status="validated", method=task["method"],
        work_identity={"groups": 1 if equivalent else 2}, inventory=[], checked_evidence=[])])
    result = module.successful_task(tmp_path, context, task, receipt, scheduler, job, prior)
    assert calls == ["replay", "validate"]
    assert result["status"] == ("validated" if equivalent else "output_mismatch")
    assert result["narrow_flagged_intervals"] == [1]
    assert result["root_context"] == {"retained": True}
    assert result["scientific_timings_admitted"] is False


def test_completed_task_may_belong_to_later_failed_panel(tmp_path):
    from benchmark_tools import verify_root_context_native_provenance as verifier
    args = task_fixture(tmp_path)
    args[-2] = args[-2].replace("JobState=COMPLETED", "JobState=FAILED").replace("ExitCode=0:0", "ExitCode=1:0")
    with pytest.raises(ValueError):
        verifier.verify(*args)
    assert verifier.verify(*args, require_completed=False)["status"] == "root_context_native_provenance_bound"


@pytest.mark.parametrize("case", ["complete", "failure", "incomplete", "raw_invalid", "overlap", "scope_change"])
def test_whole_audit_retains_every_task(tmp_path, monkeypatch, case):
    directory, context, job, scheduler = panel_fixture(tmp_path, 1 if case == "failure" else None)
    scheduler_path = tmp_path / "scheduler.txt"
    scheduler_path.write_text(scheduler)
    prior_path = tmp_path / "prior.json.gz"
    prior_path.write_bytes(gzip.compress(json.dumps(dict(runs=[])).encode()))
    monkeypatch.setattr(module, "PRIOR_SHA", module.record(prior_path)["sha256"])
    for name in ("dgx_root_context_native_plan_20260919.json", "dgx_lineage_native_plan_20260919.json",
                 "ROOT_CONTEXT_NATIVE_PROTOCOL_20260919.md"):
        target = str(module.RECIPE_ROOT / "benchmark_tools/results" / name)
        context["recipe"]["records"] = [r for r in context["recipe"]["records"] if r["path"] != target]
        context["recipe"]["records"].append(dict(module.record(RESULTS / name), path=target, kind="file"))
    monkeypatch.setattr(module, "load_context", lambda *a: context)
    monkeypatch.setattr(module, "recipe_evidence", lambda *a: [])
    called = []

    def task(*args):
        index = args[2]["index"]
        called.append(index)
        if case == "raw_invalid" and index == 1:
            raise ValueError("raw replay fails")
        return dict(index=index, status="validated", output_equivalent=True,
            scope_identity="other" if case == "scope_change" and index == 1 else "same",
            observation_start_ns=0 if case == "overlap" and index == 1 else 10*index,
            observation_finish_ns=10*index+9)

    monkeypatch.setattr(module, "successful_task", task)
    if case == "incomplete":
        (directory / "result.json").unlink()
    result = module.audit(tmp_path, RESULTS, tmp_path / "recipe.json", context["recipe_sha256"],
                          scheduler_path, job, prior_path)
    assert len(result["runs"]) == 3 and result["publication_ready"] is False
    assert result["scientific_timings_admitted"] is False
    assert result["all_tasks_validated"] == (case == "complete")
    if case == "failure":
        assert called == [0]
        assert [r["status"] for r in result["runs"]] == ["validated", "retained_failure", "retained_unrun"]
        assert result["all_outputs_equivalent"] is None
    elif case == "incomplete":
        assert called == []
        assert len(result["issues"]) == 1
    elif case == "raw_invalid":
        assert called == [0, 1, 2]
        assert result["runs"][1]["status"] == "failed_or_invalid"
        assert result["all_outputs_equivalent"] is None
    elif case in ("overlap", "scope_change"):
        assert result["issues"][0]["stage"] == "temporal_identity"
