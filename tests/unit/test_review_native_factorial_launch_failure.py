"""Synthetic refusals do not execute inference or submit scheduler jobs."""

from copy import deepcopy
import json
from pathlib import Path
import subprocess

import pytest

from benchmark_tools import review_native_factorial_launch_failure as module


@pytest.fixture
def case(tmp_path):
    root = tmp_path / "run"
    request = dict(job_id=42)
    run = dict(output_root=str(root))
    result = dict(status="factorial_attempt_failed_retained")
    verification = dict(status="verified_wrapper_failed", error_type="TimeoutError",
                        error=str(root / "measurement/ready.json"))
    terminal = dict(verified=dict(State="FAILED", ExitCode="1:0"))
    command = dict(cpus=32, timeout_s=85800, interval_s=1.)
    mask = hex(sum(1 << cpu for cpu in range(32, 96)))
    log = (f"srun: error: CPU binding outside of job step allocation, allocated CPUs are: {mask}.\n"
        "srun: error: Task launch for StepId=42.0 failed on node bizon: Unable to satisfy cpu bind request\n"
        "srun: error: Application launch failed: Unable to satisfy cpu bind request\n"
        "srun: Job step aborted\n")
    rows = [dict(JobIDRaw=identity, State=state, ExitCode=code, ElapsedRaw=elapsed,
        AllocCPUS="64", ReqMem="128G" if identity == "42" else "", NodeList="bizon")
        for identity, state, code, elapsed in [("42", "FAILED", "1:0", "137"),
            ("42.batch", "FAILED", "1:0", "137"), ("42.0", "CANCELLED", "0:64", "0")]]
    return [request, run, result, verification, terminal, command, log,
            {"abort": True}, {"release": True}, rows]


def test_exact_failure_has_no_native_cost(case):
    result = module.failure_scope(*case)
    assert result["native_outcome"] == "not_started"
    assert result["allocation_elapsed_s"] == 137
    assert result["missing_requested_cpu_ids"] == list(range(32))
    assert result["allocated_os_cpu_ids"] == list(range(32, 96))


@pytest.mark.parametrize("change", ["success", "state", "exit", "wrapper", "error_type", "ready_path",
    "result", "cpus", "cpu_bool", "timeout", "interval", "go", "cleanup", "log_job", "log_host",
    "log_error", "extra_log", "zero_mask", "large_mask", "compatible_mask", "small_allocation",
    "extra_step", "step_running", "step_exit", "step_elapsed", "account_cpu", "account_host",
    "account_mem", "account_job", "duplicate", "elapsed", "abort_integer", "cleanup_integer", "interval_bool"])
def test_refuses_other_or_contradictory_failures(case, change):
    data = deepcopy(case)
    if change == "success": data[4]["verified"].update(State="COMPLETED", ExitCode="0:0")
    elif change == "state": data[4]["verified"]["State"] = "RUNNING"
    elif change == "exit": data[4]["verified"]["ExitCode"] = "0:11"
    elif change == "wrapper": data[3]["status"] = "command_failed"
    elif change == "error_type": data[3]["error_type"] = "ValueError"
    elif change == "ready_path": data[3]["error"] = "different ready.json"
    elif change == "result": data[2]["status"] = "measurement_returned_pending_independent_review"
    elif change == "cpus": data[5]["cpus"] = 16
    elif change == "cpu_bool": data[5]["cpus"] = True
    elif change == "timeout": data[5]["timeout_s"] = 99999
    elif change == "interval": data[5]["interval_s"] = 2.
    elif change == "go": data[7] = {"go": True}
    elif change == "cleanup": data[8] = {"release": False}
    elif change == "abort_integer": data[7] = {"abort": 1}
    elif change == "cleanup_integer": data[8] = {"release": 1}
    elif change == "interval_bool": data[5]["interval_s"] = True
    elif change == "log_job": data[6] = data[6].replace("42.0", "43.0")
    elif change == "log_host": data[6] = data[6].replace("bizon", "other")
    elif change == "log_error": data[6] = data[6].replace("CPU binding outside", "Out of memory")
    elif change == "extra_log": data[6] += "another refusal\n"
    elif change in {"zero_mask", "large_mask", "compatible_mask", "small_allocation"}:
        masks = dict(zero_mask=0, large_mask=(1 << 256) - 1,
            compatible_mask=(1 << 64) - 1, small_allocation=((1 << 48) - 1) << 32)
        data[6] = data[6].replace(hex(sum(1 << cpu for cpu in range(32, 96))), hex(masks[change]))
    elif change == "extra_step": data[9].append(dict(data[9][-1], JobIDRaw="42.1"))
    elif change == "step_running": data[9][-1]["State"] = "RUNNING"
    elif change == "step_exit": data[9][-1]["ExitCode"] = "0:0"
    elif change == "step_elapsed": data[9][-1]["ElapsedRaw"] = "1"
    elif change == "account_cpu": data[9][-1]["AllocCPUS"] = "32"
    elif change == "account_host": data[9][-1]["NodeList"] = "other"
    elif change == "account_mem": data[9][0]["ReqMem"] = "256G"
    elif change == "account_job": data[9][0]["JobIDRaw"] = "43"
    elif change == "duplicate": data[9][-1] = dict(data[9][0])
    else: data[9][0]["ElapsedRaw"] = "-1"
    with pytest.raises(ValueError):
        module.failure_scope(*data)


def test_accounting_parser(case):
    text = "\n".join("|".join(row[k] for k in module.ACCOUNT_FIELDS) for row in case[-1]) + "\n"
    assert module.accounting_rows(text) == case[-1]


@pytest.mark.parametrize("raw", ["", "42|FAILED|1:0\n", "x|x|x|x|x|x|x\n" * 3])
def test_invalid_accounting_parser(raw):
    with pytest.raises(ValueError):
        module.accounting_rows(raw)


def measurement(tmp_path):
    directory = tmp_path / "measurement"
    directory.mkdir()
    for name in module.MEASUREMENT_FILES:
        (directory / name).write_text("" if name.endswith(".jsonl") else json.dumps({}))
    return directory


def test_absent_worker_and_exact_inventory(tmp_path):
    measurement(tmp_path)
    assert len(module.absent_artifacts(tmp_path)) == len(module.ABSENT_MEASUREMENT) + len(module.ABSENT_NATIVE)


@pytest.mark.parametrize("name", [*module.ABSENT_MEASUREMENT, *module.ABSENT_NATIVE, "point_000000.json",
                                   "host_processes.jsonl", "symlink"])
def test_rejects_started_or_unexpected_artifacts(tmp_path, name):
    directory = measurement(tmp_path)
    if name in module.ABSENT_NATIVE:
        (tmp_path / name).write_text("native evidence")
    elif name == "symlink":
        (directory / "step.log").unlink()
        (directory / "step.log").symlink_to(directory / "command.json")
    else:
        (directory / name).write_text("started evidence")
    with pytest.raises(ValueError):
        module.absent_artifacts(tmp_path)


@pytest.fixture
def full_case(tmp_path, monkeypatch, case):
    def store(path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value))
        return module.record(path)

    monkeypatch.setattr(module, "ROOT", tmp_path)
    request, _, result, verification, terminal, command, log, abort, cleanup, rows = case
    root = tmp_path / "panel/run_09"
    directory = root / "measurement"
    directory.mkdir(parents=True)
    session = tmp_path / "panel/sessions/run_09"
    baseline = store(tmp_path / "baseline.json",
        dict(tool_entrypoints=dict(orthohmm_python=dict(absolute_path="/private/python"))))
    plan = dict(panel_root=str(tmp_path / "panel"), baseline=baseline, helper_sources=[], evidence=[])
    plan_ref = store(tmp_path / "plan.json", plan)
    policy_ref = store(tmp_path / "policy.json", {})
    request.update(index=9, plan=plan_ref, policy=policy_ref, evidence=[])
    request_ref = store(tmp_path / "request.json", request)
    run = dict(index=9, dataset="qfo_corrected", cell="p0_c1_r1", repeat=0, output_root=str(root))
    verification.update(error=str(directory / "ready.json"), scientific_results_admitted=False,
        source_sha256=module.record(Path(module.bind_session.__globals__["ROOT"]) /
                                   "benchmark_tools/run_verified_slurm_measurement.py")["sha256"])
    result.update(job_id=42, index=9, cell=run["cell"], plan=plan_ref, request=request_ref,
        execution_scope=module.SCOPE, automatic_retry=False, next_identity_authorized=False,
        uncontended_timing=False, wrapper=verification)
    store(session / "result.json", result)
    store(root / "verification.json", verification)
    command["command"] = module.native_command(plan_ref, run, module.read(baseline))
    store(directory / "command.json", command)
    store(directory / "go.json", abort)
    store(directory / "release.json", cleanup)
    (directory / "host_processes.jsonl").write_text("")
    (directory / "step.log").write_text(log)
    terminal.update(source="fresh_accounting_after_controller_expiry")
    monkeypatch.setattr(module, "validate_plan", lambda data: [None] * 9 + [run])
    monkeypatch.setattr(module, "validate_request", lambda *args: None)
    monkeypatch.setattr(module, "verify_terminal", lambda job: terminal)
    text = "\n".join("|".join(row[k] for k in module.ACCOUNT_FIELDS) for row in rows) + "\n"
    monkeypatch.setattr(module.subprocess, "run",
        lambda argv, **kwargs: subprocess.CompletedProcess(argv, 0, text, ""))
    runtime_calls = []
    def runtime(*args):
        runtime_calls.append(args)
        return dict(status="fixture_runtime_checked", continuous_runtime_integrity_established=False)
    monkeypatch.setattr(module, "runtime_review", runtime)
    return dict(request=request_ref, destination=tmp_path / "review", root=root, runtime_calls=runtime_calls,
                result=result, verification=verification, store=store, session=session, run=run, plan=plan_ref)


def test_failure_disposition_preserves_nonadmission(full_case):
    ref = module.review(full_case["request"], full_case["destination"])
    result = module.read(ref)
    assert result["schema"] == "native_factorial_launch_failure_review_v1"
    assert result["status"] == "pre_native_cpu_binding_failure_reviewed_retained"
    assert result["scheduler_state"] == "FAILED" and result["scheduler_exit_code"] == "1:0"
    assert result["terminal_reviewed"] is True and result["next_identity_authorized"] is True
    assert result["resources"] is None and result["native_outcome"] == "not_started"
    for key in ("native_inference_started", "native_outputs_validated", "accuracy_evaluated",
        "native_command_success", "scheduler_success", "primary_resources_replayed",
        "shared_host_resources_reviewed", "scientific_timings_admitted", "eligible_for_timing_comparison",
        "original_receipts_rewritten", "inference_reexecuted", "automatic_retry", "uncontended_timing",
        "publication_ready"):
        assert result[key] is False
    assert len(full_case["runtime_calls"]) == 1


@pytest.mark.parametrize("change", ["request", "native", "runtime", "binding", "stale_plan", "destination"])
def test_invalid_full_review_never_authorizes(full_case, monkeypatch, change):
    case = full_case
    if change == "request":
        Path(case["request"]["path"]).write_text("changed")
    elif change == "native":
        (case["root"] / "native_execution.json").write_text("started")
    elif change == "runtime":
        def fail(*args): raise ValueError("runtime identity differs")
        monkeypatch.setattr(module, "runtime_review", fail)
    elif change == "binding":
        result = dict(case["result"], index=8)
        case["store"](case["session"] / "result.json", result)
    elif change == "stale_plan":
        def mutate(*args):
            Path(case["plan"]["path"]).write_text("changed")
            return {}
        monkeypatch.setattr(module, "runtime_review", mutate)
    else:
        case["destination"].mkdir()
    with pytest.raises(ValueError):
        module.review(case["request"], case["destination"])
    assert not (case["destination"] / "review.json").exists()
