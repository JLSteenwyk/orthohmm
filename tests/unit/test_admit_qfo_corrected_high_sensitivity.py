from copy import deepcopy
import pytest

from benchmark_tools.admit_qfo_corrected_high_sensitivity import (
    METHOD, admit, validate_execution, validate_metrics_execution,
)


def fixture():
    config = {"native_argv": ["/python", "-m", "orthohmm", "/inputs", "--stop", "infer"], "cwd": "/core"}
    plan = {"output_root": "/output", "methods": {METHOD: config}}
    execution = {"status": "process_succeeded_pending_native_admission", "exit_code": 0,
        "method": METHOD, "index": 0, "array_task_id": "0", "array_job_id": "100",
        "job_id": "101", "node": "bizon", "source": {"sha256": "runner"},
        "manifest": {"sha256": "plan"}, "native_argv": deepcopy(config["native_argv"]),
        "cwd": "/core", "started_epoch": 1., "finished_epoch": 4.,
        "accuracy_admitted": False, "native_outputs_validated": False,
        "log": {"path": f"/output/execution/{METHOD}/native.log", "bytes": 10},
        "timing": {"path": f"/output/execution/{METHOD}/time.txt", "bytes": 10}}
    scheduler = {"JobIDRaw": "101", "State": "COMPLETED", "ExitCode": "0:0",
                 "NodeList": "bizon", "AllocCPUS": "32"}
    return plan, execution, scheduler


def check_fixture(plan, execution, scheduler):
    validate_execution(plan, execution, scheduler, {"sha256": "plan"}, {"sha256": "runner"})


def test_completed_execution():
    check_fixture(*fixture())


@pytest.mark.parametrize("key,value", [
    ("status", "running"), ("exit_code", 1), ("method", "orthofinder_full"),
    ("index", 1), ("index", False), ("array_task_id", "1"), ("array_job_id", ""),
    ("job_id", "102"), ("node", "other"), ("source", {}), ("manifest", {}),
    ("native_argv", []), ("cwd", "/wrong"), ("started_epoch", 5.),
    ("finished_epoch", float("nan")), ("accuracy_admitted", True), ("native_outputs_validated", True),
    ("log", {"path": "/wrong", "bytes": 10}),
    ("timing", {"path": f"/output/execution/{METHOD}/time.txt", "bytes": 0}),
])
def test_reject_execution_change(key, value):
    plan, execution, scheduler = fixture()
    execution[key] = value
    with pytest.raises(ValueError):
        check_fixture(plan, execution, scheduler)


@pytest.mark.parametrize("key,value", [("State", "RUNNING"), ("ExitCode", "1:0"),
    ("NodeList", "other"), ("AllocCPUS", "16"), ("JobIDRaw", "102")])
def test_reject_scheduler_change(key, value):
    plan, execution, scheduler = fixture()
    scheduler[key] = value
    with pytest.raises(ValueError):
        check_fixture(plan, execution, scheduler)


def native_metrics():
    return {"command": ["/python", "/core/orthohmm/__main__.py", "/inputs", "--stop", "infer"],
            "cwd": "/core", "started_at_epoch_s": 2., "finished_at_epoch_s": 3.}


def test_native_module_metrics():
    plan, execution, _ = fixture()
    validate_metrics_execution(native_metrics(), plan["methods"][METHOD], execution)


@pytest.mark.parametrize("key,value", [("command", ["wrong"]), ("cwd", "/wrong"),
    ("started_at_epoch_s", 0.), ("finished_at_epoch_s", 5.), ("finished_at_epoch_s", 1.),
    ("started_at_epoch_s", float("nan"))])
def test_reject_metrics_provenance(key, value):
    plan, execution, _ = fixture()
    metrics = native_metrics()
    metrics[key] = value
    with pytest.raises(ValueError):
        validate_metrics_execution(metrics, plan["methods"][METHOD], execution)


def test_running_job_rejected_before_file_reads(tmp_path, monkeypatch):
    monkeypatch.setattr("benchmark_tools.admit_qfo_corrected_high_sensitivity.subprocess.check_output",
        lambda *args, **kwargs: "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n101|RUNNING|0:0|1:00|bizon|32\n")
    output = tmp_path / "admission.json"
    with pytest.raises(ValueError, match="COMPLETED"):
        admit(tmp_path, 101, output)
    assert not output.exists()
