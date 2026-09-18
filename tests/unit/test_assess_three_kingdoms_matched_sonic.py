from copy import deepcopy

import pytest

from benchmark_tools import assess_three_kingdoms_matched_sonic as module


@pytest.fixture
def case():
    plan = {"output_root": "/run", "native_argv": ["sonic", "-t", "32"],
            "inputs": [{"path": "/input/a.fasta", "bytes": 10, "sha256": "abc"}]}
    execution = {"status": "process_succeeded_pending_native_admission", "exit_code": 0,
                 "accuracy_admitted": False, "job_id": "21795", "node": "bizon",
                 "plan": {"sha256": "plan"}, "source": {"sha256": "runner"},
                 "native_argv": plan["native_argv"], "started_epoch": 1, "finished_epoch": 2,
                 "runtime_before": {"tree": "ok"}, "runtime_after": {"tree": "ok"},
                 "copied_inputs": [{"path": "/run/input/a.fasta", "bytes": 10, "sha256": "abc"}],
                 "native_log": {"path": "/run/native.log", "bytes": 1},
                 "timing": {"path": "/run/time.txt", "bytes": 1}}
    scheduler = {"JobIDRaw": "21795", "State": "COMPLETED", "ExitCode": "0:0",
                 "NodeList": "bizon", "AllocCPUS": "32", "ReqMem": "192G"}
    return plan, execution, scheduler, deepcopy(execution["plan"]), deepcopy(execution["source"])


def test_valid_execution(case):
    module.validate_execution(*case)


@pytest.mark.parametrize("key,value", [("State", "RUNNING"), ("ExitCode", "1:0"),
    ("JobIDRaw", "21796"), ("AllocCPUS", "16"), ("ReqMem", "64G"), ("NodeList", "other")])
def test_wrong_scheduler_rejected(case, key, value):
    case[2][key] = value
    with pytest.raises(ValueError):
        module.validate_execution(*case)


@pytest.mark.parametrize("key,value", [("status", "running"), ("exit_code", 1),
    ("accuracy_admitted", True), ("job_id", "other"), ("node", "other"),
    ("plan", {}), ("source", {}), ("native_argv", []), ("started_epoch", 3),
    ("runtime_after", {}), ("copied_inputs", []), ("timing", {"path": "/bad", "bytes": 1})])
def test_wrong_execution_rejected(case, key, value):
    case[1][key] = value
    with pytest.raises(ValueError):
        module.validate_execution(*case)


def test_pending_job_cannot_create_assessment(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
                        "JobIDRaw|State|ExitCode\n21795|PENDING|0:0\n")
    destination = tmp_path / "new"
    with pytest.raises(ValueError, match="COMPLETED"):
        module.assess(tmp_path, destination)
    assert not destination.exists()


def test_existing_destination_preserved(tmp_path):
    with pytest.raises(FileExistsError):
        module.assess(tmp_path, tmp_path)
