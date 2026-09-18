import pytest

from benchmark_tools.admit_qfo_corrected_sonic import validate_execution


def fixture():
    plan = {"native_argv": ["sonic"], "cwd": "/run", "output_root": "/run",
            "copy_inputs_to": "/run/input", "input_fastas": [{"path": "/data/a.fasta", "sha256": "a", "bytes": 1}]}
    source = {"path": "/runner", "sha256": "r", "bytes": 2}
    manifest = {"path": "/plan", "sha256": "p", "bytes": 3}
    scheduler = {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "32", "JobIDRaw": "1"}
    execution = {"status": "process_succeeded_pending_native_admission", "exit_code": 0,
                 "source": dict(source), "plan": dict(manifest), "job_id": "1", "node": "bizon",
                 "native_argv": ["sonic"], "cwd": "/run", "started_epoch": 1, "finished_epoch": 2,
                 "runtime_before": [{"verified": True}], "runtime_after": [{"verified": True}],
                 "copied_inputs": [{"path": "/run/input/a.fasta", "sha256": "a", "bytes": 1}],
                 "log": {"path": "/run/native.log", "bytes": 1}, "timing": {"path": "/run/time.txt", "bytes": 1}}
    return plan, execution, scheduler, manifest, source


def test_valid():
    validate_execution(*fixture())


@pytest.mark.parametrize("key,value", [
    ("status", "running"), ("exit_code", 1), ("job_id", "2"), ("node", "other"),
    ("source", {}), ("plan", {}), ("native_argv", ["other"]), ("cwd", "/old"),
    ("finished_epoch", 0), ("runtime_after", []), ("runtime_before", []),
    ("copied_inputs", []), ("log", {"path": "/other", "bytes": 1}),
])
def test_reject_execution_change(key, value):
    args = fixture()
    args[1][key] = value
    with pytest.raises(ValueError):
        validate_execution(*args)


@pytest.mark.parametrize("key,value", [("State", "RUNNING"), ("ExitCode", "1:0"),
                                      ("AllocCPUS", "8"), ("NodeList", "other")])
def test_reject_scheduler_change(key, value):
    args = fixture()
    args[2][key] = value
    with pytest.raises(ValueError, match="scheduler"):
        validate_execution(*args)
