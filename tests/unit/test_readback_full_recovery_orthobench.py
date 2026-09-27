from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools.readback_full_recovery_orthobench import verify_execution, admit
from benchmark_tools.run_publication_pipeline import arguments


def receipts():
    command = ["/venv/bin/python", "-I", "/launcher.py", "--input", "/input",
               "--output", "/run/native", "--cpu", "32", "--aligner", "/mafft", "--tree-builder", "/tree"]
    environment = dict(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
                       PYTHONHASHSEED="0", MAFFT_BINARIES="/helpers", PATH="/usr/bin:/bin")
    plan = dict(directory="/run", command=command, environment=environment)
    plan_record = dict(path="/run/plan.json", sha256="plan", bytes=1)
    native_started_record = dict(path="/run/native/started.json", sha256="start", bytes=1)
    native_complete_record = dict(path="/run/native/complete.json", sha256="complete", bytes=1)
    submission = dict(job_id="22326", plan=plan_record)
    execution_started = dict(job_id="22326", plan=plan_record)
    execution = dict(status="native_complete_pending_independent_readback", returncode=0,
                     job_id="22326", plan=plan_record, native_complete=native_complete_record,
                     command=["/usr/bin/time", "-v", "-o", "/run/time.txt", *command])
    native_started = dict(command=command[2:], executable=command[0], attempts=1,
        checkpoint_reuse=False, production_default_changed=False, environment=dict(environment),
        arguments=arguments(SimpleNamespace(input=Path("/input"), output=Path("/run/native"),
            cpu=32, aligner=Path("/mafft"), tree_builder=Path("/tree"))))
    native_complete = dict(status="native_complete_pending_scientific_readback",
                           started=native_started_record, accuracy_evaluated=False)
    return [plan, plan_record, submission, execution_started, execution, native_started,
            native_complete, native_started_record, native_complete_record, 22326]


def test_exact_execution_chain():
    verify_execution(*receipts())


@pytest.mark.parametrize("index,key,value", [
    (2, "job_id", "other"), (3, "job_id", "other"), (4, "status", "running"),
    (4, "returncode", 1), (4, "command", ["other"]), (4, "native_complete", {}),
    (5, "executable", "/wrong/python"), (5, "command", []), (5, "attempts", 2),
    (5, "checkpoint_reuse", True), (5, "production_default_changed", True),
    (5, "arguments", []), (5, "environment", {}),
    (6, "accuracy_evaluated", True), (6, "started", {}), (6, "status", "running")])
def test_reject_mismatched_execution_chain(index, key, value):
    rows = receipts()
    rows[index][key] = value
    with pytest.raises(ValueError):
        verify_execution(*rows)


def test_running_job_cannot_read_artifacts(monkeypatch, tmp_path):
    monkeypatch.setattr("benchmark_tools.readback_full_recovery_orthobench.subprocess.check_output",
                        lambda *a, **k: "JobIDRaw|State|ExitCode|Elapsed|AllocCPUS\n22326|RUNNING|0:0|00:01|32\n")
    with pytest.raises(ValueError, match="COMPLETED"):
        admit(tmp_path, tmp_path, 22326)


def test_other_job_rejected_before_scheduler(tmp_path):
    with pytest.raises(ValueError, match="prespecified"):
        admit(tmp_path, tmp_path, 1)
