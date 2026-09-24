from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools import admit_blast_recovery_bpo as module


@pytest.mark.parametrize("problem", [None, "pending", "cpu", "memory", "status", "job", "source",
                                    "checkpoint", "accuracy", "downstream", "time", "runtime"])
def test_preparation_contract(problem):
    scheduler = dict(JobIDRaw="fixture", State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="2", ReqMem="64G")
    runtime = dict(status="dedicated_bpo_python_runtime_verified", manifest={"fixture": "runtime"}, mapped_files=[])
    report = dict(status="recovered_bpo_prepared_pending_admission", job_id="fixture", source={"fixture": "source"},
        checkpoint={"fixture": "checkpoint"}, started_epoch=1, finished_epoch=2,
        runtime_before=runtime, runtime_after=runtime.copy(), accuracy_admitted=False, publication_ready=False,
        downstream_execution_authorized=False)
    if problem in ("pending", "cpu", "memory"):
        key, value = {"pending": ("State", "PENDING"), "cpu": ("AllocCPUS", "1"), "memory": ("ReqMem", "32G")}[problem]
        scheduler[key] = value
    elif problem:
        key, value = {"status": ("status", "failed"), "job": ("job_id", "other"), "source": ("source", {}),
            "checkpoint": ("checkpoint", {}), "accuracy": ("accuracy_admitted", True),
            "downstream": ("downstream_execution_authorized", True), "time": ("finished_epoch", float("nan")),
            "runtime": ("runtime_after", {**runtime, "manifest": {}})}[problem]
        report[key] = value
    args = (report, scheduler, {"fixture": "source"}, {"fixture": "checkpoint"}, {"fixture": "runtime"})
    if problem:
        with pytest.raises(ValueError):
            module.preparation_contract(*args)
    else:
        module.preparation_contract(*args)


def test_reviewed_source_bytes():
    assert module.record(Path(module.__file__).with_name("prepare_blast_recovery_bpo.py"))["sha256"] == module.PREPARER_SHA
    assert module.record(Path(module.__file__).with_name("admit_qfo_corrected_bpo.py"))["sha256"] == module.VALIDATOR_SHA


def test_cli_help_in_isolation(tmp_path):
    done = subprocess.run([sys.executable, "-I", str(Path(module.__file__).resolve()), "--help"],
                          cwd=tmp_path, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
