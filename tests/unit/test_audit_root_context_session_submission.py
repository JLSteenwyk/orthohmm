from datetime import datetime
from zoneinfo import ZoneInfo
import json
from pathlib import Path
import tarfile

import pytest

from benchmark_tools import audit_root_context_session_submission as module


@pytest.fixture
def evidence():
    start = int(datetime(2026, 9, 19, 19, 0, tzinfo=ZoneInfo("America/New_York")).timestamp())*10**9
    queue = dict(command=["squeue", "-h", "-p", "spark", "-o", "%i %T %j"], returncode=0, stdout="",
                 started_unix_ns=start-2*10**9, finished_unix_ns=start-10**9)
    launch = dict(command=module.command("a"*64), recipe_sha256="a"*64, local_timeout_s=1050,
        remote_timeout_s=1020, remote_kill_grace_s=10, scientific_timings_admitted=False, started_unix_ns=start-500000000)
    result = dict(source=dict(sha256="source", bytes=100), status="wait_returned", returncode=0, job_id=123,
        stdout="123\n", recipe_sha256="a"*64, scientific_timings_admitted=False, scheduler_terminal_verified=False,
        finished_unix_ns=start+301*10**9)
    recipe = dict(records=[dict(kind="file", path=str(module.RECIPE / "benchmark_tools/submit_root_context_session.py"),
                                sha256="source", bytes=100)])
    allocation = dict(ExitCode="0:0", StartTime="2026-09-19T19:00:00", EndTime="2026-09-19T19:05:00")
    return queue, launch, result, allocation, recipe


def validate(evidence):
    return module.validate(*evidence, "a"*64, 123, "America/New_York")


def test_wait_encloses_job_and_does_not_admit_timing(evidence):
    assert validate(evidence)["status"] == "bounded_session_submission_verified"
    assert validate(evidence)["scientific_timings_admitted"] is False


def test_retained_session_receipts_match_real_terminal_job():
    results = Path(module.__file__).parent / "results"
    with tarfile.open(results / "root_context_session_receipts_22020.tar.gz") as archive:
        receipts = [json.load(archive.extractfile("root_context_session_submission_v2/" + name))
                    for name in ("queue.json", "launch.json", "result.json")]
        raw = archive.extractfile("root_context_scheduler_22020/scheduler_22020.txt").read().decode()
    allocation = module.scheduler(raw, 22020, "v2")
    recipe = json.loads((results / "root_context_session_recipe_20260919.json").read_text())
    result = module.validate(*receipts, allocation, recipe,
        "a0a2fdbb98fe75b057b2acfb3d9eb94ac1bb6acab903c5a5f71a30ce32a4b627", 22020, "America/New_York")
    assert result["status"] == "bounded_session_submission_verified"


@pytest.mark.parametrize("fault", ["source", "queue", "command", "timeout", "job", "exit", "early", "late", "clock"])
def test_inconsistent_session_receipt_rejected(evidence, fault):
    queue, launch, result, allocation, recipe = evidence
    if fault == "source":
        result["source"]["sha256"] = "other"
    elif fault == "queue":
        queue["stdout"] = "999 RUNNING other"
    elif fault == "command":
        launch["command"] = []
    elif fault == "timeout":
        result["status"] = "observation_timeout"
    elif fault == "job":
        result["job_id"] = 999
    elif fault == "exit":
        result["returncode"] = 1
    elif fault == "early":
        result["finished_unix_ns"] = launch["started_unix_ns"]+100*10**9
    elif fault == "late":
        launch["started_unix_ns"] += 10*10**9
    else:
        queue["started_unix_ns"] = result["finished_unix_ns"]
    with pytest.raises(ValueError):
        validate(evidence)
