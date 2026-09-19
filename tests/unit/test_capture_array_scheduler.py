import json
import subprocess
from types import SimpleNamespace

import pytest

from benchmark_tools.capture_array_scheduler import capture, terminal_records


def record(index=0, state="COMPLETED"):
    return (f"JobId={101 + index} ArrayJobId=100 ArrayTaskId={index} JobState={state} "
            "ExitCode=0:0 Restarts=0 Requeue=0 NodeList=node OverSubscribe=NO "
            "MinMemoryNode=96G NumNodes=1 NumCPUs=20 CPUs/Task=20 Command=/tmp/run.sh\n")


def test_exact_terminal_records():
    raw = record() + record(1, "FAILED") + record(2, "RUNNING")
    assert terminal_records(raw, 100, {0, 1, 2}) == {0: record(), 1: record(1, "FAILED")}


@pytest.mark.parametrize("raw", ["100_0|COMPLETED|0:0|1|20|96G|node\n",
                                    record().replace("ArrayJobId=100", "ArrayJobId=200"),
                                    record().replace("ArrayTaskId=0", "ArrayTaskId=0-17%1")])
def test_no_accounting_or_compressed_substitution(raw):
    assert terminal_records(raw, 100, {0}) == {}


@pytest.mark.parametrize("raw", [record() + record(), record().rstrip() + " JobId=1",
                                    record().replace(" Command=/tmp/run.sh", ""),
                                    record().replace("ExitCode=0:0", "ExitCode=bad")])
def test_bad_terminal_evidence(raw):
    with pytest.raises(ValueError):
        terminal_records(raw, 100, {0})


def run_capture(tmp_path, observations, count=2):
    now = [0]
    queue = iter(observations)

    def run(argv, **kwargs):
        assert argv == ["scontrol", "show", "job", "100", "--oneliner"]
        assert kwargs["timeout"] == 10
        result = next(queue, "")
        if isinstance(result, Exception):
            raise result
        return SimpleNamespace(returncode=0, stdout=result, stderr="")

    def sleep(seconds):
        now[0] += seconds

    return capture(100, count, tmp_path / "capture", run=run, clock=lambda: now[0],
                   sleep=sleep, interval=1, max_seconds=4)


def test_expired_first_task_preserved_while_second_finishes(tmp_path):
    report = run_capture(tmp_path, [record() + record(1, "RUNNING"), record(1)])
    assert report["status"] == "complete_controller_capture"
    assert report["polls"] == 2
    assert (tmp_path / "capture/scheduler_0.txt").read_text() == record()
    assert report["scientific_timings_admitted"] is False


def test_preexisting_expiry_is_incomplete(tmp_path):
    report = run_capture(tmp_path, [record(1)])
    assert report["missing_tasks"] == [0]
    assert report["status"] == "incomplete_controller_capture"


def test_transient_timeout_does_not_mean_job_terminal(tmp_path):
    report = run_capture(tmp_path, [subprocess.TimeoutExpired("scontrol", 10), record() + record(1)])
    assert report["status"] == "complete_controller_capture"
    assert report["observation_errors"] == 1
    assert json.loads((tmp_path / "capture/poll_000000.json").read_text())["error_type"] == "TimeoutExpired"


def test_failed_job_is_retained_not_omitted(tmp_path):
    report = run_capture(tmp_path, [record(0, "FAILED")], count=1)
    assert report["status"] == "complete_controller_capture"
    assert "FAILED" in (tmp_path / "capture/scheduler_0.txt").read_text()


def test_existing_directory_never_overwritten(tmp_path):
    directory = tmp_path / "capture"
    directory.mkdir()
    with pytest.raises(FileExistsError):
        capture(100, 2, directory)


def test_first_terminal_snapshot_not_replaced(tmp_path):
    changed = record().replace("NumCPUs=20", "NumCPUs=40")
    report = run_capture(tmp_path, [record(), changed + record(1)])
    assert report["status"] == "complete_controller_capture"
    assert (tmp_path / "capture/scheduler_0.txt").read_text() == record()
    assert changed.strip() in json.loads((tmp_path / "capture/poll_000001.json").read_text())["stdout"]


def test_parse_failure_keeps_raw_evidence(tmp_path):
    malformed = record().replace(" Command=/tmp/run.sh", "")
    report = run_capture(tmp_path, [malformed], count=1)
    assert report["missing_tasks"] == [0]
    assert report["observation_errors"] == 1
    observed = json.loads((tmp_path / "capture/poll_000000.json").read_text())
    assert observed["stdout"] == malformed
    assert observed["parse_error"] == "Incomplete detailed terminal record"


def test_nonzero_command_does_not_admit_stdout(tmp_path):
    report = capture(100, 1, tmp_path / "capture", max_seconds=.01, interval=.01,
                     run=lambda *a, **k: SimpleNamespace(returncode=1, stdout=record(), stderr="failed"))
    assert report["missing_tasks"] == [0]
    assert report["observation_errors"] >= 1


@pytest.mark.parametrize("interval", [0, -1, float("nan"), float("inf")])
def test_invalid_interval(tmp_path, interval):
    with pytest.raises(ValueError):
        capture(100, 2, tmp_path / "capture", interval=interval)
