import json
import subprocess

import pytest

from benchmark_tools.replay_array_scheduler_capture import replay
from tests.unit.test_capture_array_scheduler import record, run_capture


@pytest.mark.parametrize("prefix", [[], [subprocess.TimeoutExpired("scontrol", 10)],
    [record().replace(" Command=/tmp/run.sh", "")]])
def test_replay_preserves_first_terminal_and_failed_tasks(tmp_path, prefix):
    changed = record().replace("NumCPUs=20", "NumCPUs=40")
    run_capture(tmp_path, prefix + [record(), changed + record(1, "FAILED")])
    result = replay(tmp_path / "capture", 100, 2)
    assert result["status"] == "first_terminal_scheduler_capture_replayed"
    assert result["retained"][0]["poll"] == len(prefix)
    assert result["retained"][1]["poll"] == len(prefix) + 1
    assert result["observation_errors"] == len(prefix)
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("problem", ["missing", "extra", "command", "retained", "summary", "error", "symlink"])
def test_changed_evidence_rejected(tmp_path, problem):
    run_capture(tmp_path, [record(), record(1)])
    directory = tmp_path / "capture"
    poll = directory / "poll_000000.json"
    if problem == "missing":
        poll.unlink()
    elif problem == "extra":
        (directory / "poll_000002.json").write_text(poll.read_text())
    elif problem == "command":
        value = json.loads(poll.read_text())
        value["argv"] = ["sacct"]
        poll.write_text(json.dumps(value))
    elif problem == "retained":
        (directory / "scheduler_0.txt").write_text(record().replace("NumCPUs=20", "NumCPUs=40"))
    elif problem == "symlink":
        destination = tmp_path / "moved.json"
        poll.rename(destination)
        poll.symlink_to(destination)
    else:
        path = directory / "capture.json"
        value = json.loads(path.read_text())
        if problem == "summary":
            value["retained"]["0"]["poll"] = 1
        else:
            value["observation_errors"] = 1
        path.write_text(json.dumps(value))
    with pytest.raises(ValueError):
        replay(directory, 100, 2)


def test_incomplete_capture_rejected(tmp_path):
    run_capture(tmp_path, [record(1)])
    with pytest.raises(ValueError, match="complete"):
        replay(tmp_path / "capture", 100, 2)
