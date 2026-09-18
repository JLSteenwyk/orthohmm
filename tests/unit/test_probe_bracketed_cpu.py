import hashlib
import json
from pathlib import Path

import pytest

from benchmark_tools import probe_bracketed_cpu as module


def test_unscheduled_rejected(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(KeyError):
        module.run(tmp_path / "out")
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("seconds", [0, -1, 4])
def test_unfrozen_cpu_duration_rejected(seconds):
    with pytest.raises(ValueError):
        module.burn(seconds)


def test_evaluation_checks_durations_before_screen():
    with pytest.raises(ValueError, match="native"):
        module.evaluate({}, {}, {"load": {"cpu_s": 1.9}}, None, 123, 100)


def retained():
    root = Path(__file__).resolve().parents[2]
    return root, json.loads((root / "benchmark_tools/results/dgx_bracketed_controls_21805.json").read_text())


def test_actual_controls_replay():
    root, report = retained()
    assert report["job_id"] == 21805
    assert report["all_control_expectations_met"] is True
    assert [trial["mode"] for trial in report["trials"]] == [mode for mode, _ in module.MODES]
    for name, digest in report["sources"].items():
        assert hashlib.sha256((root / "benchmark_tools" / name).read_bytes()).hexdigest() == digest
    for trial in report["trials"]:
        result = module.evaluate(*trial["host_snapshots"], trial["native"], trial["sibling"],
                                 report["job_id"], trial["clock_ticks_per_second"])
        assert result == trial["screen"]
        assert result["reasons"] == ([] if trial["mode"] == "quiet" else ["excess_unassigned_cpu"])
        assert result["controlled_workload_verified"] is False
        assert trial["expectation_met"] is True
    assert report["publication_ready"] is False
    scheduler = (root / "benchmark_tools/results/dgx_bracketed_scheduler_21805.txt").read_text()
    for field in ("JobId=21805 ", "JobState=COMPLETED ", "ExitCode=0:0", "Restarts=0", "NodeList=spark-7ff0", "CPUs/Task=2", "MinMemoryNode=2G"):
        assert field in scheduler


@pytest.mark.parametrize("change", ["duration", "scope", "bracket"])
def test_control_evidence_tampering(change):
    _, report = retained()
    trial = report["trials"][1]
    if change == "duration":
        trial["sibling"]["cpu_s"] = 0
    elif change == "scope":
        trial["sibling"]["cgroup"] = "0::/job_21805/step_0\n"
    else:
        trial["host_snapshots"][0]["finished_monotonic_ns"] = trial["native"]["snapshots"][0]["finished_monotonic_ns"] + 1
    with pytest.raises(ValueError):
        module.evaluate(*trial["host_snapshots"], trial["native"], trial["sibling"],
                        report["job_id"], trial["clock_ticks_per_second"])
