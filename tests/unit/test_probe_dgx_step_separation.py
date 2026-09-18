import json
from pathlib import Path

import pytest

from benchmark_tools import probe_dgx_step_separation as module


def sample(scope):
    return {"raw": {"cgroup_membership": "0::" + scope + "\n"}}


def test_separate_step_scopes():
    parent = sample("/system.slice/slurm/job_123/step_batch/user/task_0")
    native = sample("/system.slice/slurm/job_123/step_0/user/task_0")
    assert module.validate_scopes(parent, native, 123)["native_step"].endswith("/step_0")


@pytest.mark.parametrize("scope", ["/job_124/step_0", "/job_123/step_batch", "/job_123/step_extern",
                                    "/job_123/user", "/job_123", "/../job_123/step_0"])
def test_invalid_scopes(scope):
    with pytest.raises(ValueError):
        module.validate_scopes(sample("/job_123/step_batch"), sample(scope), 123)


def test_atomic_publication_and_no_overwrite(tmp_path):
    path = tmp_path / "ready.json"
    module.save(path, {"ready": True})
    assert module.wait_file(path) == {"ready": True}
    assert not path.with_name("ready.json.pending").exists()
    with pytest.raises(FileExistsError):
        module.save(path, {"changed": True})
    assert json.loads(path.read_text()) == {"ready": True}


def test_wait_timeout(tmp_path):
    with pytest.raises(TimeoutError):
        module.wait_file(tmp_path / "missing", seconds=.01)


def test_unscheduled_rejected(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(KeyError):
        module.run(tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_retained_dgx_trials_replay_and_bracket():
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/dgx_step_separation_probe_21800.json"
    report = json.loads(path.read_text())
    assert report["job_id"] == 21800
    assert [t["mode"] for t in report["trials"]] == ["quiet", "burst"]
    for trial in report["trials"]:
        before, after = trial["snapshots"]
        native_before, native_after = trial["native_snapshots"]
        assert native_before["finished_monotonic_ns"] < before["started_monotonic_ns"]
        assert after["finished_monotonic_ns"] < native_after["started_monotonic_ns"]
        assert module.validate_scopes(before, native_before, 21800) == trial["scopes"]
        assert module.validate_scopes(after, native_after, 21800) == trial["scopes"]
        assert all(not s["errors"] for s in (before, after, native_before, native_after))
        assert module.summarize(before, after, 100) == trial["summary"]
        assert trial["native_exit_code"] == 0
    burst = report["trials"][1]
    assert .75 <= burst["burst"]["cpu_seconds"] < 5
    assert burst["summary"]["accounted_host_busy_cpu_s"] >= .5
    assert report["controlled_workload_verified"] is False
