from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

import benchmark_tools.measure_native_hierarchy_step as module
from benchmark_tools.measure_native_interval_step import evaluate as original_evaluate


@pytest.fixture
def evidence():
    root = Path(__file__).resolve().parents[2]
    report = json.loads((root / "benchmark_tools/results/dgx_cpu_hierarchy_controls_21816.json").read_text())
    points = report["trials"][0]["points"]
    before, after = deepcopy(points[0]["host"][1]), deepcopy(points[-1]["host"][0])
    for sample, point, shift in ((before, points[0], 1000), (after, points[-1], -1000)):
        sample["started_monotonic_ns"] += shift
        sample["finished_monotonic_ns"] += shift
        native = module.interval_point(point, 21816)
        sample["raw"]["cgroup_membership"] = point["native_membership"]
        sample["optional"]["cgroup_cpu.stat"] = native["native_cpu_stat"]
    done = dict(exit_code=0, timed_out=False, snapshots=[before, after],
                started_ns=before["finished_monotonic_ns"]+1000,
                finished_ns=after["started_monotonic_ns"]-1000)
    return points, done


def test_same_native_read_and_original_screen_preserved(evidence):
    points, done = evidence
    converted = [module.interval_point(p, 21816) for p in points]
    for p, c in zip(points, converted):
        native = next(r for r in p["children"] if r["scope"] == c["native_cpu_scope"])
        assert c["host"] is p["host"]
        assert c["native_cpu_stat"] == native["raw"]
        assert c["native_read_ns"] == [native["started_ns"], native["finished_ns"]]
    result = module.evaluate(points, done, 21816)
    assert result["original_threshold_screen"] == original_evaluate(converted, done, 21816)
    assert len(result["hierarchy_intervals"]) == 1
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("fault", ["start", "finish", "missing", "scope", "gap"])
def test_bad_complete_command_evidence_rejected(evidence, fault):
    points, done = evidence
    if fault == "start":
        done["started_ns"] = points[0]["host"][0]["started_monotonic_ns"]
    elif fault == "finish":
        done["finished_ns"] = points[-1]["host"][1]["finished_monotonic_ns"]
    elif fault == "missing":
        points.pop()
    elif fault == "scope":
        points[0]["children"][0]["scope"] += "/task_0"
    else:
        for row in points[-1]["host"]:
            for key in ("started_monotonic_ns", "finished_monotonic_ns"):
                row[key] += 2_000_000_000
        for row in [*points[-1]["parent"], *points[-1]["children"]]:
            row["started_ns"] += 2_000_000_000
            row["finished_ns"] += 2_000_000_000
    with pytest.raises(ValueError):
        module.evaluate(points, done, 21816)


def setup_measure(tmp_path, monkeypatch, evidence, exit_code=0):
    points, done = evidence
    done["exit_code"] = exit_code
    done["timed_out"] = exit_code == 124
    directory = tmp_path / "measurement"
    for key, value in dict(SLURM_JOB_ID="21816", SLURM_CPUS_PER_TASK="20", SLURM_MEM_PER_NODE="98304").items():
        monkeypatch.setenv(key, value)
    monkeypatch.setattr(module.os, "uname", lambda: SimpleNamespace(nodename="spark-7ff0"))
    class Process:
        def __init__(self, command, **kwargs):
            self.finished = False
            assert command[:2] == ["srun", "--exclusive"]
        def poll(self):
            return 0 if self.finished else None
        def wait(self, timeout):
            assert (directory / "release.json").exists()
            self.finished = True
            return 0
    monkeypatch.setattr(module.subprocess, "Popen", Process)
    monkeypatch.setattr(module, "wait_file", lambda p: dict(pid=123, cgroup=points[0]["native_membership"]))
    pending = iter(points)
    monkeypatch.setattr(module, "read_hierarchy", lambda *args: next(pending))
    def sleep(seconds):
        assert (directory / "go.json").exists()
        if not (directory / "done.json").exists():
            (directory / "done.json").write_text(json.dumps(done))
    monkeypatch.setattr(module.time, "sleep", sleep)
    monkeypatch.setattr(module, "step_memory", lambda p: {"scope": p["native_cpu_scope"], "errors": []})
    return directory


@pytest.mark.parametrize("exit_code", [0, 7, 124])
def test_collector_preserves_success_failure_and_timeout(tmp_path, monkeypatch, evidence, exit_code):
    directory = setup_measure(tmp_path, monkeypatch, evidence, exit_code)
    result = module.measure(["/usr/bin/true"], directory, 21816, 20, 96*1024**3, 60, 1.)
    assert result["status"] == ("command_exited_zero" if exit_code == 0 else "command_failed")
    assert result["native"]["exit_code"] == exit_code
    assert result["native"]["timed_out"] == (exit_code == 124)
    assert result["scientific_timings_admitted"] is False
    assert json.loads((directory / "hierarchy_report.json").read_text()) == result
    assert (directory / "point_0001.json").exists()


def test_observation_error_releases_owned_worker(tmp_path, monkeypatch, evidence):
    directory = setup_measure(tmp_path, monkeypatch, evidence)
    def fail(*args):
        raise ValueError("read failed")
    monkeypatch.setattr(module, "read_hierarchy", fail)
    with pytest.raises(ValueError, match="read failed"):
        module.measure(["/usr/bin/true"], directory, 21816, 20, 96*1024**3, 60, 1.)
    assert (directory / "release.json").exists()
    assert (directory / "go.json").exists()
    assert not (directory / "hierarchy_report.json").exists()
