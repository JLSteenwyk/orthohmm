import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import measure_scaling_root_context as module
from benchmark_tools.measure_native_root_context import measure_native_run as historical


def test_long_settings_accepted_but_historical_boundary_unchanged():
    module.validate(["/bin/true"], 20, 85800, 1.)
    with pytest.raises(ValueError, match="frozen"):
        historical(["/bin/true"], Path("unused"), 1, 20, 96*1024**3, 85800, 1.)


@pytest.mark.parametrize("command,cpus,timeout,interval", [
    ([],20,85800,1.), (["relative"],20,85800,1.), (["/bin/true", ""],20,85800,1.),
    (["/bin/true"],True,85800,1.), (["/bin/true"],20,900,1.),
    (["/bin/true"],20,85800.0,1.), (["/bin/true"],20,85800,True),
    (["/bin/true"],20,85800,0.5)])
def test_other_settings_rejected(command, cpus, timeout, interval):
    with pytest.raises(ValueError): module.validate(command, cpus, timeout, interval)


@pytest.mark.parametrize("gate", [{"abort": True}, {"go": 1}, {"go": False}, {}])
def test_worker_abort_never_launches_native(tmp_path, monkeypatch, gate):
    module.save(tmp_path / "command.json", dict(command=["/bin/true"], cpus=20, timeout_s=85800, interval_s=1.))
    monkeypatch.setattr(module, "wait_file", lambda path: gate)
    monkeypatch.setattr(module, "run_command", lambda *a: pytest.fail("unobserved native launch"))
    module.worker(tmp_path)
    assert (tmp_path / "aborted_before_native.json").exists()
    assert not (tmp_path / "done.json").exists()


@pytest.mark.parametrize("code,timed_out", [(0,False), (7,False), (124,True)])
def test_worker_retains_native_status_and_exact_long_timeout(tmp_path, monkeypatch, code, timed_out):
    command = ["/bin/true"]
    module.save(tmp_path / "command.json", dict(command=command, cpus=20, timeout_s=85800, interval_s=1.))
    def wait(path):
        if path.name == "release.json":
            done = json.loads((tmp_path / "done.json").read_text())
            assert done["exit_code"] == code and done["timed_out"] is timed_out
            return {"release": True}
        return {"go": True}
    monkeypatch.setattr(module, "wait_file", wait)
    monkeypatch.setattr(module, "snapshot", lambda: {"raw": "preserved"})
    def run(actual, log, timeout):
        assert actual == command and timeout == 85800
        return code, timed_out
    monkeypatch.setattr(module, "run_command", run)
    module.worker(tmp_path)
    assert (tmp_path / "native.log").exists()


def allocation(monkeypatch):
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "20")
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "98304")
    monkeypatch.setattr(module.os, "uname", lambda: SimpleNamespace(nodename="spark-7ff0"))


@pytest.mark.parametrize("code,timed_out", [(0,False), (7,False), (124,True)])
@pytest.mark.parametrize("cycles", [1, 65])
def test_observer_retains_reports_final_reads_and_failure_status(tmp_path, monkeypatch, code, timed_out, cycles):
    allocation(monkeypatch)
    directory = tmp_path / "measurement"
    waited = []
    class Process:
        finished = False
        def poll(self): return 0 if self.finished else None
        def wait(self, timeout):
            waited.append(timeout)
            assert (directory / "release.json").exists()
            self.finished = True
            return 0
    process = Process()
    def launch(argv, **kwargs):
        assert argv[:7] == ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=20", module.sys.executable]
        assert argv[-3:] == [str(Path(module.__file__).resolve()), "--worker", str(directory)]
        return process
    monkeypatch.setattr(module.subprocess, "Popen", launch)
    monkeypatch.setattr(module, "wait_file", lambda path: {"pid": 55, "cgroup": "0::/slurm/job_123/step_0/user/task_0\n"})
    observations = []
    clock = [0.]
    class Monitor:
        def __init__(self, handle, scope):
            assert scope == "/slurm/job_123"
            self.handle = handle
        def observe(self):
            if not observations:
                assert not (directory / "go.json").exists()
            observations.append(clock[0])
            self.handle.write(json.dumps({"time": clock[0]}) + "\n")
        def summary(self, start, end):
            assert start == done["started_ns"] / 1e9
            assert end == done["finished_ns"] / 1e9
            assert (directory / "done.json").exists()
            assert observations[-1] == cycles
            return {"controlled_workload_verified": False, "snapshots": len(observations)}
    monkeypatch.setattr(module, "HostMonitor", Monitor)
    monkeypatch.setattr(module.time, "monotonic", lambda: clock[0])
    points = []
    def read(*args):
        value = {"index": len(points)}
        points.append(value)
        return value
    monkeypatch.setattr(module, "read_point", read)
    done = dict(exit_code=code, timed_out=timed_out, started_ns=100, finished_ns=1000000100)
    def sleep(duration):
        clock[0] += duration
        if clock[0] >= cycles:
            module.save(directory / "done.json", done)
    monkeypatch.setattr(module.time, "sleep", sleep)
    monkeypatch.setattr(module, "interval_point", lambda point, job: {"final": point["index"]})
    monkeypatch.setattr(module, "step_memory", lambda point: {"observed_after_final_point": point})
    monkeypatch.setattr(module, "evaluate_lineage", lambda p,d,j: {"points": len(p), "done": d})
    monkeypatch.setattr(module, "evaluate", lambda p,j: {"points": len(p)})
    result = module.measure(["/bin/true"], directory, 123, 20, 96*1024**3, 85800, 1.)
    assert result["status"] == ("command_exited_zero" if code == 0 else "command_failed")
    assert result["native"] == done and result["native_wall_s"] == 1.
    assert result["points"] == points and len(points) == cycles + 1
    assert waited == [45]
    assert [p.name for p in sorted(directory.glob("point_*.json"))] == [f"point_{i:06d}.json" for i in range(cycles + 1)]
    assert observations == ([0., 1.] if cycles == 1 else [0., 30., 60., 65.])
    summary = json.loads((directory / "host_process_summary.json").read_text())
    assert result["host_process_observation"] == summary
    assert summary["controlled_workload_verified"] is False
    assert len((directory / "host_processes.jsonl").read_text().splitlines()) == len(observations)
    context = json.loads((directory / "root_context_report.json").read_text())
    assert context["lineage_report"] == module.lineage_identity(directory)
    assert result["scientific_timings_admitted"] is False


def test_initial_observer_failure_aborts_worker_without_go(tmp_path, monkeypatch):
    allocation(monkeypatch)
    directory = tmp_path / "measurement"
    class Process:
        def poll(self): return None
        def wait(self, timeout):
            assert json.loads((directory / "go.json").read_text()) == {"abort": True}
            assert (directory / "release.json").exists()
            assert timeout == 85890
            return 0
    monkeypatch.setattr(module.subprocess, "Popen", lambda *a, **k: Process())
    monkeypatch.setattr(module, "wait_file", lambda p: {"pid": 55, "cgroup": "scope"})
    def fail(*a): raise ValueError("failed initial read")
    monkeypatch.setattr(module, "read_point", fail)
    with pytest.raises(ValueError, match="initial"):
        module.measure(["/bin/true"], directory, 123, 20, 96*1024**3, 85800, 1.)


def test_six_digit_point_names_sort_beyond_ten_thousand():
    names = [f"point_{i:06d}.json" for i in (0, 9999, 10000, 85831)]
    assert sorted(names) == names


@pytest.mark.parametrize("code", [0, 7])
def test_worker_runs_real_short_subprocess_with_long_timeout_contract(tmp_path, monkeypatch, code):
    command = [module.sys.executable, "-c", f"print('scaling worker smoke'); raise SystemExit({code})"]
    module.save(tmp_path / "command.json", dict(command=command, cpus=20, timeout_s=85800, interval_s=1.))
    monkeypatch.setattr(module, "wait_file", lambda path: {"go": True} if path.name == "go.json" else {"release": True})
    module.worker(tmp_path)
    done = json.loads((tmp_path / "done.json").read_text())
    assert done["exit_code"] == code and done["timed_out"] is False
    assert done["started_ns"] < done["finished_ns"]
    assert len(done["snapshots"]) == 2
    assert "scaling worker smoke" in (tmp_path / "native.log").read_text()
