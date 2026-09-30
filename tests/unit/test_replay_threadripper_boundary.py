import json
from pathlib import Path
import shutil
import sys

import pytest

from benchmark_tools import measure_threadripper_scaling as periodic
from benchmark_tools import replay_threadripper_boundary as module
from benchmark_tools.command_host_monitor import HostMonitor
from benchmark_tools.native_completion import completion_evidence
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_command_host_monitor import sample
from tests.unit.test_measure_native_hierarchy_step import evidence
from tests.unit.test_measure_native_lineage_step import extend_lineage
from tests.unit.test_native_root_context import add_context


def save(path, value):
    path.write_text(json.dumps(value))


@pytest.fixture
def archive(tmp_path, evidence):
    points, done = extend_lineage(evidence)
    add_context(points)
    scope = module.scoped_path(points[0]["native_membership"], 21816)
    user = scope.parent
    for point in points:
        t = point["root_context"]["host_after"]["finished_ns"] + 10
        point["thread_affinity"] = dict(scope="/sys/fs/cgroup" + str(user), allowed_cpus=list(range(32)),
            full_run_affinity_verified=False, started_ns=t, finished_ns=t+1,
            initial_tids=[123], final_tids=[123], threads=[dict(tid=123, start_ticks=1,
            cgroup=point["native_membership"], affinity=list(range(32)), outside_cpus=[])],
            errors=[], violating_tids=[], status="observed_within_affinity")
    placement = dict(host="bizon", pid=123, affinity=list(range(32)),
        topology=[dict(cpu=i, package=0, core=i) for i in range(32)],
        cgroup=points[0]["native_membership"],
        ancestors=[{"memory.max": str(128 * 1024**3)}],
        slurm=dict(SLURM_JOB_ID="21816", SLURM_CPUS_PER_TASK="64", SLURM_MEM_PER_NODE="131072"))
    end = points[1]["thread_affinity"]["finished_ns"]
    memory = dict(scope=module.interval_point(points[1], 21816)["native_cpu_scope"],
        errors=[], started_ns=end+10, finished_ns=end+20,
        raw={"memory.current": "100", "memory.peak": "200",
             "memory.events": "low 0\nhigh 0\nmax 0\noom 0\noom_kill 0\n"})
    job_scope = next(p for p in scope.parents if p.name == "job_21816")
    def job_memory(start, peak):
        return dict(scope=str(job_scope), errors=[], started_ns=start, finished_ns=start+1,
            raw={"memory.current": "100", "memory.peak": str(peak), "memory.max": str(128 * 1024**3),
                 "memory.swap.max": "max", "memory.stat": "shmem 64\n",
                 "memory.events": memory["raw"]["memory.events"]})
    before = job_memory(points[0]["host"][0]["started_monotonic_ns"] - 10, 200)
    after = job_memory(memory["finished_ns"] + 10, 300)
    finalization = dict(schema="threadripper_report_finalization_v1", status="reporting_completed",
        job_id=21816, scope=str(job_scope), scientific_timings_admitted=False,
        started_ns=after["finished_ns"]+10, finished_ns=after["finished_ns"]+20,
        reporting_wall_s=10/1e9, job_memory=job_memory(after["finished_ns"]+30, 400))
    ready = dict(pid=123, cgroup=placement["cgroup"], placement=placement)
    for index, point in enumerate(points):
        save(tmp_path / f"point_{index:06d}.json", point)
    samples = iter([sample(done["started_ns"] / 1e9 - .2), sample(done["finished_ns"] / 1e9 + .2)])
    with (tmp_path / "host_processes.jsonl").open("w") as log:
        host = HostMonitor(log, str(job_scope), sample_fn=lambda: next(samples))
        host.observe()
        host.observe()
    summary = host.summary(done["started_ns"] / 1e9, done["finished_ns"] / 1e9)
    completion = completion_evidence(points[0]["thread_affinity"], points[1]["thread_affinity"], 123, done["finished_ns"])
    assert not completion["errors"]
    measured = dict(schema="threadripper_boundary_control_v1", collector_arm="boundary",
        policy=dict(native_points=2, periodic_native_sampling=False,
                    completion_poll_interval_s=1., common_host_interval_s=30.),
        job_id=21816, native=done, native_wall_s=(done["finished_ns"]-done["started_ns"])/1e9,
        status="command_exited_zero", step_memory=memory, placement=placement,
        native_completion=completion, point_records=[record(tmp_path / f"point_{i:06d}.json") for i in range(2)],
        launched=["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1",
            "--cpus-per-task=64", "--cpu-bind=mask_cpu:0xffffffff", sys.executable,
            "-B", str(Path(periodic.__file__).resolve()), "--worker", str(tmp_path)],
        screening=module.evaluate_lineage(points, done, 21816), root_context=module.evaluate_context(points, 21816),
        job_memory=dict(before=before, after=after), host_process_observation=summary,
        scientific_timings_admitted=False, controlled_workload_verified=False,
        native_outputs_validated=False, publication_ready=False)
    for name, value in {"boundary_report": measured, "command": dict(command=["/usr/bin/true"], cpus=32,
        timeout_s=85800, interval_s=1.), "ready": ready, "done": done, "step_memory": memory,
        "go": {"go": True}, "release": {"release": True}, "native_completion": completion,
        "job_memory_before": before, "job_memory_after": after, "report_finalization": finalization,
        "host_process_summary": summary}.items():
        save(tmp_path / (name + ".json"), value)
    return tmp_path


def run(directory, **kwargs):
    return module.replay(directory, 21816, ["/usr/bin/true"], expected_launcher=sys.executable,
        expected_worker=str(Path(periodic.__file__).resolve()), **kwargs)


def test_full_raw_replay_and_relocation(archive):
    result = run(archive)
    assert result["native_outcome"] == "exited_zero"
    assert result["affinity_observation_statuses"] == ["observed_within_affinity"] * 2
    assert result["scientific_timings_admitted"] is False
    assert result["host_process_replay"]["controlled_workload_verified"] is False
    relocated = archive.parent / (archive.name + "-relocated")
    shutil.copytree(archive, relocated)
    assert run(relocated, original_directory=archive)["measured"] == result["measured"]
    with pytest.raises(ValueError, match="launch"):
        run(relocated)


@pytest.mark.parametrize("code", [7, -9, 124])
def test_native_failures_remain_failures(archive, code):
    done = json.loads((archive / "done.json").read_text())
    done["exit_code"] = code
    measured = json.loads((archive / "boundary_report.json").read_text())
    measured.update(native=done, status="command_failed")
    save(archive / "done.json", done)
    save(archive / "boundary_report.json", measured)
    result = run(archive)
    assert result["native_outcome"] == "exited_nonzero" and result["native_exit_code"] == code


@pytest.mark.parametrize("fault", ["schema", "policy", "periodic", "admission", "native_admission",
    "job", "bool_job", "command", "go", "launch", "wall", "exit_bool", "early_timeout",
    "completion", "screening", "context", "point_hash", "missing", "extra_point", "failure",
    "stale_release", "symlink", "host_stream", "host_summary", "memory_scope", "memory_peak",
    "reporting_peak", "reporting_events", "reporting_failed", "initial_affinity", "membership"])
def test_tampering_and_incomplete_evidence_rejected(archive, fault):
    path = archive / "boundary_report.json"
    measured = json.loads(path.read_text())
    if fault == "schema": measured["schema"] = "threadripper_scaling_v5"
    elif fault == "policy": measured["policy"]["native_points"] = 3
    elif fault == "periodic": measured["policy"]["periodic_native_sampling"] = True
    elif fault == "admission": measured["publication_ready"] = True
    elif fault == "native_admission": measured["native_outputs_validated"] = True
    elif fault == "job": measured["job_id"] += 1
    elif fault == "bool_job": measured["job_id"] = True
    elif fault == "command": save(archive / "command.json", {"command": ["/usr/bin/false"]})
    elif fault == "go": save(archive / "go.json", {"go": 1})
    elif fault == "launch": measured["launched"][5] = "--cpus-per-task=32"
    elif fault == "wall": measured["native_wall_s"] += 1.
    elif fault in {"exit_bool", "early_timeout"}:
        done = measured["native"]
        done.update(exit_code=False if fault == "exit_bool" else 124,
                    timed_out=fault == "early_timeout")
        measured["status"] = "command_exited_zero" if fault == "exit_bool" else "command_failed"
        save(archive / "done.json", done)
    elif fault == "completion": measured["native_completion"]["anchor_pid"] += 1
    elif fault == "screening": measured["screening"]["narrow_flagged_intervals"].append(999)
    elif fault == "context": measured["root_context"]["observations"] = 3
    elif fault == "point_hash": measured["point_records"][0]["sha256"] = "0" * 64
    elif fault == "missing": (archive / "job_memory_after.json").unlink()
    elif fault == "extra_point": save(archive / "point_000002.json", {})
    elif fault == "failure": save(archive / "boundary_failure.json", {})
    elif fault == "stale_release": save(archive / "release_freshness_failed.json", {})
    elif fault == "symlink":
        source = archive / "done.json"
        source.rename(archive / "other.json")
        source.symlink_to(archive / "other.json")
    elif fault == "host_stream": (archive / "host_processes.jsonl").write_text("")
    elif fault == "host_summary": measured["host_process_observation"]["successful_snapshots"] += 1
    elif fault in {"memory_scope", "memory_peak"}:
        memory = measured["step_memory"]
        if fault == "memory_scope": memory["scope"] += "_other"
        else: memory["raw"]["memory.peak"] = "99"
        save(archive / "step_memory.json", memory)
    elif fault.startswith("reporting_"):
        source = archive / "report_finalization.json"
        value = json.loads(source.read_text())
        if fault == "reporting_peak": value["job_memory"]["raw"]["memory.peak"] = "250"
        elif fault == "reporting_events": value["job_memory"]["raw"]["memory.events"] += "oom_group_kill 1\n"
        else: value["status"] = "reporting_failed"
        save(source, value)
    else:
        source = archive / "point_000000.json"
        value = json.loads(source.read_text())
        if fault == "initial_affinity": value["thread_affinity"]["finished_ns"] = measured["native"]["started_ns"]
        else: value["native_membership"] = "0::/wrong\n"
        save(source, value)
        measured["point_records"][0] = record(source)
    save(path, measured)
    with pytest.raises((ValueError, KeyError)):
        run(archive)


@pytest.mark.parametrize("change", ["bytes", "inventory"])
def test_mid_replay_change_fails_closed(archive, monkeypatch, change):
    original = module.evaluate_lineage
    def mutate(*args):
        value = original(*args)
        save(archive / ("done.json" if change == "bytes" else "boundary_failure.json"), {})
        return value
    monkeypatch.setattr(module, "evaluate_lineage", mutate)
    with pytest.raises(ValueError):
        run(archive)


@pytest.mark.parametrize("code", [0, 7, -9])
def test_collector_generated_raw_receipts_replay(archive, monkeypatch, code):
    from types import SimpleNamespace
    from benchmark_tools import measure_threadripper_boundary as collector
    from benchmark_tools.report_finalization import observe
    directory = archive / "generated"
    raw = lambda name: json.loads((archive / (name + ".json")).read_text())
    done = raw("done")
    done["exit_code"] = code
    ready = raw("ready")
    points = iter([raw("point_000000"), raw("point_000001")])
    memories = iter([raw("job_memory_before"), raw("job_memory_after"), raw("report_finalization")["job_memory"]])
    monkeypatch.setattr(collector.periodic, "read_point", lambda *a: next(points))
    monkeypatch.setattr(collector.periodic, "read_job_memory", lambda *a: next(memories))
    monkeypatch.setattr(collector.periodic, "step_memory", lambda *a: raw("step_memory"))
    def waiting(path):
        assert path.name == "ready.json"
        save(path, ready)
        return ready
    monkeypatch.setattr(collector, "wait_file", waiting)
    for key, value in ready["placement"]["slurm"].items():
        monkeypatch.setenv(key, value)
    monkeypatch.setattr(collector.os, "uname", lambda: SimpleNamespace(nodename="bizon"))
    class Process:
        finished = False
        def poll(self): return 0 if self.finished else None
        def wait(self, timeout):
            self.finished = True
            return 0
    monkeypatch.setattr(collector.subprocess, "Popen", lambda *a, **k: Process())
    monkeypatch.setattr(collector.time, "sleep", lambda _: save(directory / "done.json", done))
    values = iter([sample(done["started_ns"] / 1e9 - .2), sample(done["finished_ns"] / 1e9 + .2)])
    monkeypatch.setattr(collector, "enriched_snapshot", lambda: next(values))
    class Observer:
        def __init__(self, host, period):
            self.host = host
            assert period == 30.
        def start(self): pass
        def close(self): pass
        def finish(self, start, end):
            self.host.observe()
            return self.host.summary(start, end)
    monkeypatch.setattr(collector, "PeriodicHostObserver", Observer)
    finalization = raw("report_finalization")
    def reporting(*args):
        times = iter([finalization["started_ns"], finalization["finished_ns"]])
        return observe(*args, clock=lambda: next(times))
    monkeypatch.setattr(collector, "observe_reporting", reporting)
    measured = collector.measure(["/usr/bin/true"], directory, 21816, 32, 128 * 1024**3, 85800, 1.)
    result = run(directory)
    assert result["measured"] == measured
    assert result["native_exit_code"] == code
    assert result["host_process_replay"]["successful_snapshots"] == 2
    assert result["native_outputs_validated"] is False
