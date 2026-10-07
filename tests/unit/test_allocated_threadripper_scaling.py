from copy import deepcopy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import native_factorial_allocated_placement as placement
from benchmark_tools import measure_allocated_threadripper_scaling as collector
from benchmark_tools import replay_allocated_threadripper_scaling as reviewer
from benchmark_tools.native_factorial_cpu_selection import choose
from benchmark_tools.probe_dgx_step_separation import save


ALLOWED = list(range(32, 64))
ALLOCATED = ALLOWED + list(range(128, 160))
CGROUP = "0::/slurm/job_42/step_0/user/task_0\n"


def placement_value():
    topology = [dict(cpu=cpu, package=0, core=cpu % 96, numa_node=0) for cpu in ALLOCATED]
    root = Path("/sys/fs/cgroup")
    current = root / CGROUP[3:].strip().lstrip("/")
    ancestors = []
    while current != root:
        ancestors.append(dict(path=str(current), **{"memory.max": str(128*1024**3),
            "memory.swap.max": "max", "cpu.max": "max 100000",
            "cpuset.cpus.effective": "32-63,128-159"}))
        current = current.parent
    before = dict(host="bizon", pid=1, affinity=list(ALLOCATED),
        topology=[{k: row[k] for k in ("cpu", "package", "core")} for row in topology],
        cgroup=CGROUP, ancestors=ancestors, slurm=dict(SLURM_JOB_ID="42",
        SLURM_CPUS_PER_TASK="64", SLURM_MEM_PER_NODE="131072"))
    after = dict(deepcopy(before), affinity=list(ALLOWED), topology=deepcopy(before["topology"][:32]))
    return dict(schema=placement.SCHEMA, allocated=before, bound=after,
        allocated_topology=topology, selection=choose(list(ALLOCATED), topology),
        started_ns=1, finished_ns=2, sources=placement.sources(), affinity_changed=True,
        scientific_execution_authorized=False, scientific_timings_admitted=False, publication_ready=False)


def affinity_point(start=3, finish=4):
    return dict(native_membership=CGROUP,
        host=[dict(started_monotonic_ns=start)], root_context={"host_after": {"finished_ns": start}},
        thread_affinity=dict(scope="/sys/fs/cgroup/slurm/job_42/step_0/user", allowed_cpus=ALLOWED,
            full_run_affinity_verified=False, started_ns=start, finished_ns=finish,
            initial_tids=[1], final_tids=[1], threads=[dict(tid=1, start_ticks=10,
            cgroup=CGROUP, affinity=[32], outside_cpus=[])], errors=[], violating_tids=[],
            status="observed_within_affinity"))


def test_actual_selection_binds_only_own_worker(monkeypatch):
    value = placement_value()
    snapshots = iter([value["allocated"], value["bound"]])
    monkeypatch.setattr(placement, "inspect", lambda: next(snapshots))
    monkeypatch.setattr(placement, "host_topology", lambda cpus: value["allocated_topology"])
    calls = []
    monkeypatch.setattr(placement.os, "sched_setaffinity", lambda pid, cpus: calls.append((pid, cpus)))
    result = placement.bind(42)
    assert calls == [(0, ALLOWED)]
    assert placement.validate(result, 42) == ALLOWED
    assert result["selection"]["historical_fixed_mask_compatible"] is False
    assert result["scientific_execution_authorized"] is False


def test_incompatible_allocation_never_binds(monkeypatch):
    value = placement_value()
    value["allocated"]["affinity"] = list(range(32))
    monkeypatch.setattr(placement, "inspect", lambda: value["allocated"])
    monkeypatch.setattr(placement, "host_topology", lambda cpus: value["allocated_topology"])
    monkeypatch.setattr(placement.os, "sched_setaffinity", lambda *a: pytest.fail("must not bind"))
    with pytest.raises(ValueError):
        placement.bind(42)


@pytest.mark.parametrize("change", ["schema", "selection", "cpu", "pid", "job", "host", "memory", "cpuset",
    "ancestor", "topology", "before_topology", "source", "window", "bool_window", "admission", "bound_flag",
    "context", "batch", "duplicate", "siblings"])
def test_bad_placement_refused(change):
    value = placement_value()
    before, after = value["allocated"], value["bound"]
    if change == "schema":
        value["schema"] = "historical"
    elif change == "selection":
        value["selection"]["native_cpu_ids"] = list(range(32))
    elif change == "cpu":
        after["affinity"] = list(range(32))
    elif change == "pid":
        after["pid"] = 2
    elif change == "job":
        before["slurm"]["SLURM_JOB_ID"] = after["slurm"]["SLURM_JOB_ID"] = "43"
    elif change == "host":
        before["host"] = after["host"] = "other"
    elif change == "memory":
        for row in [*before["ancestors"], *after["ancestors"]]:
            row["memory.max"] = str(64*1024**3)
    elif change == "cpuset":
        before["ancestors"][0]["cpuset.cpus.effective"] = "0-31"
        after["ancestors"][0]["cpuset.cpus.effective"] = "0-31"
    elif change == "ancestor":
        before["ancestors"].pop()
        after["ancestors"].pop()
    elif change == "topology":
        value["allocated_topology"][0]["core"] = 999
    elif change == "before_topology":
        before["topology"][0]["core"] = 999
    elif change == "source":
        value["sources"][0]["sha256"] = "0"*64
    elif change == "window":
        value["finished_ns"] = 0
    elif change == "bool_window":
        value["started_ns"] = True
    elif change == "admission":
        value["scientific_execution_authorized"] = True
    elif change == "bound_flag":
        value["affinity_changed"] = False
    elif change == "context":
        after["ancestors"][0]["cpu.max"] = "3200000 100000"
    elif change == "batch":
        for sample in (before, after):
            sample["cgroup"] = sample["cgroup"].replace("step_0", "step_batch")
    elif change == "duplicate":
        before["affinity"][1] = before["affinity"][0]
    else:
        value["allocated_topology"][32]["numa_node"] = 1
    with pytest.raises(ValueError):
        placement.validate(value, 42)


@pytest.mark.parametrize("bad", [True, 0, -1, "42"])
def test_bad_job_refused(bad):
    with pytest.raises(ValueError):
        placement.validate(placement_value(), bad)


def test_selected_scope_passed_to_observer(monkeypatch):
    monkeypatch.setattr(collector, "read_root_point", lambda *args: {"root": "unchanged"})
    def observe(scope, cpus):
        assert scope == Path("/sys/fs/cgroup/slurm/job_42/step_0/user")
        assert cpus == ALLOWED
        return dict(status="incomplete", errors=[dict(error="thread exited")])
    monkeypatch.setattr(collector, "observe", observe)
    value = collector.read_point(1, CGROUP, 42, None, ALLOWED)
    assert value["root"] == "unchanged" and value["thread_affinity"]["status"] == "incomplete"
    with pytest.raises(ValueError, match="user subtree"):
        collector.read_point(1, CGROUP.replace("user/", ""), 42, None, ALLOWED)


@pytest.mark.parametrize("cpus", [[], list(range(64)), [32]*32, ALLOWED[::-1], [True]+ALLOWED[1:]])
def test_bad_selected_ids(cpus):
    with pytest.raises(ValueError):
        collector.selected_ids(cpus)


@pytest.mark.parametrize("field,value", [("schema", "threadripper_scaling_v5"), ("cpus", 64),
    ("timeout_s", 900), ("interval_s", True), ("cpu_policy", "fixed_mask"), ("sources", [])])
def test_command_amendment_refuses_tampering(field, value):
    command = collector.command_record(["/bin/sleep", "6"])
    command[field] = value
    with pytest.raises(ValueError):
        collector.validate_command(command)


def test_worker_abort_never_runs_command(tmp_path, monkeypatch):
    save(tmp_path / "command.json", collector.command_record(["/bin/sleep", "6"]))
    monkeypatch.setenv("SLURM_JOB_ID", "42")
    value = placement_value()
    monkeypatch.setattr(collector, "bind", lambda job: value)
    monkeypatch.setattr(collector.os, "getpid", lambda: 1)
    monkeypatch.setattr(collector, "wait_file", lambda path: {"abort": True})
    monkeypatch.setattr(collector, "run_command", lambda *a: pytest.fail("native not released"))
    collector.worker(tmp_path)
    ready = json.loads((tmp_path / "ready.json").read_text())
    assert ready["allocated_placement"]["selection"]["native_cpu_ids"] == ALLOWED
    assert (tmp_path / "aborted_before_native.json").exists()
    assert not (tmp_path / "done.json").exists()


@pytest.mark.parametrize("denied", [None, "budget", "stale", "descendant", "placement"])
def test_collector_release_and_cleanup(tmp_path, monkeypatch, denied):
    directory = tmp_path / "measurement"
    for key, value in {"SLURM_JOB_ID": "42", "SLURM_CPUS_PER_TASK": "64",
                       "SLURM_MEM_PER_NODE": "131072"}.items():
        monkeypatch.setenv(key, value)
    monkeypatch.setattr(collector.os, "uname", lambda: SimpleNamespace(nodename="bizon"))
    events = []
    class Process:
        finished = False
        def poll(self):
            return 0 if self.finished else None
        def wait(self, **kwargs):
            self.finished = True
            events.append("cleanup")
            return 0
    def launch(argv, **kwargs):
        assert "--cpu-bind=cores" in argv and not any("mask_cpu" in arg for arg in argv)
        return Process()
    monkeypatch.setattr(collector.subprocess, "Popen", launch)
    value = placement_value()
    if denied == "placement":
        value["selection"]["native_cpu_ids"] = list(range(32))
    monkeypatch.setattr(collector, "wait_file", lambda *a: dict(pid=1, cgroup=CGROUP,
        placement=value["bound"], allocated_placement=value))
    clock = collector.time.monotonic
    offset = [0.]
    monkeypatch.setattr(collector.time, "monotonic", lambda: clock() + offset[0])
    monkeypatch.setattr(collector, "HostMonitor", lambda *a, **k: SimpleNamespace(
        observe=lambda: None, last_started=clock()))
    started = []
    monkeypatch.setattr(collector, "PeriodicHostObserver", lambda *a: SimpleNamespace(
        start=lambda **k: started.append(True), close=lambda: None, finish=lambda *a: {}))
    monkeypatch.setattr(collector, "read_job_memory", lambda *a: {})
    count = [0]
    def read_point(*args):
        assert args[-1] == ALLOWED
        if denied == "stale":
            offset[0] += 2.
        value = affinity_point(0, 0) if count[0] == 0 else affinity_point(3, 3)
        if count[0] and denied == "descendant":
            value["thread_affinity"]["final_tids"] = [1, 2]
        count[0] += 1
        return value
    monkeypatch.setattr(collector, "read_point", read_point)
    monkeypatch.setattr(collector, "interval_point", lambda *a: {})
    monkeypatch.setattr(collector, "step_memory", lambda *a: {})
    monkeypatch.setattr(collector, "evaluate_lineage", lambda *a: {})
    monkeypatch.setattr(collector, "evaluate", lambda *a: {})
    def sleep(_):
        assert json.loads((directory / "go.json").read_text()) == {"go": True}
        events.append("native")
        save(directory / "done.json", dict(exit_code=0, timed_out=False, started_ns=1, finished_ns=2))
    monkeypatch.setattr(collector.time, "sleep", sleep)
    def guard(path):
        assert path == directory and started == [True] and not (path / "go.json").exists()
        events.append("guard")
        if denied == "budget":
            raise ValueError("insufficient time")
    def invoke():
        return collector.measure(["/native"], directory, 42, 32, 128*1024**3,
                                 85800, 1., release_guard=guard)
    if denied in {"budget", "stale", "placement"}:
        with pytest.raises(ValueError):
            invoke()
        assert json.loads((directory / "go.json").read_text()) == {"abort": True}
        assert not (directory / "done.json").exists()
        assert events[-1] == "cleanup" and "native" not in events
    elif denied == "descendant":
        with pytest.raises(ValueError, match="completion unverified"):
            invoke()
        assert not (directory / "lineage_report.json").exists()
        assert events == ["guard", "native", "cleanup"]
    else:
        result = invoke()
        assert result["schema"] == collector.SCHEMA
        assert result["allocated_placement"]["selection"]["native_cpu_ids"] == ALLOWED
        assert result["native_completion"]["status"] == "anchor_only_at_boundaries"
        assert result["scientific_timings_admitted"] is False
        assert (directory / "report_finalization.json").exists()
        assert events == ["guard", "native", "cleanup"]


def test_replay_preserves_violation_and_gap():
    value = affinity_point()
    assert reviewer.replay_affinity(value, 42, ALLOWED) == "observed_within_affinity"
    a = value["thread_affinity"]
    a["threads"][0].update(affinity=[0], outside_cpus=[0])
    a.update(violating_tids=[1], status="violation")
    assert reviewer.replay_affinity(value, 42, ALLOWED) == "violation"
    a.update(threads=[], violating_tids=[], status="incomplete", errors=[dict(error="exited")])
    assert reviewer.replay_affinity(value, 42, ALLOWED) == "incomplete"


@pytest.mark.parametrize("change", ["scope", "policy", "window", "admission", "gap", "outside",
    "status", "membership", "duplicate", "final_tids", "start_ticks", "errors"])
def test_affinity_replay_refuses_tampering(change):
    value = affinity_point()
    a = value["thread_affinity"]
    if change == "scope":
        a["scope"] += "1"
    elif change == "policy":
        a["allowed_cpus"] = list(range(32))
    elif change == "window":
        a["started_ns"] = 2
    elif change == "admission":
        a["full_run_affinity_verified"] = True
    elif change == "gap":
        a["threads"] = []
    elif change == "outside":
        a["threads"][0]["affinity"] = [0]
    elif change == "status":
        a["status"] = "violation"
    elif change == "membership":
        a["threads"][0]["cgroup"] = CGROUP.replace("step_0", "step_01")
    elif change == "duplicate":
        a["threads"].append(deepcopy(a["threads"][0]))
    elif change == "final_tids":
        a["final_tids"] = [1, 2]
    elif change == "start_ticks":
        a["threads"][0]["start_ticks"] = True
    else:
        a["errors"] = [dict(error="")]
    with pytest.raises(ValueError):
        reviewer.replay_affinity(value, 42, ALLOWED)


def synthetic_measurement(directory, monkeypatch):
    directory.mkdir()
    value = placement_value()
    ready = dict(pid=1, cgroup=CGROUP, placement=value["bound"], allocated_placement=value)
    done = dict(exit_code=0, timed_out=False, started_ns=5, finished_ns=15)
    memory = dict(scope="/slurm/job_42/step_0/user", errors=[], started_ns=30, finished_ns=31,
        raw={"memory.current": "1", "memory.peak": "2",
             "memory.events": "low 0\nhigh 0\nmax 0\noom 0\noom_kill 0\n"})
    def job_memory(start):
        return dict(scope="/slurm/job_42", errors=[], started_ns=start, finished_ns=start+1,
            raw=dict(memory["raw"], **{"memory.max": str(128*1024**3),
                "memory.swap.max": "max", "memory.stat": "shmem 0\n"}))
    before, after = job_memory(3), job_memory(40)
    screening = dict(original_screening={"original_threshold_screen": {"flagged_intervals": []}},
                     narrow_flagged_intervals=[])
    points = collector.DiskObservations(directory)
    points.append(affinity_point(3, 4))
    points.append(affinity_point(17, 19))
    points_finalization = job_memory(50)
    monkeypatch.setattr(reviewer, "evaluate_lineage", lambda *a: screening)
    monkeypatch.setattr(reviewer, "evaluate", lambda *a: {})
    monkeypatch.setattr(reviewer, "interval_point", lambda *a: dict(native_cpu_scope=memory["scope"]))
    monkeypatch.setattr(reviewer, "replay_host", lambda *a: dict(status="synthetic_host_replayed"))
    monkeypatch.setattr(reviewer, "validate_finalization", lambda *a: points_finalization)
    completion = dict(status="anchor_only_at_boundaries")
    monkeypatch.setattr(reviewer, "replay_completion", lambda *a: (
        completion, collector.record(directory / "native_completion.json")))
    command = ["/bin/sleep", "6"]
    measured = dict(schema=collector.SCHEMA, status="command_exited_zero", job_id=42,
        sources=collector.sources(), native=done, native_wall_s=1e-8, native_completion=completion,
        placement=value["bound"], allocated_placement=value,
        launched=["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=64",
                  "--cpu-bind=cores", "/usr/bin/python3", "-B", str(Path(collector.__file__).resolve()),
                  "--worker", str(directory)], point_records=points.records(), step_memory=memory,
        host_process_observation={}, job_memory=dict(before=before, after=after), screening=screening,
        scientific_timings_admitted=False, controlled_workload_verified=False, publication_ready=False)
    save(directory / "lineage_report.json", measured)
    save(directory / "command.json", collector.command_record(command))
    for name, contents in [("done", done), ("step_memory", memory), ("go", {"go": True}),
        ("release", {"release": True}), ("ready", ready), ("job_memory_before", before),
        ("job_memory_after", after), ("host_process_summary", {}), ("native_completion", completion),
        ("report_finalization", {})]:
        save(directory / f"{name}.json", contents)
    (directory / "host_processes.jsonl").write_text("")
    save(directory / "root_context_report.json", dict(status="native_root_context_measured", job_id=42,
        native_wall_s=1e-8, context={}, lineage_report=reviewer.lineage_identity(directory),
        scientific_timings_admitted=False, environmental_validity_established=False))
    return command


def test_new_replay_complete_bindings_and_accounting(tmp_path, monkeypatch):
    directory = tmp_path / "measurement"
    command = synthetic_measurement(directory, monkeypatch)
    value = reviewer.replay(directory, 42, command)
    assert value["schema"] == "allocated_threadripper_replay_v1"
    assert value["native_cpu_ids"] == ALLOWED
    assert value["job_memory_required"] is True
    assert value["affinity_observation_statuses"] == ["observed_within_affinity"]*2
    assert value["scientific_timings_admitted"] is False
    assert value["native_outputs_validated"] is False
    assert value["publication_ready"] is False


@pytest.mark.parametrize("change", ["historical", "sources", "command", "mask", "ready", "admission",
    "memory", "job_memory", "completion_missing", "host_missing", "finalization_missing",
    "abort", "stale", "wrong_expected", "wrong_job", "points", "placement_time"])
def test_new_replay_refuses_mismatched_route(tmp_path, monkeypatch, change):
    directory = tmp_path / "measurement"
    command = synthetic_measurement(directory, monkeypatch)
    target = directory / "lineage_report.json"
    value = json.loads(target.read_text())
    if change == "historical":
        value["schema"] = "threadripper_scaling_v5"
    elif change == "sources":
        value["sources"] = []
    elif change == "command":
        command = ["/bin/sleep", "7"]
    elif change == "mask":
        value["launched"][6] = "--cpu-bind=mask_cpu:0xffffffff"
    elif change == "ready":
        value["placement"]["pid"] = 2
    elif change == "admission":
        value["scientific_timings_admitted"] = True
    elif change == "memory":
        value["step_memory"]["raw"]["memory.peak"] = "0"
    elif change == "job_memory":
        value["job_memory"]["after"]["raw"]["memory.max"] = "max"
    elif change == "completion_missing":
        (directory / "native_completion.json").unlink()
    elif change == "host_missing":
        (directory / "host_processes.jsonl").unlink()
    elif change == "finalization_missing":
        (directory / "report_finalization.json").unlink()
    elif change in {"abort", "stale"}:
        save(directory / ("aborted_before_native.json" if change == "abort" else
                          "release_freshness_failed.json"), {})
    elif change == "wrong_expected":
        command = ["/bin/true"]
    elif change == "points":
        value["point_records"] = []
    elif change == "placement_time":
        value["allocated_placement"]["finished_ns"] = 4
    target.write_text(json.dumps(value))
    with pytest.raises((ValueError, FileNotFoundError)):
        reviewer.replay(directory, 43 if change == "wrong_job" else 42, command)
