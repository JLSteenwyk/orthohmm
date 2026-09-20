import json
import shutil
from types import SimpleNamespace

import pytest

from benchmark_tools import replay_scaling_root_context as module
from tests.unit.test_measure_native_hierarchy_step import evidence
from tests.unit.test_measure_native_lineage_step import extend_lineage
from tests.unit.test_native_root_context import add_context


def save(path, value):
    path.write_text(json.dumps(value))


@pytest.fixture
def archive(tmp_path, evidence):
    points, done = extend_lineage(evidence)
    add_context(points)
    end = points[-1]["root_context"]["host_after"]["finished_ns"]
    memory = dict(scope=module.interval_point(points[-1], 21816)["native_cpu_scope"], errors=[],
        started_ns=end+10, finished_ns=end+20, raw={"memory.current": "100", "memory.peak": "200",
            "memory.events": "low 0\nhigh 0\nmax 0\noom 0\noom_kill 0\n"})
    measured = dict(job_id=21816, native=done, step_memory=memory, points=points,
        scientific_timings_admitted=False, controlled_workload_verified=False, publication_ready=False,
        status="command_exited_zero", native_wall_s=(done["finished_ns"]-done["started_ns"])/1e9,
        screening=module.evaluate_lineage(points, done, 21816))
    save(tmp_path / "command.json", dict(command=["/usr/bin/true"], cpus=20, timeout_s=85800, interval_s=1.))
    return tmp_path, measured


def write(archive):
    directory, measured = archive
    save(directory / "go.json", {"go": True})
    save(directory / "release.json", {"release": True})
    save(directory / "lineage_report.json", measured)
    save(directory / "done.json", measured["native"])
    save(directory / "step_memory.json", measured["step_memory"])
    for index, point in enumerate(measured["points"]):
        save(directory / f"point_{index:06d}.json", point)
    save(directory / "root_context_report.json", dict(status="native_root_context_measured", job_id=21816,
        native_wall_s=measured["native_wall_s"], context=module.evaluate(measured["points"], 21816),
        lineage_report=module.lineage_identity(directory), scientific_timings_admitted=False,
        environmental_validity_established=False))


@pytest.mark.parametrize("code", [0, 7, -9, 124])
def test_complete_raw_replay_retains_native_failures_and_relocates(archive, code):
    directory, measured = archive
    measured["native"]["exit_code"] = code
    measured["status"] = "command_exited_zero" if code == 0 else "command_failed"
    measured["step_memory"]["raw"]["memory.events"] += "oom_group_kill 1\n"
    write(archive)
    result = module.replay(directory, 21816, ["/usr/bin/true"])
    assert result["native_outcome"] == ("exited_zero" if code == 0 else "exited_nonzero")
    assert result["native_exit_code"] == code
    assert result["memory_events"]["oom_group_kill"] == 1
    assert result["screening"] == measured["screening"]
    assert result["native_outputs_validated"] is False and result["scientific_timings_admitted"] is False
    relocated = directory.parent / (directory.name + "_relocated")
    shutil.copytree(directory, relocated)
    other = module.replay(relocated, 21816, ["/usr/bin/true"])
    assert {k:v for k,v in other.items() if k not in {"evidence", "source"}} == {
        k:v for k,v in result.items() if k not in {"evidence", "source"}}


@pytest.mark.parametrize("fault", ["timeout", "command", "cadence", "job", "bool_job", "exit_bool",
    "timeout_bool", "early_timeout", "status", "wall", "memory_scope", "memory_time", "memory_events",
    "memory_gauge", "admission", "screening", "gap", "extra", "legacy_names", "raw", "context", "hash", "symlink",
    "go", "release", "abort_marker"])
def test_invalid_records_cannot_become_valid_failures(archive, fault):
    directory, measured = archive
    if fault == "job": measured["job_id"] += 1
    elif fault == "bool_job": measured["job_id"] = True
    elif fault == "exit_bool": measured["native"]["exit_code"] = False
    elif fault == "timeout_bool": measured["native"]["timed_out"] = 0
    elif fault == "early_timeout":
        measured["native"].update(exit_code=124, timed_out=True)
        measured["status"] = "command_failed"
    elif fault == "status": measured["status"] = "command_failed"
    elif fault == "wall": measured["native_wall_s"] += 1
    elif fault == "memory_scope": measured["step_memory"]["scope"] += "/other"
    elif fault == "memory_time": measured["step_memory"]["started_ns"] = measured["native"]["finished_ns"]
    elif fault == "memory_events": measured["step_memory"]["raw"]["memory.events"] += "oom 0\n"
    elif fault == "memory_gauge": measured["step_memory"]["raw"]["memory.peak"] = "99"
    elif fault == "admission": measured["publication_ready"] = True
    elif fault == "screening": measured["screening"]["narrow_flagged_intervals"].append(999)
    write(archive)
    if fault in {"timeout", "command", "cadence"}:
        path = directory / "command.json"
        data = json.loads(path.read_text())
        key, value = {"timeout": ("timeout_s", 900), "command": ("command", ["/usr/bin/false"]),
                      "cadence": ("interval_s", True)}[fault]
        data[key] = value
        save(path, data)
    elif fault == "gap": (directory / "point_000001.json").rename(directory / "point_000002.json")
    elif fault == "extra": save(directory / "point_extra.json", {})
    elif fault == "legacy_names": (directory / "point_000001.json").rename(directory / "point_0001.json")
    elif fault == "raw": save(directory / "done.json", {})
    elif fault == "go": save(directory / "go.json", {"go": 1})
    elif fault == "release": save(directory / "release.json", {"release": False})
    elif fault == "abort_marker": save(directory / "aborted_before_native.json", {})
    elif fault in {"context", "hash"}:
        path = directory / "root_context_report.json"
        data = json.loads(path.read_text())
        if fault == "hash": data["lineage_report"]["sha256"] = "bad"
        else: data["context"]["intervals"][0]["root_minus_system_cpu_usec"] = 999
        save(path, data)
    elif fault == "symlink":
        path = directory / "done.json"
        path.rename(directory / "done_original.json")
        path.symlink_to(directory / "done_original.json")
    with pytest.raises((ValueError, KeyError)):
        module.replay(directory, 21816, ["/usr/bin/true"])


@pytest.mark.parametrize("code,timed_out,wall,accepted", [(124,True,85800.1,True),
    (124,True,1.,False), (0,True,85800.1,False), (7,False,85831.,False), (0,False,0.,False)])
def test_timeout_and_cleanup_boundary(code, timed_out, wall, accepted):
    done = dict(exit_code=code, timed_out=timed_out, started_ns=1_000_000_000,
                finished_ns=1_000_000_000+int(wall*1e9))
    measured = dict(status="command_exited_zero" if code == 0 else "command_failed", native_wall_s=wall)
    if accepted:
        assert module.native_outcome(measured, done) == ("timed_out", wall)
    else:
        with pytest.raises(ValueError): module.native_outcome(measured, done)


@pytest.mark.parametrize("change_kind", ["bytes", "added_marker"])
def test_raw_change_during_replay_rejected(archive, monkeypatch, change_kind):
    write(archive)
    original = module.evaluate_lineage
    def change(*args):
        result = original(*args)
        save(archive[0] / ("done.json" if change_kind == "bytes" else "failed_point.json"), {})
        return result
    monkeypatch.setattr(module, "evaluate_lineage", change)
    with pytest.raises(ValueError): module.replay(archive[0], 21816, ["/usr/bin/true"])


def test_long_point_inventory_uses_numeric_order_without_loading_raw_content(tmp_path):
    class Directory:
        def __truediv__(self, name): return tmp_path / name
        def glob(self, pattern): return reversed([tmp_path / f"point_{i:06d}.json" for i in range(10002)])
    paths = module.point_inventory(Directory())
    assert len(paths) == 10002 and paths[-1].name == "point_010001.json"


@pytest.mark.parametrize("code", [0, 7, -9])
def test_collector_generated_reports_replay_with_real_evaluators(archive, monkeypatch, code):
    from benchmark_tools import measure_scaling_root_context as collector
    root, source = archive
    directory = root / "composed"
    done = dict(source["native"], exit_code=code, timed_out=False)
    for key, value in {"SLURM_JOB_ID": "21816", "SLURM_CPUS_PER_TASK": "20", "SLURM_MEM_PER_NODE": "98304"}.items():
        monkeypatch.setenv(key, value)
    monkeypatch.setattr(collector.os, "uname", lambda: SimpleNamespace(nodename="spark-7ff0"))
    class Process:
        finished = False
        def poll(self): return 0 if self.finished else None
        def wait(self, timeout):
            assert (directory / "release.json").exists()
            self.finished = True
            return 0
    monkeypatch.setattr(collector.subprocess, "Popen", lambda *a, **k: Process())
    monkeypatch.setattr(collector, "wait_file", lambda p: {"pid": 123, "cgroup": source["points"][0]["native_membership"]})
    pending = iter(source["points"])
    monkeypatch.setattr(collector, "read_point", lambda *a: next(pending))
    monkeypatch.setattr(collector, "step_memory", lambda p: source["step_memory"])
    monkeypatch.setattr(collector.time, "sleep", lambda delay: save(directory / "done.json", done))
    measured = collector.measure(["/usr/bin/true"], directory, 21816, 20, 96*1024**3, 85800, 1.)
    result = module.replay(directory, 21816, ["/usr/bin/true"])
    assert result["measured"] == measured
    assert result["native_exit_code"] == code
    assert result["screening"] == measured["screening"]
    assert result["scientific_timings_admitted"] is False
