import errno
import json
import sys

import pytest

from benchmark_tools import probe_host_counters as module


def sample(start, cpu):
    text = "cpu " + " ".join(map(str, cpu)) + "\ncpu0 0 0\n"
    return dict(started_monotonic_ns=start, finished_monotonic_ns=start + 100,
        raw=dict(proc_stat=text, boot_id="boot", online_cpus="0-19", cgroup_membership="0::/test"),
        cpu_ticks=module.parse_cpu(text))


def test_busy_counter_does_not_double_count_guest_or_include_steal():
    first = sample(0, [0] * 10)
    last = sample(1_000_000_000, [100, 10, 20, 500, 4, 3, 2, 7, 30, 5])
    report = module.summarize(first, last, 100)
    assert report["accounted_host_busy_cpu_s"] == 1.35
    assert report["accounted_host_busy_average_cores"] == 1.35
    assert report["cpu_delta_ticks"]["steal"] == 7
    assert report["controlled_workload_verified"] is False


def test_iowait_decrease_retained_not_clamped():
    first = sample(0, [100] * 10)
    last = sample(1_000_000_000, [110, 110, 110, 110, 90, 110, 110, 110, 110, 110])
    result = module.summarize(first, last, 100)
    assert result["iowait_decreased"] is True
    assert result["cpu_delta_ticks"]["iowait"] == -10


@pytest.mark.parametrize("problem", ["boot_id", "online_cpus", "cgroup_membership", "parsed", "time", "counter", "hz"])
def test_invalid_comparisons(problem):
    first, last = sample(0, [100] * 10), sample(1_000_000_000, [110] * 10)
    if problem in ("boot_id", "online_cpus", "cgroup_membership"):
        last["raw"][problem] = "changed"
    elif problem == "parsed":
        last["cpu_ticks"]["user"] += 1
    elif problem == "time":
        last["started_monotonic_ns"] = 1
    elif problem == "counter":
        last = sample(1_000_000_000, [90] * 10)
    with pytest.raises(ValueError):
        module.summarize(first, last, True if problem == "hz" else 100)


@pytest.mark.parametrize("text", ["", "cpu 1 2", "cpu " + "-1 " * 10, "cpu " + "x " * 10])
def test_bad_cpu_rows(text):
    with pytest.raises(ValueError):
        module.parse_cpu(text)


@pytest.mark.parametrize("text", ["0::relative", "0::/../test", "1:cpu:/test", "0::/a\n0::/b"])
def test_bad_membership(text):
    with pytest.raises(ValueError):
        module.parse_group(text)


@pytest.mark.skipif(sys.platform != "linux", reason="Linux /proc integration; executed in linux-native-diagnostics")
def test_live_read_only_probe_and_no_overwrite(tmp_path):
    output = tmp_path / "probe.json"
    result = module.run(output, .01)
    assert result["publication_ready"] is False
    assert len(result["snapshots"]) == 2
    assert result["summary"]["controlled_workload_verified"] is False
    with pytest.raises(FileExistsError):
        module.run(output, .01)


@pytest.fixture
def synthetic_reader(monkeypatch):
    membership = ["0::/test", "0::/test"]
    required = {
        "/proc/stat": "cpu " + " ".join(["100"] * 10) + "\ncpu0 0 0\n",
        "/sys/devices/system/cpu/online": "0-19\n",
        "/proc/sys/kernel/random/boot_id": "boot\n",
    }
    optional = {
        **{f"/proc/pressure/{r}": "some total=1\n" for r in ("cpu", "memory", "io")},
        **{f"/sys/fs/cgroup/test/{n}": "synthetic\n"
           for n in ("cpu.stat", "memory.current", "memory.peak", "memory.events")},
    }
    reads = []
    def read(path, *args, **kwargs):
        name = str(path)
        reads.append(name)
        if name == "/proc/self/cgroup":
            return membership.pop(0)
        value = required.get(name, optional.get(name))
        if value is None:
            raise AssertionError("Unexpected filesystem read: " + name)
        if isinstance(value, Exception):
            raise value
        return value
    clock = iter([100, 200])
    monkeypatch.setattr(module.Path, "read_text", read)
    monkeypatch.setattr(module.time, "monotonic_ns", lambda: next(clock))
    return required, optional, membership, reads


def test_snapshot_preserves_raw_optional_and_read_bracket(synthetic_reader):
    required, optional, membership, reads = synthetic_reader
    result = module.snapshot()
    assert result["started_monotonic_ns"] == 100
    assert result["finished_monotonic_ns"] == 200
    assert result["raw"] == dict(proc_stat=required["/proc/stat"], boot_id="boot\n",
                                online_cpus="0-19\n", cgroup_membership="0::/test")
    assert result["cpu_ticks"] == dict.fromkeys(module.CPU_FIELDS, 100)
    assert len(result["optional"]) == len(optional) == 7
    assert result["optional"]["cgroup_memory.peak"] == "synthetic\n"
    assert result["errors"] == []
    assert not membership
    assert reads.count("/proc/self/cgroup") == 2


@pytest.mark.parametrize("error", [FileNotFoundError(errno.ENOENT, "missing"),
                                   PermissionError(errno.EACCES, "denied")])
def test_optional_read_failure_is_recorded_not_fabricated(synthetic_reader, error):
    _, optional, _, _ = synthetic_reader
    optional["/sys/fs/cgroup/test/memory.peak"] = error
    result = module.snapshot()
    assert "cgroup_memory.peak" not in result["optional"]
    assert len(result["optional"]) == 6
    assert result["errors"] == [dict(field="cgroup_memory.peak", type=type(error).__name__, errno=error.errno)]


@pytest.mark.parametrize("path", ["/proc/stat", "/sys/devices/system/cpu/online",
                                  "/proc/sys/kernel/random/boot_id"])
def test_missing_required_host_file_fails_closed(synthetic_reader, path):
    required, _, _, _ = synthetic_reader
    required[path] = FileNotFoundError(errno.ENOENT, "missing")
    with pytest.raises(FileNotFoundError):
        module.snapshot()


@pytest.mark.parametrize("fault", ["changed", "invalid", "initial_invalid", "cpu"])
def test_snapshot_rejects_changed_scope_or_malformed_evidence(synthetic_reader, fault):
    required, _, membership, _ = synthetic_reader
    if fault == "cpu":
        required["/proc/stat"] = "cpu 1 2\n"
    elif fault == "initial_invalid":
        membership[0] = "1:cpu:/test"
    else:
        membership[1] = "0::/different" if fault == "changed" else "0::relative"
    with pytest.raises(ValueError):
        module.snapshot()


@pytest.mark.parametrize("interval", [0, -1, 61, float("inf"), float("nan")])
def test_invalid_interval_never_reads_host(monkeypatch, tmp_path, interval):
    def forbidden():
        raise AssertionError("Host read was not expected")
    monkeypatch.setattr(module, "snapshot", forbidden)
    output = tmp_path / "invalid.json"
    with pytest.raises(ValueError):
        module.run(output, interval)
    assert not output.exists()


def test_run_orchestration_with_synthetic_snapshots(monkeypatch, tmp_path):
    samples = iter([sample(0, [100] * 10), sample(1_000_000_000, [110] * 10)])
    sleeps = []
    monkeypatch.setattr(module, "snapshot", lambda: next(samples))
    monkeypatch.setattr(module.time, "sleep", sleeps.append)
    monkeypatch.setattr(module.os, "sysconf", lambda name: 100 if name == "SC_CLK_TCK" else None)
    monkeypatch.setattr(module.platform, "node", lambda: "synthetic-host")
    monkeypatch.setattr(module.platform, "release", lambda: "synthetic-kernel")
    output = tmp_path / "synthetic.json"
    result = module.run(output, .01)
    assert sleeps == [.01]
    assert result["hostname"] == "synthetic-host"
    assert result["kernel"] == "synthetic-kernel"
    assert result["clock_ticks_per_second"] == 100
    assert result["summary"]["accounted_host_busy_cpu_s"] == .5
    assert result["summary"]["controlled_workload_verified"] is False
    assert result["publication_ready"] is False
    assert json.loads(output.read_text()) == result
    with pytest.raises(FileExistsError):
        module.run(output, .01)


def test_source_mutation_prevents_success_receipt(monkeypatch, tmp_path):
    samples = iter([sample(0, [100] * 10), sample(1_000_000_000, [110] * 10)])
    source = iter([b"first source", b"changed source"])
    monkeypatch.setattr(module, "snapshot", lambda: next(samples))
    monkeypatch.setattr(module.time, "sleep", lambda interval: None)
    monkeypatch.setattr(module.os, "sysconf", lambda name: 100)
    monkeypatch.setattr(module.platform, "node", lambda: "synthetic-host")
    monkeypatch.setattr(module.platform, "release", lambda: "synthetic-kernel")
    monkeypatch.setattr(module.Path, "read_bytes", lambda path: next(source))
    output = tmp_path / "changed.json"
    with pytest.raises(ValueError, match="source changed"):
        module.run(output, .01)
    assert not output.exists()
