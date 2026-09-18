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


def test_live_read_only_probe_and_no_overwrite(tmp_path):
    output = tmp_path / "probe.json"
    result = module.run(output, .01)
    assert result["publication_ready"] is False
    assert len(result["snapshots"]) == 2
    assert result["summary"]["controlled_workload_verified"] is False
    with pytest.raises(FileExistsError):
        module.run(output, .01)
