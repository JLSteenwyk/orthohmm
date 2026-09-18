import json
from pathlib import Path

import pytest

from benchmark_tools.screen_bracketed_cpu import screen, usage
from benchmark_tools.probe_host_counters import parse_cpu


def sample(start, scope, busy=0, usage_usec=0):
    raw = "cpu " + " ".join(map(str, [busy, 0, 0, 1000, 0, 0, 0, 0, 0, 0])) + "\n"
    return dict(started_monotonic_ns=start, finished_monotonic_ns=start + 100,
        raw=dict(proc_stat=raw, cgroup_membership=f"0::/job_123/step_{scope}\n", boot_id="boot", online_cpus="0-19"),
        cpu_ticks=parse_cpu(raw), errors=[], optional={"cgroup_cpu.stat": f"usage_usec {usage_usec}\n"})


@pytest.fixture
def samples():
    return [sample(0, "batch"), sample(1000, "0"),
            sample(10_000_010_000, "0", busy=10000, usage_usec=100_000_000),
            sample(10_000_020_000, "batch", busy=10000)]


def evaluate(samples):
    return screen(*samples, 123, 100, 2000, 10_000_002_000)


def test_quiet_arithmetic_not_admission(samples):
    result = evaluate(samples)
    assert result["screen_passed"] is True
    assert result["signed_unassigned_cpu_s"] == 0
    assert result["controlled_workload_verified"] is False
    assert result["scientific_timings_admitted"] is False


def test_completed_burst_is_retained(samples):
    samples[-1] = sample(10_000_020_000, "batch", busy=10300)
    result = evaluate(samples)
    assert result["reasons"] == ["excess_unassigned_cpu"]
    assert result["signed_unassigned_cpu_s"] == 3


def test_negative_residual_is_not_clipped(samples):
    samples[2]["optional"]["cgroup_cpu.stat"] = "usage_usec 101000000\n"
    result = evaluate(samples)
    assert result["reasons"] == ["negative_accounting_discrepancy"]
    assert result["signed_unassigned_cpu_s"] == -1


@pytest.mark.parametrize("failure", ["left", "right", "boot", "scope", "error", "counter", "work"])
def test_invalid_evidence(samples, failure):
    if failure == "left":
        samples[0]["finished_monotonic_ns"] = 2000
    elif failure == "right":
        samples[-1]["started_monotonic_ns"] = samples[2]["started_monotonic_ns"]
    elif failure == "boot":
        samples[1]["raw"]["boot_id"] = "other"
    elif failure == "scope":
        samples[1]["raw"]["cgroup_membership"] = "0::/job_999/step_0\n"
    elif failure == "error":
        samples[0]["errors"] = ["missing"]
    elif failure == "counter":
        samples[1]["optional"]["cgroup_cpu.stat"] = "usage_usec 100000001\n"
    else:
        samples[1]["finished_monotonic_ns"] = 2001
    with pytest.raises(ValueError):
        evaluate(samples)


@pytest.mark.parametrize("raw", ["", "usage_usec 1\nusage_usec 2\n", "usage_usec -1\n", "usage_usec 1 extra\n"])
def test_invalid_usage(raw):
    with pytest.raises(ValueError):
        usage({"optional": {"cgroup_cpu.stat": raw}})


def test_existing_smokes_cannot_be_retroactively_screened():
    root = Path(__file__).resolve().parents[2]
    report = json.loads((root / "benchmark_tools/results/dgx_native_window_rejections_20260918.json").read_text())
    assert [row["index"] for row in report["trials"]] == [0, 1, 2]
    for row in report["trials"]:
        host_before, before, after, host_after = row["snapshots"]
        assert host_before["finished_monotonic_ns"] - before["started_monotonic_ns"] == row["left_overlap_ns"]
        assert row["left_overlap_ns"] > 0
        with pytest.raises(ValueError, match="fully bracket"):
            screen(host_before, before, after, host_after, row["job_id"], report["clock_ticks_per_second"],
                   row["work_started_ns"], row["work_finished_ns"])


def test_steal_is_not_native_cpu(samples):
    samples[-1]["raw"]["proc_stat"] = "cpu 10000 0 0 1000 0 0 0 1 0 0\n"
    samples[-1]["cpu_ticks"] = parse_cpu(samples[-1]["raw"]["proc_stat"])
    assert evaluate(samples)["reasons"] == ["host_steal_time"]


def test_small_signed_discrepancy_is_reported(samples):
    samples[2]["optional"]["cgroup_cpu.stat"] = "usage_usec 100100000\n"
    result = evaluate(samples)
    assert result["screen_passed"] is True
    assert result["signed_unassigned_cpu_s"] == pytest.approx(-.1)
