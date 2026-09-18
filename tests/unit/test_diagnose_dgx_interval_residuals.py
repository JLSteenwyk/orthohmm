import json
from pathlib import Path

import pytest

from benchmark_tools.diagnose_dgx_interval_residuals import decompose, delta_stat, diagnose, parse_stat
from benchmark_tools.probe_host_counters import parse_cpu


def sample(start, ticks, usage):
    raw = f"cpu {ticks} 0 0 1000 0 0 0 0 0 0\n"
    return dict(started_monotonic_ns=start, finished_monotonic_ns=start+100,
                raw=dict(proc_stat=raw, cgroup_membership="0::/job_123/step_batch\n",
                         boot_id="boot", online_cpus="0-19"),
                cpu_ticks=parse_cpu(raw), errors=[],
                optional={"cgroup_cpu.stat": f"usage_usec {usage}\nuser_usec {usage}\nsystem_usec 0\n"})


def point(i):
    start = i * 1_000_000_000
    return dict(host=[sample(start, i*100, i*1000), sample(start+400, i*100+2+i, i*1000+1000)],
                native_membership="0::/job_123/step_0/task_0\n", native_cpu_scope="/job_123/step_0",
                native_read_ns=[start+200, start+300], ticks_per_second=100,
                native_cpu_stat=f"usage_usec {i*1000000}\nuser_usec {i*900000}\nsystem_usec {i*100000}\n")


def test_exact_window_algebra_and_negative_inner_residual_retained():
    row = decompose(point(0), point(1), 123)
    assert row["host_outer_busy_cpu_s"] == 1.03
    assert row["host_inner_busy_cpu_s"] == .98
    assert row["read_window_busy_cpu_s"] == pytest.approx(.05)
    assert row["inner_host_minus_native_cpu_s"] == pytest.approx(-.02)
    assert row["observer_leaf_delta_cpu_s"]["usage_usec"] == .002
    assert row["outer_host_minus_native_minus_observer_cpu_s"] == pytest.approx(.028)
    assert row["native_step_delta_cpu_s"] == dict(usage_usec=1., user_usec=.9, system_usec=.1)
    assert row["scientific_timings_admitted"] is False


@pytest.mark.parametrize("raw", ["usage_usec 1\n", "usage_usec 1\nusage_usec 2\n",
                                 "usage_usec -1\nuser_usec 0\nsystem_usec 0\n",
                                 "usage_usec nan\nuser_usec 0\nsystem_usec 0\n"])
def test_invalid_cpu_stats_rejected(raw):
    with pytest.raises(ValueError):
        parse_stat(raw)


def test_user_counter_decrease_rejected_even_when_usage_increases():
    with pytest.raises(ValueError):
        delta_stat("usage_usec 5\nuser_usec 4\nsystem_usec 1\n",
                   "usage_usec 6\nuser_usec 3\nsystem_usec 3\n")


@pytest.mark.parametrize("fault", ["scope", "errors", "raw", "observer_decrease"])
def test_original_validity_gates_preserved(fault):
    left, right = point(0), point(1)
    if fault == "scope":
        right["native_cpu_scope"] = "/job_123/step_batch"
    elif fault == "errors":
        right["host"][0]["errors"] = ["failed"]
    elif fault == "raw":
        right["host"][0]["cpu_ticks"]["user"] += 1
    else:
        left["host"][0]["optional"]["cgroup_cpu.stat"] = "usage_usec 999999\nuser_usec 999999\nsystem_usec 0\n"
    with pytest.raises(ValueError):
        decompose(left, right, 123)


def test_real_counter_control_all_intervals_independent_tick_arithmetic():
    root = Path(__file__).resolve().parents[2]
    report = json.loads((root / "benchmark_tools/results/dgx_interval_controls_21806.json").read_text())
    for trial in report["trials"]:
        for left, right in zip(trial["points"], trial["points"][1:]):
            row = decompose(left, right, report["job_id"])
            def busy(snapshot):
                ticks = [int(v) for v in snapshot["raw"]["proc_stat"].splitlines()[0].split()[1:]]
                return sum(ticks[i] for i in (0, 1, 2, 5, 6))
            mass = (busy(right["host"][1]) - busy(right["host"][0])
                    + busy(left["host"][1]) - busy(left["host"][0])) / left["ticks_per_second"]
            assert row["read_window_busy_cpu_s"] == pytest.approx(mass)


def test_stored_screen_mismatch_rejected(monkeypatch):
    import benchmark_tools.diagnose_dgx_interval_residuals as module
    monkeypatch.setattr(module, "evaluate", lambda *args: {"flagged_intervals": [3]})
    with pytest.raises(ValueError, match="Stored screening"):
        diagnose(dict(points=[], native={}, job_id=123, screening={"flagged_intervals": []}))


def test_retained_smoke_decomposition_preserves_adverse_evidence():
    root = Path(__file__).resolve().parents[2]
    report = json.loads((root / "benchmark_tools/results/dgx_interval_residual_decomposition_20260918.json").read_text())
    assert report["scientific_timings_admitted"] is False
    assert [r["original_flagged_intervals"] for r in report["runs"]] == [[], [3], [3]]
    for run in report["runs"]:
        for row in run["intervals"]:
            residual = row["original_screen"]["signed_unassigned_cpu_s"]
            assert row["inner_host_minus_native_cpu_s"] == pytest.approx(residual-row["read_window_busy_cpu_s"])
            assert row["outer_host_minus_native_minus_observer_cpu_s"] == pytest.approx(
                residual-row["observer_leaf_delta_cpu_s"]["usage_usec"])
    for run in report["runs"][1:]:
        row = run["intervals"][3]
        assert row["read_window_busy_cpu_s"] == pytest.approx(.01)
        assert row["observer_leaf_delta_cpu_s"]["usage_usec"] < .004
        assert row["inner_host_minus_native_cpu_s"] > .25
