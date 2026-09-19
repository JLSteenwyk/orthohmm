from copy import deepcopy
import json

import pytest

from benchmark_tools.audit_dgx_pressure import audit, parse_pressure, pressure_summary


def raw(some=100, full=50):
    return (f"some avg10=0.00 avg60=0.00 avg300=0.00 total={some}\n"
            f"full avg10=0.00 avg60=0.00 avg300=0.00 total={full}\n")


def samples():
    before = dict(started_monotonic_ns=0, finished_monotonic_ns=10,
                  raw=dict(boot_id="boot", online_cpus="0-19", cgroup_membership="observer"),
                  optional={f"host_{r}_pressure": raw() for r in ("cpu", "io", "memory")}, errors=[])
    after = deepcopy(before)
    after.update(started_monotonic_ns=1_000_000_000, finished_monotonic_ns=1_000_000_010)
    after["optional"] = {f"host_{r}_pressure": raw(200, 75) for r in ("cpu", "io", "memory")}
    return [before, after]


def test_cumulative_totals_not_moving_average():
    result = pressure_summary(samples(), "io")
    assert result["stall_usec"] == {"some": 100, "full": 25}
    assert result["midpoint_percent"] == {"some": .01, "full": .0025}
    assert result["midpoint_span_s"] == 1


def test_system_cpu_full_not_interpreted():
    assert parse_pressure(raw(), "cpu") == {"some": 100}
    assert pressure_summary(samples(), "cpu")["cpu_full_undefined"] is True


@pytest.mark.parametrize("value", ["", raw() + raw(), raw().replace("total=100", "total=-1"),
    raw().replace("avg10=0.00", "avg10=nan"), raw().replace("avg60=0.00", "avg60=101"),
    raw().replace("avg300=0.00", "avg10=0.00"), raw().replace("total=100", "total=1.1")])
def test_invalid_psi_rejected(value):
    with pytest.raises(ValueError):
        parse_pressure(value, "io")


@pytest.mark.parametrize("fault", ["missing", "read_error", "boot", "scope", "time", "decrease", "reset_then_recover"])
def test_invalid_sample_series_rejected(fault):
    points = samples()
    if fault == "missing":
        del points[1]["optional"]["host_memory_pressure"]
    elif fault == "read_error":
        points[1]["errors"] = [dict(field="host_memory_pressure")]
    elif fault in ("boot", "scope"):
        points[1]["raw"]["boot_id" if fault == "boot" else "cgroup_membership"] = "changed"
    elif fault == "time":
        points[1]["started_monotonic_ns"] = 0
    elif fault == "decrease":
        points[1]["optional"]["host_memory_pressure"] = raw(10, 5)
    else:
        middle = deepcopy(points[0])
        middle.update(started_monotonic_ns=100, finished_monotonic_ns=110)
        middle["optional"]["host_memory_pressure"] = raw(0, 0)
        points.insert(1, middle)
    with pytest.raises(ValueError):
        pressure_summary(points, "memory")


def test_unrelated_optional_read_failure_not_silently_applied_to_psi():
    points = samples()
    points[0]["errors"] = [dict(field="cgroup_memory.peak")]
    assert pressure_summary(points, "io")["stall_usec"]["some"] == 100


def test_cpu_some_only_is_valid_but_memory_requires_full():
    value = raw().splitlines()[0]
    assert parse_pressure(value, "cpu") == {"some": 100}
    with pytest.raises(ValueError):
        parse_pressure(value, "memory")


def test_terminal_gate_precedes_measurement_access(tmp_path):
    (tmp_path / "accounting_terminal.txt").write_text("21838_0|RUNNING|0:0|00:00:01|20|96G|spark-7ff0\n")
    with pytest.raises(ValueError, match="nonterminal"):
        audit(tmp_path)


def test_complete_archive_retains_partial_and_missing_resource(tmp_path):
    (tmp_path / "accounting_terminal.txt").write_text("".join(
        f"21838_{i}|COMPLETED|0:0|00:00:01|20|96G|spark-7ff0\n" for i in range(18)))
    for i in range(18):
        directory = tmp_path / "frontier_overhead_v1" / f"run_{i:02d}" / "measurement"
        directory.mkdir(parents=True)
        points = samples()
        if i == 0:
            del points[-1]["optional"]["host_io_pressure"]
        (directory / "point_0000.json").write_text(json.dumps({"host": points}))
        (directory / "done.json").write_text(json.dumps(dict(
            exit_code=0, started_ns=100, finished_ns=2_000_000_000 if i == 1 else 900_000_000)))
    result = audit(tmp_path)
    assert len(result["runs"]) == 18
    assert len(result["records"]) == 37
    assert result["runs"][0]["resources"]["io"]["status"] == "unavailable"
    assert "stall_usec" not in result["runs"][0]["resources"]["io"]
    assert result["runs"][0]["resources"]["memory"]["stall_usec"]["some"] == 100
    assert result["runs"][1]["recorded_window_encloses_native"] is False
    assert result["runs"][2]["recorded_window_encloses_native"] is True
    assert result["scientific_timings_admitted"] is False
