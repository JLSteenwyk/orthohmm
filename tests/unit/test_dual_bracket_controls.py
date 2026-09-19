import hashlib
from pathlib import Path

import pytest

from benchmark_tools import run_dual_bracket_controls as module
from tests.unit.test_native_pressure_controls import row as pressure_row, rows as pressure_rows


def row(mode="contended"):
    value = pressure_row(mode)
    value["brackets"] = dict(native_pressure=value["pressure"], narrow=dict(
        screen_passed=mode != "contended", reasons=["excess_unassigned_cpu"] if mode == "contended" else []))
    value["points"] = [dict(host=p["host"], native_pressure=p) for p in value["points"]]
    return value


@pytest.mark.parametrize("mode", ["quiet", "native-only", "contended"])
def test_valid_work_enclosure_and_pressure_witness(monkeypatch, mode):
    value = row(mode)
    monkeypatch.setattr(module, "compare", lambda *a, **kw: value["brackets"])
    assert module.validate_trial(value, 123)["status"] == "control_injection_validated"


@pytest.mark.parametrize("fault", ["replay", "delay", "competitor_after", "competitor_before", "exit"])
def test_bad_enclosure_or_execution_rejected(monkeypatch, fault):
    value = row()
    monkeypatch.setattr(module, "compare", lambda *a, **kw: {} if fault == "replay" else value["brackets"])
    if fault == "delay":
        value["native"]["finished_ns"] = 2_900_000_000
    elif fault == "competitor_after":
        value["competitor"]["finished_ns"] = 4_000_000_000
    elif fault == "competitor_before":
        value["competitor"]["started_ns"] = 0
    elif fault == "exit":
        value["worker_exit_code"] = 1
    with pytest.raises(ValueError):
        module.validate_trial(value, 123)


def rows():
    return [{**r, "brackets": dict(native_pressure=r["pressure"], narrow=dict(
        screen_passed=r["mode"] != "contended", reasons=["excess_unassigned_cpu"] if r["mode"] == "contended" else []))}
            for r in pressure_rows()]


def test_summary_requires_all_cpu_and_pressure_checks():
    result = module.summarize(rows())
    assert result["all_control_checks_passed"]
    assert not result["scientific_timings_admitted"]


@pytest.mark.parametrize("fault", ["missing", "order", "failure", "quiet_flag", "no_response", "steal", "pressure"])
def test_bad_panel_not_passed(fault):
    values = rows()
    if fault == "missing":
        values.pop()
    elif fault == "order":
        values.reverse()
    elif fault == "failure":
        values[0] = dict(block=0, mode="quiet", status="failed")
    elif fault == "quiet_flag":
        values[0]["brackets"]["narrow"] = dict(screen_passed=False, reasons=["excess_unassigned_cpu"])
    elif fault == "no_response":
        values[2]["brackets"]["narrow"] = dict(screen_passed=True, reasons=[])
    elif fault == "steal":
        values[2]["brackets"]["narrow"]["reasons"].append("host_steal_time")
    else:
        values[2]["brackets"]["native_pressure"]["native_stall_usec"]["cpu"]["some"] = 1
    if fault in ("missing", "order"):
        with pytest.raises(ValueError):
            module.summarize(values)
    else:
        assert not module.summarize(values)["all_control_checks_passed"]


def test_reader_failure_releases_worker_and_restores_affinity(tmp_path, monkeypatch):
    directory = tmp_path / "trial"
    affinity = []
    monkeypatch.setattr(module.os, "sched_getaffinity", lambda _: {1, 2})
    monkeypatch.setattr(module.os, "sched_setaffinity", lambda _, cpus: affinity.append(set(cpus)))
    class Process:
        def __init__(self, *a, **kw):
            pass
        def wait(self, timeout):
            assert (directory / "go.json").exists()
            assert (directory / "release.json").exists()
            return 0
    monkeypatch.setattr(module.subprocess, "Popen", Process)
    monkeypatch.setattr(module, "wait_file", lambda _: dict(pid=42, cpu=1, membership="test"))
    def fail(*args):
        raise ValueError("counter failure")
    monkeypatch.setattr(module, "read_point", fail)
    with pytest.raises(ValueError, match="counter failure"):
        module.trial(directory, 0, "quiet", 123)
    assert affinity == [{2}, {1, 2}]


def test_protocol_frozen():
    path = Path(module.__file__).parent / "results/DUAL_BRACKET_CONTROL_PROTOCOL_20260919.md"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == module.PROTOCOL_SHA
