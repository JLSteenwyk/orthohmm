import json

import pytest

from benchmark_tools import measure_native_root_context as module
from benchmark_tools.measure_native_scaling_run import measure_run
from tests.unit.test_measure_native_scaling_run import fixture


@pytest.mark.parametrize("status", ["command_exited_zero", "command_failed"])
def test_native_workflow_retains_report_contract_and_input_checks(fixture, monkeypatch, status):
    run, order, root = fixture
    calls = []

    def collect(command, directory, job, cpus, memory, timeout, interval):
        assert (root / "preparation.json").exists()
        assert command[-3:] == run["native_argv"]
        assert (job, cpus, memory, timeout, interval) == (123, 20, 96*1024**3, 900, 1.)
        calls.append(command)
        return dict(status=status, sentinel="original lineage report"), dict(context="separate report")

    monkeypatch.setattr(module, "measure", collect)
    result = measure_run(run, order, [], lambda _: ["a.fa"], module.measure_native_run,
                         123, timeout_s=900, runtime_checker=lambda _: [])
    assert len(calls) == 1
    assert result["status"] == status
    assert result["measurement"]["sentinel"] == "original lineage report"
    assert result["after"]["native_order"] == ["a.fa"]
    assert json.loads((root / "verification.json").read_text())["measurement"] == result["measurement"]


@pytest.mark.parametrize("field,value", [("cpus", 19), ("memory_bytes", 1), ("timeout_s", 60),
    ("interval_s", 2), ("interval_s", True), ("monitor_host", False), ("host_interval_s", 1)])
def test_collection_drift_rejected_before_work(tmp_path, monkeypatch, field, value):
    kwargs = dict(cpus=20, memory_bytes=96*1024**3, timeout_s=900, interval_s=1., monitor_host=True, host_interval_s=30.)
    kwargs[field] = value
    monkeypatch.setattr(module, "measure", lambda *a, **k: pytest.fail("Unexpected work"))
    with pytest.raises(ValueError, match="frozen"):
        module.measure_native_run(["/usr/bin/true"], tmp_path, 123, **kwargs)


def test_supplement_failure_retained_as_workflow_failure(fixture, monkeypatch):
    run, order, root = fixture

    def fail(*args):
        raise ValueError("invalid supplementary context")

    monkeypatch.setattr(module, "measure", fail)
    result = measure_run(run, order, [], lambda _: ["a.fa"], module.measure_native_run,
                         123, timeout_s=900, runtime_checker=lambda _: [])
    assert result["status"] == "verified_wrapper_failed"
    assert result["error"] == "invalid supplementary context"
    assert result["after"]["native_order"] == ["a.fa"]
