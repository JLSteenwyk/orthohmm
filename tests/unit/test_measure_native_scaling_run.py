import json
from pathlib import Path

import pytest

from benchmark_tools.measure_native_scaling_run import measure_run
from benchmark_tools.snapshot_orthohmm_input_order import record


@pytest.fixture
def fixture(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    source = inputs / "a.fa"
    source.write_text(">a\nMALW\n")
    root = tmp_path / "run"
    run = {"measurement_directory": str(root / "measurement"), "cwd": str(tmp_path),
           "native_method": "orthohmm_high_sensitivity", "native_argv": ["/python", "-m", "orthohmm"],
           "configuration": {"output": str(root / "native"), "metrics": str(root / "metrics.json")},
           "dataset": {"input_directory": str(inputs), "proteomes": 1, "inputs": [record(source)]}}
    order = {"input_directory": str(inputs), "proteomes": 1, "native_order": ["a.fa"],
             "inputs_in_native_order": [record(source)]}
    return run, order, root


def invoke(run, order, collector, enumerator=lambda _: ["a.fa"], checker=lambda _: []):
    return measure_run(run, order, [], enumerator, collector, 123, runtime_checker=checker)


def test_preparation_precedes_exact_collector_call(fixture):
    run, order, root = fixture
    events = []
    def check(specs):
        events.append("runtime")
        return []
    def collect(command, directory, job, cpus, memory, timeout, interval, **kwargs):
        events.append("collect")
        assert (root / "preparation.json").exists()
        assert (root / "native").is_dir()
        assert not directory.exists()
        assert command[-3:] == run["native_argv"]
        assert (job, cpus, memory, timeout, interval) == (123, 20, 96 * 1024 ** 3, 86400, 1.)
        return {"status": "command_exited_zero"}
    result = invoke(run, order, collect, checker=check)
    assert events == ["runtime", "collect", "runtime"]
    assert result["status"] == "command_exited_zero"
    assert not result["scientific_results_admitted"]


def test_changed_enumeration_prevents_launch(fixture):
    run, order, root = fixture
    result = invoke(run, order, lambda *a, **k: pytest.fail("launched"), enumerator=lambda _: [])
    assert result["status"] == "verified_wrapper_failed"
    assert (root / "verification.json").exists()
    assert not (root / "native").exists()


def test_native_failure_retained_and_inputs_rechecked(fixture):
    run, order, root = fixture
    result = invoke(run, order, lambda *a, **k: {"status": "command_failed", "exit_code": 7})
    assert result["status"] == "command_failed"
    assert result["after"]["native_order"] == ["a.fa"]
    assert json.loads((root / "verification.json").read_text())["measurement"]["exit_code"] == 7


def test_input_mutation_after_success_invalidates_run(fixture):
    run, order, root = fixture
    def collect(*args, **kwargs):
        Path(run["dataset"]["inputs"][0]["path"]).write_text("changed")
        return {"status": "command_exited_zero"}
    result = invoke(run, order, collect)
    assert result["status"] == "runtime_changed_or_unverifiable"


def test_orthofinder_copies_checked_after_inference(fixture):
    run, order, root = fixture
    run["native_method"] = "orthofinder_full"
    run["configuration"].pop("metrics")
    run["configuration"].update(copy_inputs_from=run["dataset"]["input_directory"],
                                copy_inputs_to=str(root / "native/input"))
    result = invoke(run, order, lambda *a, **k: {"status": "command_exited_zero"})
    assert len(result["before"]["original_inputs"]) == 1
    assert len(result["after"]["copied_inputs"]) == 1


def test_wrong_cwd_does_not_create_run(fixture):
    run, order, root = fixture
    run["cwd"] = "/"
    with pytest.raises(ValueError, match="working directory"):
        invoke(run, order, lambda *a, **k: None)
    assert not root.exists()
