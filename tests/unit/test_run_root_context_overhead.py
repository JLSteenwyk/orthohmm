import hashlib
import json
from pathlib import Path

import pytest

from benchmark_tools import run_root_context_overhead as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record

FOLDER = Path(module.__file__).parent
PLAN = FOLDER / "results/dgx_root_context_overhead_plan_20260919.json"


def plan():
    return json.loads(PLAN.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def recipe(tmp_path):
    paths = [*FOLDER.glob("*.py"), PLAN, FOLDER / "run_dgx_root_context_overhead.sh",
        FOLDER / "results/dgx_root_context_native_plan_20260919.json",
        FOLDER / "results/ROOT_CONTEXT_OVERHEAD_PROTOCOL_20260919.md"]
    path = tmp_path / "recipe.json"
    path.write_text(json.dumps(dict(records=[dict(record(p), kind="file") for p in paths])))
    return path


def test_pinned_plan_and_source_select(recipe):
    assert module.select(PLAN, recipe, sha(recipe)) == plan()


@pytest.mark.parametrize("fault", ["missing", "changed", "duplicate", "digest"])
def test_source_drift_rejected(recipe, fault):
    value = json.loads(recipe.read_text())
    if fault == "missing":
        value["records"].pop()
    elif fault == "changed":
        value["records"][0]["sha256"] = "0"*64
    elif fault == "duplicate":
        value["records"].append(value["records"][0])
    recipe.write_text(json.dumps(value))
    with pytest.raises(ValueError):
        module.select(PLAN, recipe, "0"*64 if fault == "digest" else sha(recipe))


def install_panel(tmp_path, monkeypatch, failure_index=None, fault="exception"):
    output = tmp_path / "panel"
    value = plan()
    for task in value["runs"]:
        task["run"]["measurement_directory"] = str(output / f"run_{task['index']:02d}" / "measurement")
    monkeypatch.setattr(module, "OUTPUT_ROOT", output)
    monkeypatch.setattr(module, "select", lambda *args: value)
    monkeypatch.setattr(module, "preflight", lambda *args: 123)
    calls = []

    def execute(plan, task, *args):
        index = task["index"]
        calls.append(index)
        if index == failure_index and fault == "exception":
            raise RuntimeError("cleanup failed")
        native = dict(exit_code=0, timed_out=False)
        status = "command_exited_zero"
        if index == failure_index:
            if fault == "wrapper":
                status = "command_failed"
            elif fault == "native":
                native["exit_code"] = 1
            elif fault == "bool_exit":
                native["exit_code"] = False
            elif fault == "timeout":
                native["timed_out"] = True
        result = dict(status=status, measurement=dict(native=native, native_wall_s=1.))
        directory = Path(task["run"]["measurement_directory"])
        directory.mkdir(parents=True)
        (directory.parent / "verification.json").write_text(json.dumps(result))
        (directory / "lineage_report.json").write_text("{}")
        if (task["arm"] == "root_context" and not (index == failure_index and fault == "missing_root")) or (
                index == failure_index and fault == "extra_root"):
            (directory / "root_context_report.json").write_text("{}")
        return result

    monkeypatch.setattr(module, "execute_task", execute)
    return value, calls, output


@pytest.mark.parametrize("failure_index", [None, 0, 7, 17])
def test_serial_panel_retains_every_outcome(tmp_path, monkeypatch, failure_index):
    value, calls, output = install_panel(tmp_path, monkeypatch, failure_index)
    result = module.panel(PLAN, tmp_path / "recipe", "sha")
    assert calls == list(range(18 if failure_index is None else failure_index+1))
    assert len(result["tasks"]) == 18
    assert result["native_outputs_validated"] is False and result["scientific_timings_admitted"] is False
    assert result["status"] == ("root_context_overhead_completed" if failure_index is None else "root_context_overhead_stopped")
    for index, row in enumerate(result["tasks"]):
        assert row["task"] == value["runs"][index]
        assert json.loads((output / f"task_{index:02d}.json").read_text()) == row
        checkpoint = json.loads((output / f"progress_{index:02d}.json").read_text())
        assert checkpoint["tasks"] == result["tasks"][:index+1]
        if failure_index is not None and index > failure_index:
            assert row["status"] == "not_run_after_failure" and row["failed_index"] == failure_index
            assert not (output / f"run_{index:02d}").exists()


@pytest.mark.parametrize("fault,index", [("wrapper", 0), ("native", 0), ("bool_exit", 0),
    ("timeout", 0), ("missing_root", 1), ("extra_root", 0)])
def test_failed_execution_or_wrong_arm_artifacts_stop_panel(tmp_path, monkeypatch, fault, index):
    _, calls, _ = install_panel(tmp_path, monkeypatch, index, fault)
    result = module.panel(PLAN, tmp_path / "recipe", "sha")
    assert calls == list(range(index+1))
    assert result["tasks"][index]["status"] == "failed"
    assert all(row["status"] == "not_run_after_failure" for row in result["tasks"][index+1:])


@pytest.mark.parametrize("index", [0, 1, 4, 6])
@pytest.mark.parametrize("fails", [False, True])
def test_exact_collector_and_environment_restoration(tmp_path, monkeypatch, index, fails):
    value = plan()
    task = value["runs"][index]
    task["run"]["cwd"] = str(tmp_path)
    monkeypatch.setattr(module, "OUTPUT_ROOT", tmp_path / "panel")
    monkeypatch.setenv("PYTHONPATH", "old_path")
    monkeypatch.setenv("CONDA_PREFIX", "old_conda")
    before, cwd = dict(module.os.environ), Path.cwd()
    called = []

    def measure(run, order, specifications, enumerate_native, collector, job, **kwargs):
        called.append(collector)
        assert collector is (module.measure_lineage if task["arm"] == "lineage" else module.measure_root)
        assert job == 123 and kwargs == dict(timeout_s=900)
        assert order == value["order"]
        assert specifications[-1] == (str(tmp_path / "recipe"), "sha")
        assert Path.cwd() == tmp_path
        assert "CONDA_PREFIX" not in module.os.environ
        assert module.os.environ["PYTHONPYCACHEPREFIX"] == str(tmp_path / "panel" / f"cache_{index:02d}")
        assert module.os.environ["PATH"] == module.os.pathsep.join(value["environment_paths"][run["environment_role"]])
        if fails:
            raise RuntimeError("measurement failed")
        return {"status": "test"}

    monkeypatch.setattr(module, "measure_run", measure)
    if fails:
        with pytest.raises(RuntimeError):
            module.execute_task(value, task, tmp_path / "recipe", "sha", 123)
    else:
        assert module.execute_task(value, task, tmp_path / "recipe", "sha", 123) == {"status": "test"}
    assert len(called) == 1 and Path.cwd() == cwd and dict(module.os.environ) == before


def test_source_failure_after_first_task_retains_unrun(tmp_path, monkeypatch):
    value, calls, _ = install_panel(tmp_path, monkeypatch)
    selections = []

    def select(*args):
        selections.append(1)
        if len(selections) == 3:
            raise ValueError("source changed after task")
        return value

    monkeypatch.setattr(module, "select", select)
    result = module.panel(PLAN, tmp_path / "recipe", "sha")
    assert calls == [0]
    assert result["tasks"][0]["status"] == "failed"
    assert all(row["status"] == "not_run_after_failure" for row in result["tasks"][1:])


def test_existing_panel_and_dangling_cache_refused(tmp_path, monkeypatch):
    value, calls, output = install_panel(tmp_path, monkeypatch)
    output.mkdir()
    with pytest.raises(FileExistsError):
        module.panel(PLAN, tmp_path / "recipe", "sha")
    assert calls == []
    monkeypatch.undo()
    monkeypatch.setattr(module, "OUTPUT_ROOT", output)
    (output / "cache_00").symlink_to(tmp_path / "missing")
    with pytest.raises(ValueError, match="cache"):
        module.execute_task(value, value["runs"][0], tmp_path / "recipe", "sha", 123)
