from copy import deepcopy
import os
from pathlib import Path

import pytest

from benchmark_tools import measure_root_context_scaling as module

PLAN = Path(__file__).resolve().parents[2] / "benchmark_tools/results/dgx_root_context_scaling_plan_20260920.json"


@pytest.mark.parametrize("index", range(27))
def test_every_task_uses_its_exact_native_dataset_order(index):
    plan, task, order = module.load_task(PLAN, index)
    assert task == plan["runs"][index]
    assert order["proteomes"] == task["proteomes"]
    assert order["input_directory"] == task["run"]["dataset"]["input_directory"]
    assert sorted(order["inputs_in_native_order"], key=lambda r:r["path"]) == sorted(
        task["run"]["dataset"]["inputs"], key=lambda r:r["path"])


@pytest.mark.parametrize("index", [-1, 27, True, False, "0", 0.0])
def test_invalid_indices_rejected_before_read(index):
    with pytest.raises(ValueError, match="index"):
        module.load_task(Path("missing"), index)


def test_modified_plan_rejected(tmp_path):
    path = tmp_path / "plan.json"
    path.write_text("{}")
    with pytest.raises(ValueError):
        module.load_task(path, 0)


def configured(tmp_path, monkeypatch, index=0):
    plan, task, order = deepcopy(module.load_task(PLAN, index))
    native_cwd = tmp_path / "native"
    native_cwd.mkdir()
    task["run"]["cwd"] = str(native_cwd)
    monkeypatch.setattr(module, "load_task", lambda *a: (plan, task, order))
    monkeypatch.setattr(module, "OUTPUT_ROOT", tmp_path / "outputs")
    monkeypatch.setattr(module, "native_enumerator", lambda spec: ("enumerator", spec))
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("CONDA_PREFIX", "retain-me")
    monkeypatch.setenv("PATH", "/original-path")
    monkeypatch.delenv("PYTHONDONTWRITEBYTECODE", raising=False)
    return plan, task, order


@pytest.mark.parametrize("index", range(27))
def test_all_tasks_forward_exact_work_and_restore_environment(tmp_path, monkeypatch, index):
    plan, task, order = configured(tmp_path, monkeypatch, index)
    before = dict(os.environ)
    expected = {"status": "retained_native_status", "exit_code": 7}
    calls = []

    def measured(run, actual_order, specs, enumerator, collector, job, **kwargs):
        calls.append(index)
        assert run == task["run"] and actual_order == order
        assert specs == [(r["path"], r["sha256"]) for r in plan["runtime_manifests"]] + [
            (str(tmp_path / "recipe.json"), "a"*64)]
        assert enumerator == ("enumerator", plan["enumerator"])
        assert collector is module.measure_native_run and job == 123
        assert kwargs == dict(cpus=20, memory_gib=96, timeout_s=85800)
        assert Path.cwd() == tmp_path / "native"
        assert "CONDA_PREFIX" not in os.environ
        assert os.environ["PYTHONDONTWRITEBYTECODE"] == "1"
        assert os.environ["PYTHONPYCACHEPREFIX"] == str(tmp_path / "outputs" / f"cache_{index:02d}")
        assert os.environ["PATH"] == os.pathsep.join(plan["environment_paths"][run["environment_role"]])
        return expected

    monkeypatch.setattr(module, "measure_run", measured)
    assert module.measure_task(PLAN, index, Path("recipe.json"), "a"*64, 123) is expected
    assert calls == [index]
    assert dict(os.environ) == before and Path.cwd() == tmp_path


@pytest.mark.parametrize("error", [RuntimeError, KeyboardInterrupt])
def test_exception_and_interrupt_restore_environment(tmp_path, monkeypatch, error):
    configured(tmp_path, monkeypatch)
    before = dict(os.environ)
    def fail(*a, **k):
        raise error("retained failure")
    monkeypatch.setattr(module, "measure_run", fail)
    with pytest.raises(error):
        module.measure_task(PLAN, 0, Path("recipe.json"), "a"*64, 123)
    assert dict(os.environ) == before and Path.cwd() == tmp_path


@pytest.mark.parametrize("kind", ["file", "directory", "dangling_link"])
def test_existing_cache_rejected_without_environment_changes(tmp_path, monkeypatch, kind):
    configured(tmp_path, monkeypatch)
    prefix = tmp_path / "outputs/cache_00"
    prefix.parent.mkdir()
    if kind == "file": prefix.write_text("existing")
    elif kind == "directory": prefix.mkdir()
    else: prefix.symlink_to(tmp_path / "absent")
    before = dict(os.environ)
    monkeypatch.setattr(module, "measure_run", lambda *a, **k: pytest.fail("unexpected measurement"))
    with pytest.raises(ValueError, match="cache"):
        module.measure_task(PLAN, 0, Path("recipe.json"), "a"*64, 123)
    assert dict(os.environ) == before and Path.cwd() == tmp_path


@pytest.mark.parametrize("job,sha", [(0, "a"*64), (-1, "a"*64), (True, "a"*64), (1, "bad"), (1, None)])
def test_bad_job_or_recipe_rejected_before_loading(job, sha, monkeypatch):
    monkeypatch.setattr(module, "load_task", lambda *a: pytest.fail("premature plan load"))
    with pytest.raises(ValueError):
        module.measure_task(PLAN, 0, Path("recipe.json"), sha, job)
