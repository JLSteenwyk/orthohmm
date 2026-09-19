import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import run_root_context_native as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record

FOLDER = Path(module.__file__).parent


def plan():
    return json.loads((FOLDER / "results/dgx_root_context_native_plan_20260919.json").read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def recipe(tmp_path):
    paths = [*FOLDER.glob("*.py"), FOLDER / "run_dgx_root_context_native.sh",
             FOLDER / "results/dgx_root_context_native_plan_20260919.json",
             FOLDER / "results/dgx_lineage_native_plan_20260919.json",
             FOLDER / "results/ROOT_CONTEXT_NATIVE_PROTOCOL_20260919.md"]
    path = tmp_path / "recipe.json"
    path.write_text(json.dumps(dict(records=[dict(record(p), kind="file") for p in paths])))
    return path


def test_exact_plan_and_recipe_selection(recipe):
    assert module.select(FOLDER / "results/dgx_root_context_native_plan_20260919.json", recipe, sha(recipe)) == plan()


@pytest.mark.parametrize("fault", ["missing", "changed", "duplicate", "digest"])
def test_source_binding_rejects_drift(recipe, fault):
    data = json.loads(recipe.read_text())
    if fault == "missing":
        data["records"].pop()
    elif fault == "changed":
        data["records"][0]["sha256"] = "0"*64
    elif fault == "duplicate":
        data["records"].append(data["records"][0])
    recipe.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        module.select(FOLDER / "results/dgx_root_context_native_plan_20260919.json", recipe,
                      "0"*64 if fault == "digest" else sha(recipe))


@pytest.mark.parametrize("failure", [None, "wrapper", "exception", "native"])
def test_serial_panel_stops_and_retains_unrun(tmp_path, monkeypatch, failure):
    value = plan()
    for task in value["runs"]:
        task["run"]["measurement_directory"] = str(tmp_path / "root_context_native_v1" / f"run_{task['index']:02d}" / "measurement")
    monkeypatch.setattr(module, "ROOT", tmp_path)
    monkeypatch.setattr(module, "select", lambda *args: value)
    monkeypatch.setattr(module, "preflight", lambda *args: 123)
    calls = []

    def execute(plan, task, *args):
        index = task["index"]
        calls.append(index)
        if index == 1 and failure == "exception":
            raise RuntimeError("cleanup failed")
        result = dict(status="command_failed" if index == 1 and failure == "wrapper" else "command_exited_zero",
            measurement=dict(native=dict(exit_code=1 if index == 1 and failure == "native" else 0, timed_out=False),
                             native_wall_s=1.0))
        root = Path(task["run"]["measurement_directory"]).parent
        root.mkdir()
        (root / "verification.json").write_text(json.dumps(result))
        return result

    monkeypatch.setattr(module, "execute_task", execute)
    result = module.panel(tmp_path / "plan.json", tmp_path / "recipe.json", "sha")
    assert calls == ([0, 1, 2] if failure is None else [0, 1])
    assert len(result["tasks"]) == 3
    assert result["native_outputs_validated"] is False and result["scientific_timings_admitted"] is False
    if failure is None:
        assert result["status"] == "root_context_native_completed"
    else:
        assert result["status"] == "root_context_native_stopped"
        assert result["tasks"][1]["status"] == "failed"
        assert result["tasks"][2]["status"] == "not_run_after_failure"
        assert not (tmp_path / "root_context_native_v1/run_02").exists()
    checkpoint = json.loads((tmp_path / "root_context_native_v1/progress_02.json").read_text())
    assert checkpoint["tasks"] == result["tasks"]


@pytest.mark.parametrize("fails", [False, True])
def test_task_environment_and_cwd_are_restored(tmp_path, monkeypatch, fails):
    value = plan()
    task = value["runs"][0]
    task["run"]["cwd"] = str(tmp_path)
    monkeypatch.setattr(module, "ROOT", tmp_path)
    monkeypatch.setenv("PYTHONPATH", "original")
    monkeypatch.setenv("OMP_NUM_THREADS", "original_threads")
    monkeypatch.setenv("CONDA_PREFIX", "original_conda")
    original_cwd, original_env = Path.cwd(), dict(module.os.environ)
    calls = []

    def measure(run, order, runtime, enumerator, collector, job, **kwargs):
        calls.append(run)
        assert Path.cwd() == tmp_path
        assert module.os.environ["PYTHONPATH"] == value["environment_overrides"]["PYTHONPATH"]
        assert "CONDA_PREFIX" not in module.os.environ
        assert collector is module.measure_native_run and kwargs == dict(timeout_s=900)
        if fails:
            raise RuntimeError("native failure")
        return dict(status="command_exited_zero")

    monkeypatch.setattr(module, "measure_run", measure)
    if fails:
        with pytest.raises(RuntimeError):
            module.execute_task(value, task, tmp_path / "recipe.json", "sha", 123)
    else:
        module.execute_task(value, task, tmp_path / "recipe.json", "sha", 123)
    assert len(calls) == 1
    assert Path.cwd() == original_cwd and dict(module.os.environ) == original_env


def test_preflight_environment_limits(tmp_path, monkeypatch):
    value = dict(launcher_python=module.sys.executable)
    monkeypatch.setattr(module.os, "uname", lambda: SimpleNamespace(nodename="spark-7ff0"))
    monkeypatch.setattr(module.os, "sched_getaffinity", lambda pid: set(range(20)))
    monkeypatch.setattr(module.sys, "dont_write_bytecode", True)
    monkeypatch.setattr(module.sys, "pycache_prefix", str(tmp_path / "absent"))
    for key, val in dict(SLURM_CPUS_PER_TASK="20", SLURM_MEM_PER_NODE="98304", SLURM_JOB_ID="123",
                         PYTHONHASHSEED="0", PYTHONNOUSERSITE="1").items():
        monkeypatch.setenv(key, val)
    for key in ("LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT", "PYTHONPATH"):
        monkeypatch.delenv(key, raising=False)
    assert module.preflight(value) == 123
    monkeypatch.setenv("PYTHONPATH", "/unexpected")
    with pytest.raises(ValueError):
        module.preflight(value)


def test_existing_panel_not_reused(tmp_path, monkeypatch):
    monkeypatch.setattr(module, "ROOT", tmp_path)
    monkeypatch.setattr(module, "select", lambda *args: plan())
    monkeypatch.setattr(module, "preflight", lambda *args: 123)
    (tmp_path / "root_context_native_v1").mkdir()
    with pytest.raises(FileExistsError):
        module.panel(tmp_path / "plan", tmp_path / "recipe", "sha")
