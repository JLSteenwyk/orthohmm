import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import run_dual_overhead_panel as module


@pytest.fixture
def context(tmp_path, monkeypatch):
    base = Path(module.__file__).resolve().parent
    plan_path = base / "results/dgx_dual_overhead_plan_20260919.json"
    plan = module.read_pinned(plan_path, module.PLAN_SHA)
    paths = [*base.glob("*.py"), plan_path, base / "results/DUAL_COLLECTOR_OVERHEAD_PROTOCOL_20260919.md"]
    recipe = dict(records=[dict(kind="file", path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths])
    auth = module.authorization("recipe_sha")
    recipe_path, auth_path = tmp_path / "recipe.json", tmp_path / "auth.json"
    values = {plan_path: plan, recipe_path: recipe, auth_path: auth}
    monkeypatch.setattr(module, "read_pinned", lambda p, sha: values[p])
    return plan_path, recipe_path, "recipe_sha", auth_path, "auth_sha", plan, recipe, auth


@pytest.mark.parametrize("index", [0, 1, 17])
def test_select_actual_complete_plan(context, index):
    result = module.select(*context[:5], index)
    assert result[2]["index"] == index
    assert result[0]["scientific_timings_admitted"] is False


@pytest.mark.parametrize("index", [-1, 18, True, 1.0])
def test_invalid_index(context, index):
    with pytest.raises(ValueError, match="integer task"):
        module.select(*context[:5], index)


@pytest.mark.parametrize("fault", ["scope", "boolean", "source", "protocol", "inventory", "admission"])
def test_rejects_unpinned_or_changed_scope(context, fault):
    plan, recipe, auth = context[5:]
    if fault == "scope":
        auth["allowed_indices"].pop()
    elif fault == "boolean":
        auth["execution_authorized"] = 1
    elif fault == "source":
        recipe["records"][0]["sha256"] = "0"*64
    elif fault == "protocol":
        plan["protocol_sha256"] = "0"*64
    elif fault == "inventory":
        plan["runs"].pop()
    else:
        plan["scientific_timings_admitted"] = True
    with pytest.raises(ValueError):
        module.select(*context[:5], 0)


@pytest.mark.parametrize("index", [0, 1])
def test_launch_selects_intended_arm_and_retains_receipt(context, tmp_path, monkeypatch, index):
    plan = context[5]
    task = plan["runs"][index]
    directory = tmp_path / "native"
    directory.mkdir()
    task["run"]["measurement_directory"] = str(directory / "measurement")
    plan["cache_directory"] = str(tmp_path / "cache")
    env = dict(module.os.environ)
    for key in ("LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT", "PYTHONPATH"):
        env.pop(key, None)
    env.update(SLURM_CPUS_PER_TASK="20", SLURM_MEM_PER_NODE="98304", SLURM_ARRAY_TASK_ID=str(index),
               SLURM_JOB_ID="123", PYTHONHASHSEED="0")
    monkeypatch.setattr(module.os, "environ", env)
    monkeypatch.setattr(module.os, "uname", lambda: SimpleNamespace(nodename="spark-7ff0"))
    monkeypatch.setattr(module.os, "chdir", lambda path: None)
    monkeypatch.setattr(module.sys, "executable", plan["launcher_python"])
    monkeypatch.setattr(module.sys, "dont_write_bytecode", True)
    monkeypatch.setattr(module.sys, "pycache_prefix", str(tmp_path / "cache" / f"cache_{index}"))
    monkeypatch.setattr(module, "native_enumerator", lambda value: "enumerator")
    calls = []

    def measure(*args, **kwargs):
        calls.append((args, kwargs))
        return dict(status="command_exited_zero")

    monkeypatch.setattr(module, "measure_run", measure)
    module.launch(*context[:5], index)
    collector = calls[0][0][4]
    if index == 1:
        assert collector is module.periodic
    else:
        assert collector.func is module.boundary and collector.keywords == {"native_pressure": True}
    assert calls[0][1] == {"timeout_s": 900}
    receipt = json.loads((directory / "overhead_task.json").read_text())
    assert receipt["task"] == task
    assert receipt["scientific_timings_admitted"] is False
