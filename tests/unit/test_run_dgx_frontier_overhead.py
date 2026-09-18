import copy
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import run_dgx_frontier_overhead as module

ROOT = Path(__file__).resolve().parents[2] / "benchmark_tools"
PLAN = ROOT / "results/dgx_frontier_overhead_plan_20260918.json"


def test_transferred_recipe_and_authorization_match_frozen_sources():
    manifest_path = ROOT / "results/dgx_frontier_overhead_recipe_v1_20260918.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest["external_symlinks"] == []
    rows = [row for row in manifest["records"] if row["kind"] == "file"]
    assert len(rows) == 39
    for row in rows:
        path = ROOT / Path(row["path"]).name
        if not path.exists():
            path = ROOT / "results" / path.name
        data = path.read_bytes()
        assert len(data) == row["bytes"]
        assert hashlib.sha256(data).hexdigest() == row["sha256"]
    auth = json.loads((ROOT / "results/dgx_frontier_overhead_authorization_20260918.json").read_text())
    assert auth["recipe_sha256"] == hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    assert auth["plan_sha256"] == module.PLAN_SHA
    assert auth["scientific_execution_authorized"] is False
    assert auth["allowed_indices"] == list(range(18))


def save(path, value):
    path.write_text(json.dumps(value))
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def pinned(tmp_path):
    paths = [Path(module.__file__).resolve(), Path(module.periodic.__code__.co_filename).resolve(),
             Path(module.boundary.__code__.co_filename).resolve(), PLAN]
    manifest = {"records": [dict(kind="file", path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths]}
    recipe = tmp_path / "recipe.json"
    recipe_sha = save(recipe, manifest)
    auth = dict(purpose="native_frontier_incremental_overhead", execution_authorized=True,
                scientific_execution_authorized=False, plan_sha256=module.PLAN_SHA,
                recipe_sha256=recipe_sha, allowed_indices=list(range(18)))
    authorization = tmp_path / "auth.json"
    auth_sha = save(authorization, auth)
    return [PLAN, authorization, auth_sha, recipe, recipe_sha]


@pytest.mark.parametrize("index", range(18))
def test_every_frozen_task_selects_with_exact_authorization(pinned, index):
    plan, auth, row = module.select(*pinned, index)
    assert row == plan["runs"][index]
    assert row["index"] == index
    assert auth["scientific_execution_authorized"] is False


@pytest.mark.parametrize("index", [-1, 18, True, 1.5, "1"])
def test_invalid_index_rejected(pinned, index):
    with pytest.raises(ValueError, match="integer index"):
        module.select(*pinned, index)


@pytest.mark.parametrize("key,value", [("purpose", "scientific_scaling"),
    ("scientific_execution_authorized", True), ("execution_authorized", False),
    ("recipe_sha256", "0"*64), ("plan_sha256", "0"*64), ("allowed_indices", [0]),
    ("execution_authorized", 1), ("scientific_execution_authorized", 0),
    ("allowed_indices", [False, *range(1, 18)])])
def test_scope_or_recipe_mismatch_rejected(pinned, key, value):
    auth = json.loads(pinned[1].read_text())
    auth[key] = value
    pinned[2] = save(pinned[1], auth)
    with pytest.raises(ValueError, match="scope/recipe"):
        module.select(*pinned, 0)


def test_missing_collector_in_recipe_rejected(pinned):
    manifest = json.loads(pinned[3].read_text())
    manifest["records"].pop(1)
    pinned[4] = save(pinned[3], manifest)
    auth = json.loads(pinned[1].read_text())
    auth["recipe_sha256"] = pinned[4]
    pinned[2] = save(pinned[1], auth)
    with pytest.raises(ValueError, match="missing from pinned recipe"):
        module.select(*pinned, 0)


@pytest.mark.parametrize("index", [0, 1])
@pytest.mark.parametrize("status", ["command_exited_zero", "command_failed"])
def test_launch_uses_correct_collector_and_retains_failure(pinned, tmp_path, monkeypatch, index, status):
    plan, auth, row = module.select(*pinned, index)
    row = copy.deepcopy(row)
    directory = tmp_path / "run"
    row["run"]["measurement_directory"] = str(directory / "measurement")
    monkeypatch.setattr(module, "select", lambda *args: (plan, auth, row))
    monkeypatch.setattr(module.os, "uname", lambda: SimpleNamespace(nodename="spark-7ff0"))
    monkeypatch.setattr(module.os, "chdir", lambda path: None)
    monkeypatch.setattr(module.sys, "dont_write_bytecode", True)
    for key, value in dict(SLURM_JOB_ID="123", SLURM_CPUS_PER_TASK="20", SLURM_MEM_PER_NODE="98304", PYTHONHASHSEED="0").items():
        monkeypatch.setenv(key, value)
    for key in ("LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"):
        monkeypatch.delenv(key, raising=False)
    original_exists = Path.exists
    monkeypatch.setattr(Path, "exists", lambda path: False if str(path) == "/etc/ld.so.preload"
                        or str(path).startswith(str(module.ROOT / "frontier_overhead_v1")) else original_exists(path))
    calls = []
    def measure(*args, **kwargs):
        calls.append((args, kwargs))
        directory.mkdir()
        return {"status": status}
    monkeypatch.setattr(module, "measure_run", measure)
    assert module.launch(*pinned, index) == {"status": status}
    args, kwargs = calls[0]
    assert args[0] == row["run"]
    assert args[1] == plan["order"]
    assert args[2][-1] == (str(pinned[3]), pinned[4])
    assert args[4] is (module.boundary if row["mode"] == "boundary" else module.periodic)
    assert args[5] == 123
    assert kwargs == {"timeout_s": 900}
    evidence = json.loads((directory / "overhead_task.json").read_text())
    assert evidence["task"] == row
    assert evidence["authorization"] == auth
    assert evidence["scientific_timings_admitted"] is False
