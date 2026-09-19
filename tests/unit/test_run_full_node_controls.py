import hashlib
from pathlib import Path

import pytest

from benchmark_tools import run_full_node_controls as module


@pytest.mark.parametrize("fail_index,missed", [(None, False), (2, False), (None, True)])
def test_panel_retains_fixed_order_and_failures(tmp_path, monkeypatch, fail_index, missed):
    calls = []

    def trial(directory, mode, job):
        calls.append(mode)
        if len(calls)-1 == fail_index:
            raise RuntimeError("injected trial failure")
        return dict(mode=mode, status="workload_validated", job_id=job,
                    positive_control_detected=not missed if mode == "contended" else None)

    monkeypatch.setattr(module, "trial", trial)
    result = module.panel(tmp_path / "panel", 17)
    assert calls == [mode for modes in module.ORDER for mode in modes]
    assert len(result["trials"]) == 9
    assert len(list((tmp_path / "panel").glob("trial_*/panel_trial.json"))) == 9
    assert result["summary"]["all_workloads_valid"] == (fail_index is None)
    assert result["summary"]["all_positive_controls_detected"] == (fail_index is None and not missed)
    assert result["scientific_timings_admitted"] is False
    if fail_index is not None:
        assert result["trials"][fail_index]["error"] == "injected trial failure"


def test_summary_rejects_incomplete_inventory():
    with pytest.raises(ValueError, match="inventory/order"):
        module.summarize([])


def test_existing_panel_never_reused(tmp_path):
    with pytest.raises(FileExistsError):
        module.panel(tmp_path, 17)


@pytest.fixture
def launch_context(tmp_path, monkeypatch):
    base = Path(module.__file__).resolve().parent
    paths = [*base.glob("*.py"), base / "results/FULL_NODE_CPU_CONTROL_PROTOCOL_20260919.md",
             base / "results/dgx_dual_native_plan_20260919.json"]
    recipe = dict(records=[dict(path=str(path), kind="file",
                               sha256=hashlib.sha256(path.read_bytes()).hexdigest()) for path in paths])
    plan = dict(launcher_python=module.sys.executable)
    monkeypatch.setattr(module, "read_pinned", lambda path, sha: plan if sha == module.PLAN_SHA else recipe)
    monkeypatch.setattr(module.platform, "node", lambda: "spark-7ff0")
    monkeypatch.setattr(module.os, "sched_getaffinity", lambda pid: set(range(20)))
    monkeypatch.setattr(module.sys, "dont_write_bytecode", True)
    monkeypatch.setattr(module.sys, "pycache_prefix", str(tmp_path / "absent_cache"))
    for key, value in dict(SLURM_CPUS_PER_TASK="20", SLURM_MEM_PER_NODE="98304",
                           SLURM_JOB_ID="17", PYTHONHASHSEED="0").items():
        monkeypatch.setenv(key, value)
    for key in ("LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT", "PYTHONPATH"):
        monkeypatch.delenv(key, raising=False)
    return recipe


def test_preflight_checks_pinned_environment(tmp_path, launch_context):
    _, job = module.preflight(tmp_path / "recipe.json", "fixture")
    assert job == 17


@pytest.mark.parametrize("fault", ["cpus", "memory", "loader", "pythonpath", "source", "cache", "protocol"])
def test_preflight_rejects_launch_drift(tmp_path, monkeypatch, launch_context, fault):
    if fault == "cpus":
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "1")
    elif fault == "memory":
        monkeypatch.setenv("SLURM_MEM_PER_NODE", "1")
    elif fault == "loader":
        monkeypatch.setenv("LD_PRELOAD", "/unexpected.so")
    elif fault == "pythonpath":
        monkeypatch.setenv("PYTHONPATH", "/unexpected")
    elif fault == "source":
        launch_context["records"][0]["sha256"] = "0"*64
    elif fault == "cache":
        Path(module.sys.pycache_prefix).mkdir()
    else:
        monkeypatch.setattr(module, "PROTOCOL_SHA", "0"*64)
    with pytest.raises(ValueError):
        module.preflight(tmp_path / "recipe.json", "fixture")
