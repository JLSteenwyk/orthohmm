import hashlib
import json
from pathlib import Path

import pytest

from benchmark_tools import run_root_context_controls as module

MANAGER = "/user.slice/user-1000.slice/user@1000.service"


@pytest.mark.parametrize("fail_index,missed", [(None, False), (0, False), (3, False), (11, False), (None, True)])
def test_panel_retains_every_condition_and_stops_after_failure(tmp_path, monkeypatch, fail_index, missed):
    calls = []

    def trial(directory, mode, job, index, manager):
        calls.append(index)
        assert manager == MANAGER
        if index == fail_index:
            raise RuntimeError("cleanup not confirmed")
        return dict(status="root_context_workload_validated",
                    positive_control_response=not missed if mode == "user-contended" else None)

    monkeypatch.setattr(module, "trial", trial)
    output = tmp_path / "panel"
    result = module.panel(output, 17, MANAGER)
    assert calls == list(range(12 if fail_index is None else fail_index + 1))
    assert [(row["block"], row["mode"]) for row in result["trials"]] == [
        (block, mode) for block, modes in enumerate(module.ORDER) for mode in modes]
    assert len(list(output.glob("trial_*/panel_trial.json"))) == 12
    assert json.loads((output / "progress_11.json").read_text())["trials"] == result["trials"]
    assert result["summary"]["all_workloads_valid"] == (fail_index is None)
    assert result["summary"]["all_positive_controls_detected"] == (not missed and (fail_index is None or fail_index == 11))
    assert result["scientific_timings_admitted"] is False
    if fail_index is not None:
        assert result["status"] == "root_context_controls_stopped"
        assert result["trials"][fail_index]["error"] == "cleanup not confirmed"
        assert all(row["status"] == "not_run_after_failure" for row in result["trials"][fail_index+1:])
    else:
        assert result["status"] == "root_context_controls_completed"


def test_existing_directory_never_reused(tmp_path):
    with pytest.raises(FileExistsError):
        module.panel(tmp_path, 17, MANAGER)


def test_incomplete_summary_rejected():
    with pytest.raises(ValueError, match="inventory"):
        module.summarize([])


@pytest.fixture
def preflight_context(tmp_path, monkeypatch):
    base = Path(module.__file__).resolve().parent
    paths = [base / "results/ROOT_CPU_CONTEXT_CONTROL_PROTOCOL_20260919.md",
             base / "run_dgx_root_context_controls.sh"]
    recipe = dict(records=[dict(path=str(path), kind="file",
        sha256=hashlib.sha256(path.read_bytes()).hexdigest()) for path in paths])
    calls = []
    monkeypatch.setattr(module, "runtime_preflight", lambda *args: calls.append(args) or ({}, 17))
    monkeypatch.setattr(module, "read_pinned", lambda *args: recipe)
    monkeypatch.setenv("PYTHONNOUSERSITE", "1")
    monkeypatch.setattr(module.subprocess, "check_output", lambda *args, **kwargs: MANAGER + "\n")
    original = Path.is_dir
    monkeypatch.setattr(Path, "is_dir", lambda path: True if str(path) == "/sys/fs/cgroup" + MANAGER else original(path))
    return recipe, calls


def test_preflight_delegates_runtime_and_checks_root_bindings(tmp_path, preflight_context):
    recipe, calls = preflight_context
    assert module.preflight(tmp_path / "recipe.json", "sha") == ({}, 17, MANAGER)
    assert calls == [(tmp_path / "recipe.json", "sha")]


@pytest.mark.parametrize("fault", ["protocol", "script", "hash", "user_site", "manager", "missing_scope"])
def test_preflight_rejects_drift(tmp_path, monkeypatch, preflight_context, fault):
    recipe, _ = preflight_context
    if fault in ("protocol", "script"):
        recipe["records"][fault == "script"]["sha256"] = "0" * 64
    elif fault == "hash":
        monkeypatch.setattr(module, "PROTOCOL_SHA", "0" * 64)
    elif fault == "user_site":
        monkeypatch.delenv("PYTHONNOUSERSITE")
    elif fault == "manager":
        monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: "/system.slice\n")
    else:
        monkeypatch.setattr(Path, "is_dir", lambda path: False)
    with pytest.raises(ValueError):
        module.preflight(tmp_path / "recipe.json", "sha")


def test_shell_has_frozen_resource_limits():
    source = Path(module.__file__).with_name("run_dgx_root_context_controls.sh").read_text()
    for option in ("--cpus-per-task=20", "--mem=96G", "--exclusive", "--time=00:15:00", "--no-requeue"):
        assert "#SBATCH " + option in source
    assert "-m benchmark_tools.run_root_context_controls" in source
