import json
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

from benchmark_tools.run_qfo_corrected_replay import run, validate_completion
from benchmark_tools.prepare_qfo_corrected_replay import command_for
from benchmark_tools.prepare_ob_candidate_neighborhood import record

STAGES = ["multipass", "multipass_refined", "strict_profiles", "strict_profiles_refined"]


def fixture():
    worker = {"status": "corrected_checked_replay_worker_returned", "calls": [
        {"stage": stage, "status": "checked", "exit_code": 0}
        for stage in ("initial", "multipass", "profile_base", "profile_expanded")]}
    replay = {"counts": {"genes": 984137, "profiles_built": 10}, "stages": [{"label": s} for s in STAGES]}
    return worker, replay


def test_complete_corrected_replay():
    validate_completion(*fixture(), STAGES)


@pytest.mark.parametrize("problem", ["worker_failed", "missing_call", "call_failed", "wrong_stage",
    "call_exit", "old_genes", "zero_profiles", "nan_profiles", "boolean_profiles", "missing_stage", "reordered_stage"])
def test_reject_incomplete_replay(problem):
    worker, replay = fixture()
    if problem == "worker_failed":
        worker["status"] = "failed"
    elif problem == "missing_call":
        worker["calls"].pop()
    elif problem == "call_failed":
        worker["calls"][0]["status"] = "failed"
    elif problem == "wrong_stage":
        worker["calls"][0]["stage"] = "profile_base"
    elif problem == "call_exit":
        worker["calls"][0]["exit_code"] = 1
    elif problem == "old_genes":
        replay["counts"]["genes"] = 976504
    elif problem == "zero_profiles":
        replay["counts"]["profiles_built"] = 0
    elif problem == "nan_profiles":
        replay["counts"]["profiles_built"] = float("nan")
    elif problem == "boolean_profiles":
        replay["counts"]["profiles_built"] = True
    elif problem == "missing_stage":
        replay["stages"].pop()
    else:
        replay["stages"].reverse()
    with pytest.raises(ValueError):
        validate_completion(worker, replay, STAGES)


def test_requires_scheduled_allocation(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError, match="scheduled"):
        run(tmp_path, tmp_path / "missing.json", "not-read")
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("outcome", ["equal", "different", "missing_gene", "process_failed"])
def test_parent_preserves_outcomes_without_score_transfer(tmp_path, monkeypatch, outcome):
    # Parent orchestration fixture, not a biological or frozen-runtime replay.
    monkeypatch.setenv("SLURM_JOB_ID", "fixture")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "32")
    output = tmp_path / "output"
    launcher = tmp_path / "benchmarks/work/publication_qfo_replay_native_v1"
    names = tmp_path / "gene_names.txt"
    names.write_text("a\nb\nc\n")
    native = tmp_path / "native.txt"
    native.write_text("OG0: a b\nOG1: c\n")
    admission = tmp_path / "admission.json"
    admission.write_text(json.dumps({"content": {"native_groups": record(native)}}))
    checkpoint = tmp_path / "checkpoint/manifest.json"
    plan = {"output_root": str(output), "input_fastas": [{"path": str(tmp_path / "inputs/s.fa")}],
        "checkpoint_manifest": {"path": str(checkpoint), "sha256": "fixture"}, "runtime": {},
        "native_command": command_for(Path(sys.executable), launcher, output, tmp_path / "inputs", checkpoint.parent, "fixture"),
        "cwd": str(launcher), "expected_stages": STAGES,
        "environment_overrides": {"PYTHONPATH": str(launcher), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                                  "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}}
    plan_path = tmp_path / "plan.json"
    plan_path.write_text(json.dumps(plan))
    monkeypatch.setattr("benchmark_tools.checked_replay_payload_worker.corrected_evidence",
        lambda *args: (plan, record(plan_path), record(admission), record(names)))
    monkeypatch.setattr("benchmark_tools.verify_qfo_replay_launcher.verify", lambda *args: {})
    calls = []
    def subprocess_run(command, **kwargs):
        calls.append(command)
        if outcome == "process_failed":
            return SimpleNamespace(returncode=1)
        worker, replay = fixture()
        text = "a b\nc\n" if outcome == "equal" else "a\nb c\n"
        if outcome == "missing_gene":
            text = "a b\n"
        for i, stage in enumerate(replay["stages"]):
            path = output / f"stage{i}.txt"
            path.write_text(text)
            stage["output"] = record(path)
        (output / "checked_worker.json").write_text(json.dumps(worker))
        (output / "replay.json").write_text(json.dumps(replay))
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr("benchmark_tools.run_qfo_corrected_replay.subprocess.run", subprocess_run)
    monkeypatch.setattr("benchmark_tools.run_qfo_corrected_replay.subprocess.check_output", lambda *a, **k: "fixture-commit\n")
    if outcome in ("process_failed", "missing_gene"):
        with pytest.raises((ValueError, RuntimeError)):
            run(tmp_path, plan_path, record(plan_path)["sha256"])
        result = json.loads((output / "results.json").read_text())
        assert result["status"] == "failed"
    else:
        result = run(tmp_path, plan_path, record(plan_path)["sha256"])
        assert result["status"] == "corrected_checked_replay_complete_pending_admission"
        assert result["native_partition_comparison"]["partition_equal"] is (outcome == "equal")
    assert result["accuracy_evaluated"] is False
    assert len(calls) == 1
