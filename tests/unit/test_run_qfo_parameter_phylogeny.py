import copy
import json
from pathlib import Path

import pytest

from benchmark_tools import run_qfo_parameter_phylogeny as module
from benchmark_tools.prepare_orthobench_factorial import plan_cells


def panel():
    return {
        "status": "corrected_qfo_candidate_neighborhood_admitted_unscored",
        "accuracy_evaluated": False, "publication_ready": False,
        "arms": [{"label": label, "status": "candidate_prepared_unscored",
                  "candidate_arm": {"candidate_partition": {"path": label + ".txt"},
                                    "membership_constraints": {"path": label + ".json"}}}
                 for label in ("control", *module.VARIANTS)]}


@pytest.mark.parametrize("index", range(4))
def test_select_exact_prespecified_arm(index):
    admission = panel()
    before = copy.deepcopy(admission)
    arm = module.select_arm(admission, index)
    assert arm["label"] == module.VARIANTS[index]
    assert arm["partition"] == admission["arms"][index + 1]["candidate_arm"]["candidate_partition"]
    assert admission == before


@pytest.mark.parametrize("problem", ["status", "accuracy", "ready", "order", "missing", "incomplete"])
def test_reject_changed_panel(problem):
    admission = panel()
    if problem == "status":
        admission["status"] = "prepared"
    elif problem == "accuracy":
        admission["accuracy_evaluated"] = True
    elif problem == "ready":
        admission["publication_ready"] = True
    elif problem == "order":
        admission["arms"].reverse()
    elif problem == "missing":
        admission["arms"].pop()
    else:
        admission["arms"][1]["status"] = "failed"
    with pytest.raises(ValueError):
        module.select_arm(admission, 0)


@pytest.mark.parametrize("index", [-1, 4, True, "0", 0.0])
def test_reject_invalid_index(index):
    with pytest.raises(ValueError):
        module.select_arm(panel(), index)


def allocation(monkeypatch):
    for key, value in {"SLURM_JOB_ID": "test", "SLURM_CPUS_PER_TASK": "32",
                       "SLURM_JOB_NODELIST": "bizon", "SLURM_ARRAY_TASK_ID": "0"}.items():
        monkeypatch.setenv(key, value)


@pytest.mark.parametrize("key,value", [("SLURM_JOB_ID", ""), ("SLURM_CPUS_PER_TASK", "2"),
                                      ("SLURM_JOB_NODELIST", "spark-7ff0"), ("SLURM_ARRAY_TASK_ID", "1")])
def test_allocation_gate_precedes_source_access(tmp_path, monkeypatch, key, value):
    allocation(monkeypatch)
    monkeypatch.setenv(key, value)
    monkeypatch.setattr(module, "verify_sources", lambda *a: pytest.fail("Accessed sources"))
    with pytest.raises(ValueError, match="scheduled"):
        module.run(tmp_path, 0)
    assert not (tmp_path / "benchmarks").exists()


@pytest.mark.parametrize("outcome", ["success", "failed", "exception", "source_changed"])
def test_execution_retains_evidence_and_restores_cwd(tmp_path, monkeypatch, outcome):
    allocation(monkeypatch)
    launcher = tmp_path / "launcher"
    launcher.mkdir()
    executor = tmp_path / "executor"
    original = next(c for c in plan_cells(tmp_path / "baseline", tmp_path / "fastas",
                    executor / "benchmark_tools/replay_phylogeny.py", 32) if c["label"] == "p1_c1_r1")
    before = copy.deepcopy(original)
    arm = module.select_arm(panel(), 0)
    manifest = {"input_fastas": [], "environment_overrides": {"OMP_NUM_THREADS": "1"}}
    checks = []
    def verify(*args):
        checks.append(args)
        if outcome == "source_changed" and len(checks) == 2:
            raise ValueError("Source changed")
        return arm, manifest, original, launcher, executor, {}, [], {}
    monkeypatch.setattr(module, "verify_sources", verify)
    monkeypatch.setattr(module, "native_command", lambda cell, *a: (cell["argv"], []))
    monkeypatch.setattr(module, "execution_environment", lambda *a: ({}, {}))
    calls = []
    def execute(dataset, order, env, evidence, inputs, provenance):
        calls.append(dataset)
        assert Path.cwd() == launcher
        assert order == ["candidate_norm_low"]
        assert env == {"OMP_NUM_THREADS": "1", "PYTHONPATH": str(launcher)}
        argv = dataset["methods"][order[0]]["argv"]
        assert argv[argv.index("--species-tree-mode") + 1] == "infer"
        assert "--species-tree" not in argv
        assert argv[-2:] == ["--checkpoint-source", original["argv"][original["argv"].index("--output-directory") + 1]]
        if outcome == "exception":
            raise RuntimeError("Execution interrupted")
        return {"failed_methods": [order[0]] if outcome == "failed" else []}
    monkeypatch.setattr(module, "execute", execute)
    cwd = Path.cwd()
    if outcome == "success":
        module.run(tmp_path, 0)
    else:
        with pytest.raises((RuntimeError, ValueError)):
            module.run(tmp_path, 0)
    assert Path.cwd() == cwd
    assert original == before
    output = tmp_path / "benchmarks/results/qfo_parameter_phylogeny_v1/norm_low"
    report = json.loads((output / "postflight.json").read_text())
    assert report["status"] == ("complete_pending_native_validation" if outcome == "success" else "failed")
    assert report["accuracy_evaluated"] is False
    assert report["native_outputs_validated"] is False
    assert len(checks) == (1 if outcome == "exception" else 2)
    with pytest.raises(FileExistsError):
        module.run(tmp_path, 0)
    assert len(calls) == 1
