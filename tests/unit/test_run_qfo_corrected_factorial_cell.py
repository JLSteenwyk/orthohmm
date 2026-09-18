import json
from pathlib import Path

import pytest

from benchmark_tools import run_qfo_corrected_factorial_cell as module
from benchmark_tools.admit_qfo_corrected_candidates import PARAMETERS
from benchmark_tools.prepare_orthobench_factorial import plan_cells
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def fixture(root):
    executor = root / "benchmarks/work/publication_qfo_corrected_candidates_v1"
    output = root / "benchmarks/results/qfo_corrected_factorial_v1"
    fasta = root / "inputs"
    arms = {}
    for label in ("p0_c0", "p0_c1", "p1_c0", "p1_c1"):
        expanded = label.endswith("c1")
        arm = {"candidate_expansion": expanded, "content_audit": {"label": label},
               "candidate_partition": {"path": str(output / "candidates" / label / "orthohmm_working_res/orthohmm_edges_clustered.txt")}}
        if expanded:
            arm.update(membership_constraints={}, expansion={"parameters": dict(PARAMETERS)})
        arms[label] = arm
    manifest = {"status": "corrected_qfo_four_candidate_arms_prepared_unscored", "accuracy_computed": False,
        "candidate_arms": arms, "job_id": "20", "input_fastas": [{"path": str(fasta / "s.fasta")}],
        "cells": plan_cells(output, fasta, executor / "benchmark_tools/replay_phylogeny.py", 32),
        "core_root": str(root / "benchmarks/work/publication_method_native_v2"),
        "launcher_root": str(root / "benchmarks/work/publication_qfo_replay_native_v1"),
        "environment_overrides": {"PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}}
    admission = {"status": "corrected_qfo_candidates_admitted", "accuracy_evaluated": False, "publication_ready": False,
        "scheduler": {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "2", "JobIDRaw": "20"},
        "cells": manifest["cells"], "candidate_arms": {k: v["content_audit"] for k, v in arms.items()}}
    return admission, manifest, output, executor


@pytest.mark.parametrize("index", [0, 1, 2, 3])
def test_four_corrected_reconciliation_cells(tmp_path, index):
    admission, manifest, output, executor = fixture(tmp_path)
    cell, observed, source = module.select_cell(admission, manifest, tmp_path, index)
    assert cell["label"] == ("p0_c0_r1", "p0_c1_r1", "p1_c0_r1", "p1_c1_r1")[index]
    assert observed == output and source == executor
    assert ("--membership-constraints" in cell["argv"]) == cell["candidate_expansion"]


@pytest.mark.parametrize("problem", ["status", "accuracy", "missing_arm", "audit", "cells", "index", "bool_index", "fastas"])
def test_reject_wrong_candidate_binding(tmp_path, problem):
    admission, manifest, *_ = fixture(tmp_path)
    index = 0
    if problem == "status":
        admission["status"] = "prepared"
    elif problem == "accuracy":
        admission["accuracy_evaluated"] = True
    elif problem == "missing_arm":
        admission["candidate_arms"].pop("p1_c1")
    elif problem == "audit":
        admission["candidate_arms"]["p0_c0"] = {}
    elif problem == "cells":
        admission["cells"] = []
    elif problem == "index":
        index = 4
    elif problem == "bool_index":
        index = False
    else:
        manifest["input_fastas"].append({"path": "/foreign/input.fasta"})
    with pytest.raises(ValueError):
        module.select_cell(admission, manifest, tmp_path, index)


@pytest.mark.parametrize("outcome", ["success", "failed", "preflight"])
def test_execution_handoff_preserves_failure_and_cwd(tmp_path, monkeypatch, outcome):
    admission, manifest, output, executor = fixture(tmp_path)
    launcher = Path(manifest["launcher_root"])
    launcher.mkdir(parents=True)
    admission_path, env_path = tmp_path / "admission.json", tmp_path / "environment.json"
    admission_path.write_text("{}")
    env_path.write_text("{}")
    prepared = tmp_path / "manifest.json"
    prepared.write_text("{}")
    admission["prepared_manifest"] = record(prepared)
    manifest["input_fastas"] = []
    cell = manifest["cells"][1]
    checks, calls = [], []
    def verify(*args):
        checks.append(args)
        return admission, manifest, cell, output, executor, {"JobIDRaw": "fixture"}
    monkeypatch.setattr(module, "verify_admission", verify)
    monkeypatch.setattr(module, "native_command", lambda *a: (cell["argv"], []))
    monkeypatch.setattr(module, "read_frozen", lambda *a: {})
    monkeypatch.setattr(module, "verify_environment", lambda *a: None)
    monkeypatch.setattr(module, "execution_environment", lambda *a: ({}, {}))
    monkeypatch.setenv("SLURM_JOB_ID", "fixture")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "32")
    monkeypatch.delenv("SLURM_ARRAY_TASK_ID", raising=False)
    def execute(dataset, order, env, evidence, inputs, provenance):
        calls.append(provenance)
        assert Path.cwd() == launcher
        assert order == [cell["label"]]
        assert dataset["methods"][cell["label"]]["argv"] == cell["argv"]
        assert env["PYTHONPATH"] == str(launcher) and env["OMP_NUM_THREADS"] == "1"
        evidence.mkdir(parents=True)
        return {"failed_methods": [cell["label"]] if outcome == "failed" else []}
    monkeypatch.setattr(module, "execute", execute)
    cwd = Path.cwd()
    if outcome == "failed":
        with pytest.raises(RuntimeError, match="retained"):
            module.run(tmp_path, admission_path, "fixture", "1", env_path, 0)
    else:
        result = module.run(tmp_path, admission_path, "fixture", "1", env_path, 0, outcome == "preflight")
        assert result["accuracy_evaluated"] is False
    assert Path.cwd() == cwd
    assert len(checks) == (1 if outcome == "preflight" else 2)
    assert len(calls) == (0 if outcome == "preflight" else 1)
    if outcome != "preflight":
        postflight = json.loads((output / "execution" / cell["label"] / "postflight.json").read_text())
        assert bool(postflight["failed_methods"]) is (outcome == "failed")
        assert postflight["native_outputs_validated"] is False


def test_allocation_refused_before_input_access(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError, match="32-CPU"):
        module.run(tmp_path, tmp_path / "absent", "missing", "1", tmp_path / "absent", 0)
    monkeypatch.setenv("SLURM_JOB_ID", "fixture")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "32")
    monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "1")
    with pytest.raises(ValueError, match="Array task"):
        module.run(tmp_path, tmp_path / "absent", "missing", "1", tmp_path / "absent", 0)
