from copy import deepcopy

import pytest

from benchmark_tools import admit_qfo_corrected_candidates as module
from benchmark_tools.prepare_orthobench_factorial import plan_cells


def fixture(root):
    executor, fasta = root / "executor", root / "inputs"
    output = root / "benchmarks/results/qfo_corrected_factorial_v1"
    arms = {}
    for label in module.ARMS:
        expanded = label.endswith("c1")
        arm = {"candidate_expansion": expanded, "candidate_partition": {
            "path": str(output / "candidates" / label / "orthohmm_working_res/orthohmm_edges_clustered.txt")}}
        if expanded:
            arm.update(membership_constraints={}, expansion={"parameters": deepcopy(module.PARAMETERS)})
        arms[label] = arm
    manifest = {"status": "corrected_qfo_four_candidate_arms_prepared_unscored", "accuracy_computed": False,
        "candidate_arms": arms, "job_id": "20", "cells": plan_cells(output, fasta, executor / "benchmark_tools/replay_phylogeny.py", 32),
        "core_root": str(root / "benchmarks/work/publication_method_native_v2"),
        "launcher_root": str(root / "benchmarks/work/publication_qfo_replay_native_v1"),
        "environment_overrides": {"PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}}
    scheduler = {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "2", "JobIDRaw": "20"}
    return manifest, scheduler, root, executor, fasta


@pytest.mark.parametrize("problem", [None, "old_release", "scored", "missing_arm", "job", "scheduler", "node", "cpus",
    "cell_command", "cell_count", "core", "launcher", "environment", "expansion", "candidate_path", "constraints", "parameters"])
def test_candidate_manifest_gate(tmp_path, problem):
    args = fixture(tmp_path)
    manifest, scheduler, *_ = args
    if problem == "old_release":
        manifest["status"] = "qfo_four_candidate_arms_prepared_unscored"
    elif problem == "scored":
        manifest["accuracy_computed"] = True
    elif problem == "missing_arm":
        manifest["candidate_arms"].pop("p1_c1")
    elif problem == "job":
        manifest["job_id"] = "21"
    elif problem == "scheduler":
        scheduler["State"] = "RUNNING"
    elif problem == "node":
        scheduler["NodeList"] = "other"
    elif problem == "cpus":
        scheduler["AllocCPUS"] = "32"
    elif problem == "cell_command":
        manifest["cells"][1]["argv"][-1] = "changed"
    elif problem == "cell_count":
        manifest["cells"].pop()
    elif problem == "core":
        manifest["core_root"] = "development"
    elif problem == "launcher":
        manifest["launcher_root"] = "development"
    elif problem == "environment":
        manifest["environment_overrides"]["OMP_NUM_THREADS"] = "2"
    elif problem == "expansion":
        manifest["candidate_arms"]["p0_c0"]["candidate_expansion"] = True
    elif problem == "candidate_path":
        manifest["candidate_arms"]["p0_c0"]["candidate_partition"]["path"] = "other"
    elif problem == "constraints":
        manifest["candidate_arms"]["p0_c1"].pop("membership_constraints")
    elif problem == "parameters":
        manifest["candidate_arms"]["p0_c1"]["expansion"]["parameters"]["min_norm"] = .02
    if problem:
        with pytest.raises(ValueError):
            module.validate_manifest(*args)
    else:
        assert module.validate_manifest(*args) == tmp_path / "benchmarks/results/qfo_corrected_factorial_v1"


@pytest.mark.parametrize("state", ["RUNNING", "PENDING", "FAILED", "CANCELLED"])
def test_scheduler_precedes_partial_manifest(tmp_path, monkeypatch, state):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        f"JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n20|{state}|0:0|00:01:00|bizon|2\n")
    output = tmp_path / "admission.json"
    with pytest.raises(ValueError, match="COMPLETED"):
        module.admit(tmp_path, tmp_path / "absent.json", "missing", "20", output)
    assert not output.exists()


def test_existing_report_never_overwritten(tmp_path):
    output = tmp_path / "admission.json"
    output.write_text("preserved")
    with pytest.raises(FileExistsError):
        module.admit(tmp_path, tmp_path / "absent.json", "missing", "20", output)
    assert output.read_text() == "preserved"
