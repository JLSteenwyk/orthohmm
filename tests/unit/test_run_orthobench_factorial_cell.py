from copy import deepcopy
from pathlib import Path

import pytest

from benchmark_tools.prepare_orthobench_factorial import plan_cells
from benchmark_tools.run_orthobench_factorial_cell import select_cell
from benchmark_tools import run_orthobench_factorial_cell as runner


def fixture():
    root = Path("/experiment")
    launcher = Path("/pinned/benchmark_tools/replay_phylogeny.py")
    return {"status": "prepared_not_reconciled", "accuracy_computed": False,
            "launcher": {"path": str(launcher)}, "fasta_inputs": [{"path": "/inputs/a.fasta"}],
            "candidate_arms": {"p0_c0": {"candidate_partition": {"path": "/experiment/candidates/p0_c0/orthohmm_working_res/orthohmm_edges_clustered.txt"}}},
            "cells": plan_cells(root, Path("/inputs"), launcher, 32)}


def test_selects_exact_four_reconciliation_cells():
    manifest = fixture()
    for index, label in enumerate(("p0_c0_r1", "p0_c1_r1", "p1_c0_r1", "p1_c1_r1")):
        cell, root, launcher = select_cell(manifest, index)
        assert cell["label"] == label
        assert root == Path("/experiment") and launcher == Path("/pinned")
        assert "--official-benchmark" not in cell["argv"]
        assert ("--membership-constraints" in cell["argv"]) == cell["candidate_expansion"]


@pytest.mark.parametrize("change", ["cpu", "constraints", "reference", "incomplete", "scored", "index"])
def test_rejects_changed_design(change):
    manifest = deepcopy(fixture())
    index = 0
    if change == "cpu":
        argv = manifest["cells"][1]["argv"]
        argv[argv.index("--cpu") + 1] = "64"
    elif change == "constraints":
        manifest["cells"][3]["argv"][-1] = "/other/constraints.json"
    elif change == "reference":
        manifest["cells"][1]["argv"].extend(["--official-benchmark", "/reference"])
    elif change == "incomplete":
        manifest["status"] = "preparing"
    elif change == "scored":
        manifest["accuracy_computed"] = True
    else:
        index = 4
    with pytest.raises(ValueError):
        select_cell(manifest, index)


@pytest.mark.parametrize("outcome", ["success", "failed_method", "exception"])
def test_execution_restores_verification_directory(tmp_path, monkeypatch, outcome):
    original = tmp_path / "verification"
    launcher = tmp_path / "launcher"
    original.mkdir()
    launcher.mkdir()
    monkeypatch.chdir(original)
    manifest = fixture()
    manifest["environment_overrides"] = {}
    monkeypatch.setattr(runner.sys, "argv", ["runner", "--manifest", "manifest.json",
                        "--manifest-sha256", "hash", "--environment-manifest", "env.json",
                        "--environment-sha256", "hash", "--preparation-job", "1", "--index", "0"])
    monkeypatch.setattr(runner, "read_frozen", lambda *args: manifest)
    cell = manifest["cells"][1]
    monkeypatch.setattr(runner, "select_cell", lambda *args: (cell, tmp_path, launcher))
    checks = []

    def verify(*args):
        checks.append(Path.cwd())
        assert Path.cwd() == original
        return {}

    monkeypatch.setattr(runner, "verify_prepared", verify)
    monkeypatch.setattr(runner, "verify_environment", verify)
    monkeypatch.setattr(runner, "execution_environment", lambda *args: ({}, {}))
    monkeypatch.setattr(runner, "file_provenance", lambda *args: {})

    def execute(*args):
        assert Path.cwd() == launcher
        assert args[-1]["cwd"] == str(launcher)
        if outcome == "exception":
            raise RuntimeError("execution error")
        return {"failed_methods": [cell["label"]] if outcome == "failed_method" else []}

    monkeypatch.setattr(runner, "execute", execute)
    if outcome == "exception":
        with pytest.raises(RuntimeError, match="execution error"):
            runner.main()
    elif outcome == "failed_method":
        with pytest.raises(SystemExit, match="Reconciliation failed"):
            runner.main()
    else:
        runner.main()
    assert Path.cwd() == original
    assert len(checks) == (2 if outcome == "exception" else 4)
