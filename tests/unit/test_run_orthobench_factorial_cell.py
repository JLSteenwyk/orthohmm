from copy import deepcopy
from pathlib import Path

import pytest

from benchmark_tools.prepare_orthobench_factorial import plan_cells
from benchmark_tools.run_orthobench_factorial_cell import select_cell


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
