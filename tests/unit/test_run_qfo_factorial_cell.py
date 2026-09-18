from copy import deepcopy
from pathlib import Path
import sys

import pytest

from benchmark_tools.prepare_orthobench_factorial import plan_cells
from benchmark_tools.run_qfo_factorial_cell import select_cell, native_command


def manifest(tmp_path):
    executor, output, fasta = (tmp_path / s for s in ("executor", "output", "fasta"))
    arms = {f"p{p}_c{c}": {"candidate_partition": {"path": str(output / "candidates" / f"p{p}_c{c}" /
            "orthohmm_working_res/orthohmm_edges_clustered.txt")}} for p in (0, 1) for c in (0, 1)}
    return {"status": "qfo_four_candidate_arms_prepared_unscored", "accuracy_computed": False,
            "source": {"path": str(executor / "benchmark_tools/prepare_qfo_factorial.py")},
            "input_fastas": [{"path": str(fasta / "one.fasta")}], "candidate_arms": arms,
            "cells": plan_cells(output, fasta, executor / "benchmark_tools/replay_phylogeny.py", 32)}


def test_all_four_reconciliation_cells(tmp_path):
    data = manifest(tmp_path)
    assert [select_cell(data, i)[0]["label"] for i in range(4)] == ["p0_c0_r1", "p0_c1_r1", "p1_c0_r1", "p1_c1_r1"]
    assert select_cell(data, 0)[1] == tmp_path / "output"


@pytest.mark.parametrize("mutation", ["status", "scored", "arm", "cpu", "constraints", "extra_flag", "fasta"])
def test_changed_design_rejected(tmp_path, mutation):
    data = deepcopy(manifest(tmp_path))
    if mutation == "status":
        data["status"] = "preparing"
    elif mutation == "scored":
        data["accuracy_computed"] = True
    elif mutation == "arm":
        del data["candidate_arms"]["p1_c1"]
    elif mutation == "cpu":
        argv = data["cells"][1]["argv"]
        argv[argv.index("--cpu") + 1] = "64"
    elif mutation == "constraints":
        data["cells"][3]["argv"] = data["cells"][3]["argv"][:-2]
    elif mutation == "extra_flag":
        data["cells"][1]["argv"].append("--unconstrained-membership")
    else:
        data["input_fastas"].append({"path": "/other/two.fasta"})
    with pytest.raises(ValueError):
        select_cell(data, 0)


@pytest.mark.parametrize("index", [-1, 4, True, 0.5])
def test_invalid_index_rejected(tmp_path, index):
    with pytest.raises(ValueError):
        select_cell(manifest(tmp_path), index)


def test_launcher_relocation_changes_only_source_path(tmp_path, monkeypatch):
    import benchmark_tools.run_qfo_factorial_cell as module
    executor, launcher = tmp_path / "executor", tmp_path / "frozen"
    cell = select_cell(manifest(tmp_path), 3)[0]
    monkeypatch.setattr(module, "record", lambda path: {"path": str(path), "sha256": Path(path).name, "bytes": 123})
    argv, sources = native_command(cell, launcher, executor)
    assert argv == [sys.executable, str(launcher / "benchmark_tools/replay_phylogeny.py"), *cell["argv"][2:]]
    assert len(sources) == 2
    monkeypatch.setattr(module, "record", lambda path: {"path": str(path), "sha256": str(path), "bytes": 123})
    with pytest.raises(ValueError, match="differs"):
        native_command(cell, launcher, executor)
