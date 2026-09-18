from pathlib import Path

import pytest

from benchmark_tools.prepare_qfo_sequence_search_control import STAGING_SHA, validate_inputs, search_plan, prepare
from benchmark_tools.prepare_sequence_search_control import search_command


def fixture():
    inputs = [{"path": f"/corrected/species_{i:02d}.fa"} for i in range(78)]
    return {"input_fastas": inputs}, {"status": "corrected_staged_inventory_verified_pending_execution_freeze",
        "total_sequences": 984137, "proteomes": 78, "inputs": [{"sha256": STAGING_SHA}],
        "files": [{"file": r} for r in inputs]}


def test_corrected_plan_preserves_all_target_commands():
    stage, inventory = fixture()
    inputs = validate_inputs(stage, inventory)
    plan = search_plan(inputs, Path("/diamond"), Path("/queries.fa"), Path("/new"))
    assert len(plan) == 78
    for i, row in enumerate(plan):
        assert row["index"] == i
        assert row["target_fasta"] == inputs[i]
        target = Path(f"/new/target_{i:02d}")
        assert row["search"] == search_command("/diamond", "/queries.fa", target / "target", target / "hits.tsv")
        assert row["makedb"] == ["/diamond", "makedb", "--in", inputs[i]["path"], "--db", str(target / "target"), "--threads", "32"]


@pytest.mark.parametrize("problem", ["status", "genes", "species", "hash", "mismatch", "duplicate", "parent"])
def test_wrong_inventory(problem):
    stage, inventory = fixture()
    if problem == "status":
        inventory["status"] = "historical"
    elif problem == "genes":
        inventory["total_sequences"] = 976504
    elif problem == "species":
        inventory["proteomes"] = 12
    elif problem == "hash":
        inventory["inputs"][0]["sha256"] = "wrong"
    elif problem == "mismatch":
        inventory["files"].pop()
    elif problem == "duplicate":
        stage["input_fastas"][0]["path"] = stage["input_fastas"][1]["path"]
    else:
        stage["input_fastas"][0]["path"] = "/historical/species.fa"
    with pytest.raises(ValueError):
        validate_inputs(stage, inventory)


def test_existing_destination_not_reused(tmp_path):
    with pytest.raises(FileExistsError):
        prepare(tmp_path, tmp_path)
