from pathlib import Path

import pytest

from benchmark_tools.run_qfo_recovered_assessment import select_stage, command_for, checked_record, STAGES
from benchmark_tools.prepare_qfo_recovered_pairs import record


def fixture():
    return {"status": "four_stage_pairs_prepared_unscored", "accuracy_evaluated": False,
            "stages": [{"index": i, "stage": s, "participant": f"ohmm_checked_v2_{i}",
                        "retained_pairs": 9, "total_pairs": 10, "removed_mapping_pairs": 1,
                        "filtered_pairs": {"path": f"/pairs/{i}"}} for i, s in enumerate(STAGES)]}


@pytest.mark.parametrize("index", [-1, 4, True, "1"])
def test_bad_index_rejected(index):
    with pytest.raises(ValueError):
        select_stage(fixture(), index)


@pytest.mark.parametrize("problem", [None, "missing", "count", "identity", "scored", "order"])
def test_inventory_and_counts(problem):
    pairs = fixture()
    if problem == "missing":
        pairs["stages"].pop()
    elif problem == "count":
        pairs["stages"][0]["removed_mapping_pairs"] = 0
    elif problem == "identity":
        pairs["stages"][0]["participant"] = "historical"
    elif problem == "scored":
        pairs["accuracy_evaluated"] = True
    elif problem == "order":
        pairs["stages"].reverse()
    if problem:
        with pytest.raises(ValueError):
            select_stage(pairs, 0)
    else:
        assert select_stage(pairs, 0)["stage"] == "multipass"


def test_six_challenges_no_resume_and_local_configuration():
    manifest = {"pipeline": "/pipeline", "execution_config": {"path": "/config"}}
    command = command_for(Path("/root"), fixture()["stages"][0], manifest, Path("/work"), Path("/result"))
    assert "-resume" not in command
    assert command[command.index("--challenges_ids") + 1] == "GO EC VGNC SwissTrees TreeFam-A FAS"
    assert command[command.index("-c") + 1] == "/config"
    with pytest.raises(ValueError, match="Darwin"):
        command_for(Path("/root"), fixture()["stages"][0], manifest, Path("/" + "x" * 160), Path("/result"))


def test_drift_rejected(tmp_path):
    p = tmp_path / "input"
    p.write_text("before")
    evidence = record(p)
    checked_record(evidence)
    p.write_text("after")
    with pytest.raises(ValueError, match="checksum"):
        checked_record(evidence)
