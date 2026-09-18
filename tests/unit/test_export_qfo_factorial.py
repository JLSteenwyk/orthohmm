import json
from pathlib import Path

import pytest

from benchmark_tools.export_qfo_factorial import CELLS, extract, export
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

RESULTS = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


def read(index):
    date = "20260918" if index % 2 else "20260917"
    return json.loads((RESULTS / f"qfo_factorial_assessment_{CELLS[index]}_{date}.json").read_text())


@pytest.mark.parametrize("index", range(8))
def test_real_admitted_cells(index):
    result = extract(read(index), index)
    assert result["cell"] == CELLS[index]
    assert len(result["scores"]) == 6
    assert result["total_pairs"] == result["retained_pairs"] + result["removed_mapping_pairs"]
    assert result["secondary_mean"] == pytest.approx(sum(result["scores"].values()) / 6)
    if index == 7:
        assert result["retained_pairs"] == 5646139
        assert result["scores"]["SwissTrees"] == pytest.approx(.7967509658653852)


@pytest.mark.parametrize("problem", ["cell", "index", "admission", "semantics", "pairs", "boolean", "participant",
    "missing", "nonfinite", "recall", "score", "mean", "status"])
def test_reject_misbound_or_corrupt_fresh_results(problem):
    report = read(7)
    conversion = report["conversion"]
    endpoint = report["assessment"]["endpoints"]["SwissTrees"]
    if problem == "cell":
        conversion["cell"] = CELLS[0]
    elif problem == "index":
        conversion["index"] = 0
    elif problem == "admission":
        report["accuracy_admitted"] = False
    elif problem == "semantics":
        conversion["semantics"] = "cross-species group-derived clique pairs"
    elif problem == "pairs":
        conversion["retained_pairs"] += 1
    elif problem == "boolean":
        conversion["removed_mapping_pairs"] = True
    elif problem == "participant":
        endpoint["native_participant"]["participant_id"] = "wrong"
    elif problem == "missing":
        report["assessment"]["endpoints"].pop("FAS")
    elif problem == "nonfinite":
        endpoint["native_participant"]["metric_y"] = float("nan")
    elif problem == "recall":
        endpoint["native_participant"]["metric_x"] = 1.1
    elif problem == "score":
        endpoint["score"] += .1
    elif problem == "mean":
        report["assessment"]["secondary_six_metric_mean"] += .1
    elif problem == "status":
        conversion["status"] = "failed"
    with pytest.raises(ValueError):
        extract(report, 7)


def test_reuse_requires_original_identity():
    report = read(0)
    report["original_participant"] = "different"
    with pytest.raises(ValueError):
        extract(report, 0)


def test_export_binds_files_and_preserves_precision(tmp_path):
    entries = []
    for index, cell in enumerate(CELLS):
        report = read(index)
        pairs = tmp_path / (cell + "_pairs.json")
        pairs.write_text(json.dumps(report["stage"] if index in (0, 4) else report["conversion"]))
        report["pairs_manifest"] = record(pairs)
        admission = tmp_path / (cell + ".json")
        admission.write_text(json.dumps(report))
        entries.append({"cell": cell, "admission": record(admission)})
    inventory = tmp_path / "inventory.json"
    inventory.write_text(json.dumps({"cells": entries}))
    sha = record(inventory)["sha256"]
    output = tmp_path / "table"
    result = export(inventory, sha, output)
    assert len(result["cells"]) == 8
    assert result["cells"][7]["scores"]["SwissTrees"] == .7967509658653852
    for item in result["outputs"]:
        check(item)
    assert "not corrected-input results" in (output / "scores.md").read_text()
    with pytest.raises(FileExistsError):
        export(inventory, sha, output)
    with pytest.raises(ValueError):
        export(inventory, "0" * 64, tmp_path / "bad_hash")
    pairs.write_text("{}")
    with pytest.raises(ValueError):
        export(inventory, sha, tmp_path / "changed_conversion")
    assert not (tmp_path / "changed_conversion").exists()
