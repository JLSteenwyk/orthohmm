import copy
import json

import pytest

from benchmark_tools.export_qfo_corrected_factorial import CELLS, CONVERTERS, ENDPOINTS, export, extract
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def fixture(index):
    kind, _, semantics = CONVERTERS[index % 2]
    participant = "ohmm_qfo_corrected_factorial_" + CELLS[index]
    marker = {"path": "/fixture/recheck", "bytes": 1, "sha256": "fixture"}
    conversion = {"status": f"corrected_factorial_{kind}_pairs_prepared_unscored", "index": index,
        "cell": CELLS[index], "participant": participant, "semantics": semantics,
        "accuracy_evaluated": False, "publication_ready": False, "job_id": "123",
        "total_pairs": 2, "retained_pairs": 2, "removed_mapping_pairs": 0, "expected_pairs": 2,
        "pairs": marker, "filtered_pairs": marker, "native_admission_recheck": marker,
        "checked_records": [marker]}
    return {"status": "corrected_factorial_assessment_admitted", "index": index, "cell": CELLS[index],
        "accuracy_admitted": True, "publication_ready": False, "conversion": conversion,
        "conversion_scheduler": {"JobIDRaw": "123", "State": "COMPLETED", "ExitCode": "0:0",
                                 "NodeList": "bizon", "AllocCPUS": "2"},
        "assessment": {"participant": participant, "secondary_six_metric_mean": .5,
            "endpoints": {e: {"native_participant": {"participant_id": participant,
                "metric_x": .5 if e in ("VGNC", "SwissTrees", "TreeFam-A") else 10, "metric_y": .5},
                "score": .5, "score_semantics": "fixture statistic"} for e in ENDPOINTS}}}


@pytest.mark.parametrize("index", range(8))
def test_all_cells(index):
    report = fixture(index)
    row = extract(report, report["conversion"])
    assert row["cell"] == CELLS[index]
    assert row["secondary_mean"] == .5
    assert row["prediction_semantics"] == CONVERTERS[index % 2][2]


@pytest.mark.parametrize("problem", ["historical", "unadmitted", "cell", "index", "conversion",
    "loss", "semantics", "running", "participant", "missing_endpoint", "nan", "bool", "recall", "score", "mean"])
def test_corruption(problem):
    report = fixture(7)
    conversion = copy.deepcopy(report["conversion"])
    endpoint = report["assessment"]["endpoints"]["SwissTrees"]
    if problem == "historical":
        report["status"] = "original_factorial_assessment_admitted"
    elif problem == "unadmitted":
        report["accuracy_admitted"] = False
    elif problem == "cell":
        report["cell"] = CELLS[0]
    elif problem == "index":
        report["index"] = True
    elif problem == "conversion":
        conversion["total_pairs"] = 3
    elif problem in ("loss", "semantics"):
        key, value = ("removed_mapping_pairs", 1) if problem == "loss" else ("semantics", CONVERTERS[0][2])
        conversion[key] = report["conversion"][key] = value
    elif problem == "running":
        report["conversion_scheduler"]["State"] = "RUNNING"
    elif problem == "participant":
        endpoint["native_participant"]["participant_id"] = "historical"
    elif problem == "missing_endpoint":
        report["assessment"]["endpoints"].pop("FAS")
    elif problem in ("nan", "bool", "recall"):
        endpoint["native_participant"]["metric_x"] = {"nan": float("nan"), "bool": True, "recall": 1.1}[problem]
    elif problem == "score":
        endpoint["score"] = .7
    elif problem == "mean":
        report["assessment"]["secondary_six_metric_mean"] = .7
    with pytest.raises(ValueError):
        extract(report, conversion)


def test_empty_export(tmp_path):
    output = tmp_path / "empty"
    result = export([], output)
    assert result["admitted_cells"] == 0
    assert len(result["cells"]) == 8
    assert all(r["scores"]["SwissTrees"] is None for r in result["cells"])
    assert "not admitted" in (output / "scores.md").read_text()
    with pytest.raises(FileExistsError):
        export([], output)
    broken = tmp_path / "broken"
    broken.symlink_to(tmp_path / "absent")
    with pytest.raises(FileExistsError):
        export([], broken)


def test_partial_export_binding(tmp_path):
    report = fixture(7)
    pairs = tmp_path / "pairs.json"
    pairs.write_text(json.dumps(report["conversion"]))
    report["pairs_manifest"] = record(pairs)
    admission = tmp_path / "admission.json"
    admission.write_text(json.dumps(report))
    source = (admission, record(admission)["sha256"])
    result = export([source], tmp_path / "partial")
    assert result["admitted_cells"] == 1
    assert result["cells"][7]["scores"]["SwissTrees"] == .5
    assert result["cells"][0]["scores"]["SwissTrees"] is None
    for item in result["outputs"]:
        check(item)
    with pytest.raises(ValueError):
        export([source, source], tmp_path / "duplicate")
    with pytest.raises(ValueError):
        export([(admission, "0" * 64)], tmp_path / "wrong_hash")
    pairs.write_text("{}")
    with pytest.raises(ValueError):
        export([source], tmp_path / "corrupt")
    assert not (tmp_path / "corrupt").exists()
