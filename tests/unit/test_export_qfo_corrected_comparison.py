import csv
import json

import pytest

from benchmark_tools.export_qfo_corrected_comparison import ENDPOINTS, export, extract
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def fixture():
    participant = "qfo_corrected_proteinortho"
    endpoints = {e: {"score": .5, "score_semantics": "fixture",
        "native_participant": {"participant_id": participant, "metric_x": .5 if e in ("VGNC", "SwissTrees", "TreeFam-A") else 10,
                               "metric_y": .5}} for e in ENDPOINTS}
    report = {"status": "corrected_comparator_assessment_admitted", "accuracy_admitted": True,
        "publication_ready": False, "method": "proteinortho", "assessment": {
            "participant": participant, "endpoints": endpoints, "secondary_six_metric_mean": .5}}
    conversion = {"status": "corrected_comparator_pairs_prepared_unscored", "method": "proteinortho",
        "participant": participant, "total_pairs": 12, "retained_pairs": 10, "removed_mapping_pairs": 2,
        "semantics": "native post-clustering pairs"}
    return report, conversion


def test_native_arithmetic():
    report, conversion = fixture()
    result = extract(report, conversion)
    assert result["scores"] == dict.fromkeys(ENDPOINTS, .5)
    assert result["details"]["SwissTrees"] == {"recall": .5, "precision": .5, "statistic": "F1"}
    assert result["retained_pairs"] == 10


@pytest.mark.parametrize("problem", ["unadmitted", "old_participant", "missing_endpoint", "mixed_participant",
    "nan", "negative", "recall_above_one", "wrong_f1", "wrong_mean", "counts", "boolean_count", "unknown_method"])
def test_reject_invalid_scores(problem):
    report, conversion = fixture()
    endpoints = report["assessment"]["endpoints"]
    if problem == "unadmitted":
        report["accuracy_admitted"] = False
    elif problem == "old_participant":
        report["assessment"]["participant"] = "proteinortho_native"
    elif problem == "missing_endpoint":
        endpoints.pop("FAS")
    elif problem == "mixed_participant":
        endpoints["GO"]["native_participant"]["participant_id"] = "another_method"
    elif problem == "nan":
        endpoints["GO"]["native_participant"]["metric_y"] = float("nan")
    elif problem == "negative":
        endpoints["GO"]["native_participant"]["metric_x"] = -1
    elif problem == "recall_above_one":
        endpoints["SwissTrees"]["native_participant"]["metric_x"] = 1.01
    elif problem == "wrong_f1":
        endpoints["VGNC"]["score"] = .8
    elif problem == "wrong_mean":
        report["assessment"]["secondary_six_metric_mean"] = .8
    elif problem == "counts":
        conversion["total_pairs"] = 11
    elif problem == "boolean_count":
        conversion["removed_mapping_pairs"] = True
    else:
        report["method"] = "unreviewed_method"
    with pytest.raises(ValueError):
        extract(report, conversion)


def test_export_explicit_missing_rows_and_hashes(tmp_path):
    report, conversion = fixture()
    pairs = tmp_path / "pairs.json"
    pairs.write_text(json.dumps(conversion))
    report["pairs_manifest"] = record(pairs)
    source = tmp_path / "admission.json"
    source.write_text(json.dumps(report))
    output = tmp_path / "table"
    result = export([(source, record(source)["sha256"])], output)
    assert len(result["methods"]) == 8
    missing = [r for r in result["methods"] if r["status"] == "not_admitted"]
    assert len(missing) == 7 and all(set(r["scores"].values()) == {None} for r in missing)
    with (output / "scores.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == 8 and rows[0]["SwissTrees F1"] == ""
    assert "pending" in (output / "scores.md").read_text()
    for item in result["outputs"]:
        assert record(item["path"]) == item
    with pytest.raises(FileExistsError):
        export([(source, record(source)["sha256"])], output)
    with pytest.raises(ValueError, match="Duplicate"):
        export([(source, record(source)["sha256"])] * 2, tmp_path / "duplicates")
    assert not (tmp_path / "duplicates").exists()
    pairs.write_text("{}")
    with pytest.raises(ValueError):
        export([(source, record(source)["sha256"])], tmp_path / "changed")
