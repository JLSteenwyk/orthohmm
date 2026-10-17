import copy
import csv
import json
import sys

import pytest

from benchmark_tools import export_hillclimb_diagnostics as export


def qfo():
    scores = {"QfO: " + key: 0.8 for key in export.ENDPOINTS["QfO"]}
    details = {key: {"similarity": 0.8, "assessed_relations": 10.0} if key in ("EC", "GO", "FAS")
               else {"precision": 1.0, "recall": 2 / 3} for key in export.ENDPOINTS["QfO"]}
    cell = dict(dataset="QfO", status="independent_content_readback_verified", scores=scores, details=details,
        coverage={"genes": 12, "species": 3, "groups": 4, "pairs": 50})
    return dict(methods=[dict(key="test:arm", label="Arm", scores=scores, cells=[cell])])


def retained_ob():
    prediction = {"path": "groups", "bytes": 10, "sha256": "identity"}
    cell = dict(dataset="OrthoBench", status="retained_registered_result",
        retained_details={"retained_ob_metadata": {"prediction": prediction}})
    summary = dict(methods=[dict(key="orthofinder_3_1_5_full", label="Full", cells=[cell],
        scores={"OrthoBench: OrthoBench F1": 0.75})])
    retained = dict(rows=[dict(key="orthofinder_3_1_5_full", orthobench_retained_evidence=dict(
        prediction_provenance=prediction, f_score_percent=75, precision_percent=90, recall_percent=65,
        coverage_percent=82.1, genes_in_orthogroups=123, exact_refogs=4))])
    return summary, retained


def test_qfo_f1_and_functional_means_keep_different_diagnostics():
    current, history = export.build(qfo(), "comparison", {}, "retained")
    assert len(current) == 8 and len(history) == 6
    assert current[0]["cell_status"] == "no_admitted_endpoint"
    vgnc = next(r for r in current if r["endpoint"] == "VGNC F")
    assert vgnc["statistic"] == "F1" and vgnc["precision"] == 1
    assert vgnc["native_pairs"] == 50 and vgnc["native_genes"] == 12
    go = next(r for r in current if r["endpoint"] == "GO")
    assert go["statistic"] == "conditional_similarity"
    assert go["assessed_relations"] == 10
    assert go["precision"] is None and go["recall"] is None
    assert "retained_pairs" not in go


def test_retained_ob_percent_conversion_binding_and_semantic_annotation():
    source, retained = retained_ob()
    current, _ = export.build(source, "comparison", retained, "retained")
    row = current[0]
    assert row["score"] == 0.75 and row["precision"] == 0.9 and row["recall"] == 0.65
    assert row["benchmark_reported_gene_coverage"] == 0.821
    assert row["benchmark_genes_in_orthogroups"] == 123
    assert row["native_genes"] is None
    assert row["exact_refogs"] == 4
    assert "later debug" in row["semantics"]
    assert row["detail_source_document"] == "retained"
    assert row["detail_source_location"] == ["rows", 0, "orthobench_retained_evidence"]


def test_retained_ob_supplemental_prediction_binding_supported():
    source, retained = retained_ob()
    evidence = retained["rows"][0]["orthobench_retained_evidence"]
    prediction = evidence.pop("prediction_provenance")
    retained["rows"][0]["orthobench_supplemental_readback"] = {"prediction": prediction}
    assert export.build(source, "comparison", retained, "retained")[0][0]["score"] == 0.75


@pytest.mark.parametrize("change", ["prediction", "score", "duplicate", "missing"])
def test_unbound_retained_ob_evidence_rejected(change):
    source, retained = retained_ob()
    if change == "prediction":
        retained["rows"][0]["orthobench_retained_evidence"]["prediction_provenance"] = {"wrong": 1}
    elif change == "score":
        retained["rows"][0]["orthobench_retained_evidence"]["f_score_percent"] = 50
    elif change == "duplicate":
        retained["rows"].append(copy.deepcopy(retained["rows"][0]))
    else:
        retained["rows"] = []
    with pytest.raises(ValueError):
        export.build(source, "comparison", retained, "retained")


def test_failed_attempts_are_retained_but_cannot_supply_current_scores():
    source = qfo()
    cell = source["methods"][0]["cells"][0]
    cell["status"] = "failed"
    source["methods"][0]["scores"] = {}
    current, history = export.build(source, "comparison", {}, "retained")
    assert len(history) == 6 and all(row["score"] is None for row in history)
    assert all(row["cell_status"] == "failed" for row in history)
    assert all(row["cell_status"] == "no_admitted_endpoint" for row in current)


def test_scored_endpoint_requires_admitted_ownership():
    source = qfo()
    source["methods"][0]["cells"][0]["status"] = "failed"
    with pytest.raises(ValueError, match="lacks admitted"):
        export.build(source, "comparison", {}, "retained")


def test_current_table_uses_latest_admission_not_an_average_or_precision_winner():
    source = qfo()
    later = copy.deepcopy(source["methods"][0]["cells"][0])
    later["details"]["VGNC F"] = {"precision": 0.8, "recall": 0.8}
    source["methods"][0]["cells"].append(later)
    current, history = export.build(source, "comparison", {}, "retained")
    assert len(history) == 12
    row = next(r for r in current if r["endpoint"] == "VGNC F")
    assert row["cell_index"] == 1 and row["precision"] == 0.8


def test_retained_qfo_pair_counts_not_imputed_to_native_coverage():
    source = qfo()
    cell = source["methods"][0]["cells"][0]
    cell["status"] = "retained_registered_result"
    cell["retained_details"] = {"submitted_pairs": 60, "retained_pairs": 50,
        "score_details": {key.removesuffix(" F"): value for key, value in cell.pop("details").items()}}
    cell.pop("coverage")
    row = export.build(source, "comparison", {}, "retained")[0][1]
    assert row["submitted_pairs"] == 60 and row["retained_pairs"] == 50
    assert row["native_pairs"] is None


def test_three_kingdoms_reference_coverage_and_pairs_remain_distinct():
    score = {"ThreeKingdoms: Three Kingdoms F1": 0.8}
    cell = dict(dataset="ThreeKingdoms", status="independent_content_readback_verified", scores=score,
        details={"f_score": 0.8, "precision": 1, "recall": 2 / 3, "reference_gene_coverage": 1,
            "reference_genes": 5, "reference_genes_in_prediction": 5, "true_positive_gene_pairs": 4,
            "false_positive_gene_pairs": 0, "false_negative_gene_pairs": 2, "predicted_orthogroups": 9},
        coverage={"genes": 100, "groups": 9, "pairs": 30, "species": 3})
    source = {"methods": [{"key": "test", "label": "Test", "scores": score, "cells": [cell]}]}
    row = export.build(source, "comparison", {}, "retained")[0][-1]
    assert row["native_pairs"] == 30 and row["true_positive_gene_pairs"] == 4
    assert row["native_genes"] == 100 and row["reference_genes"] == 5
    assert row["reference_gene_coverage"] == 1


@pytest.mark.parametrize("value", [True, -1, 1.01, float("nan"), float("inf"), "0.8"])
def test_bad_unit_values_rejected(value):
    with pytest.raises(ValueError):
        export.unit(value)


@pytest.mark.parametrize("value", [True, -1, 1.5, float("nan"), float("inf"), "10"])
def test_bad_counts_rejected(value):
    with pytest.raises(ValueError):
        export.count(value)


def test_functional_precision_is_rejected():
    source = qfo()
    source["methods"][0]["cells"][0]["details"]["GO"]["precision"] = 0.9
    with pytest.raises(ValueError, match="Functional similarity"):
        export.build(source, "comparison", {}, "retained")


def test_score_disagreement_rejected():
    source = qfo()
    source["methods"][0]["scores"] = dict(source["methods"][0]["scores"], **{"QfO: GO": 0.7})
    with pytest.raises(ValueError, match="does not reproduce"):
        export.build(source, "comparison", {}, "retained")


def test_cli_writes_both_tables_with_na_and_preserves_existing_export(tmp_path, monkeypatch):
    comparison, retained, output = tmp_path / "comparison.json", tmp_path / "retained.json", tmp_path / "export"
    comparison.write_text(json.dumps(qfo()))
    retained.write_text("{}")
    monkeypatch.setattr(sys, "argv", ["export", "--comparison", str(comparison), "--retained-scores", str(retained),
        "--output", str(output)])
    export.main()
    with (output / "endpoint_diagnostics.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == 8 and rows[0]["score"] == "NA"
    assert rows[1]["precision"] == "1.0" and rows[4]["precision"] == "NA"
    report = json.loads((output / "report.json").read_text())
    assert len(report["history"]) == 6
    assert not report["defaults_promoted"] and not report["goal_completed"]
    with pytest.raises(FileExistsError):
        export.main()
