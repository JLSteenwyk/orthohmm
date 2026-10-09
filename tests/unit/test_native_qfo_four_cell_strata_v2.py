"""Exercise the actual legacy vocabulary using invented scientific records."""

import copy
import json
from pathlib import Path

import pytest

from benchmark_tools import export_native_qfo_four_cell_strata_v2 as current
from benchmark_tools import readback_native_qfo_four_cell_strata_v2 as independent
from tests.unit.test_native_qfo_four_cell_strata import fixture, ROOT


def legacy_fixture():
    docs = fixture()
    for row in docs["distance"]["rows"]:
        if row["cell"] == current.legacy.CELLS[1]:
            row["prediction_semantics"] = "native_pair"
    return docs


def projected(docs):
    return dict(schema="native_qfo_four_cell_strata_v2", **current.legacy.SCOPE, **current.project(docs))


def test_real_vocabulary_regression_preserves_inputs_and_old_arithmetic():
    docs = legacy_fixture()
    before = copy.deepcopy(docs)
    with pytest.raises(ValueError, match="Changed prior score identity"):
        current.legacy.project(docs)
    result = projected(docs)
    assert docs == before
    assert [len(result[k]) for k in ("family_rows", "rows", "differences")] == [72, 92, 69]
    normalized, mapping = current.compatibility_view(docs)
    assert result["compatibility"] == mapping and mapping["row_count"] == 3
    assert {r["stratum"] for r in mapping["rows"]} == set(docs["distance"]["bins"])
    old = current.legacy.project(normalized)
    for field in ("family_rows", "rows", "differences", "cells", "bins"):
        assert result[field] == old[field]
    records = independent.verify(result, docs)
    assert [len(r) for r in records] == [72, 92, 69]
    assert all(row["F1"] is None for row in result["rows"] if row["families"] == 0)


@pytest.mark.parametrize("change", ["fixed_label", "distance_label", "distance_r0", "cell", "stratum", "duplicate",
                                   "prior_score", "prior_unit", "admission", "counts", "timing"])
def test_adapter_refuses_unexpected_originals(change):
    docs = legacy_fixture()
    r1 = next(r for r in docs["distance"]["rows"] if r["cell"] == current.legacy.CELLS[1])
    if change == "fixed_label":
        next(r for r in docs["fixed"]["rows"] if r["cell"] == current.legacy.CELLS[1])["prediction_semantics"] = "native_pair"
    elif change == "distance_label": r1["prediction_semantics"] = "resolved_native_pairs"
    elif change == "distance_r0": docs["distance"]["rows"][0]["prediction_semantics"] = "native_pair"
    elif change == "cell": r1["cell"] = current.legacy.CELLS[3]
    elif change == "stratum": r1["stratum"] = "not_a_bin"
    elif change == "duplicate": docs["distance"]["rows"][1] = copy.deepcopy(docs["distance"]["rows"][0])
    elif change == "prior_score": r1["F1"] += .01
    elif change == "prior_unit": docs["distance"]["differences"][0]["F1_pp"] /= 100
    elif change == "admission": docs["profile"]["cells"][0]["admission"]["sha256"] = "1"*64
    elif change == "counts": docs["profile"]["cells"][0]["families"][0]["counts_without_prior"]["TP"] += 1
    elif change == "timing": docs["snapshot"]["rows"][1]["timing_admitted"] = True
    with pytest.raises(ValueError): current.project(docs)


@pytest.mark.parametrize("change", ["schema", "mapping", "mapping_identity", "mapped_count", "modified_input", "semantics",
                                   "score", "empty", "difference", "prior_identity", "prior_difference", "source_label"])
def test_independent_metadata_boundary_and_arithmetic(change):
    docs = legacy_fixture()
    result = projected(docs)
    if change == "schema": result["schema"] = "native_qfo_four_cell_strata_v1"
    elif change == "mapping": result["compatibility"]["rows"][0]["retained"] = "other"
    elif change == "mapping_identity": result["compatibility"]["rows"][0]["cell"] = current.legacy.CELLS[3]
    elif change == "mapped_count": result["compatibility"]["row_count"] = 0
    elif change == "modified_input": result["compatibility"]["original_input_bytes_changed"] = True
    elif change == "semantics": result["rows"][0]["prediction_semantics"] = "native_pair"
    elif change == "score": result["rows"][0]["F1"] += .01
    elif change == "empty": next(r for r in result["rows"] if not r["families"])["F1"] = 0
    elif change == "difference": result["differences"][0]["TPR"] += .01
    elif change == "prior_identity": docs["distance"]["rows"][0]["family_members"] = []
    elif change == "prior_difference": docs["distance"]["differences"][0]["reference"] = current.legacy.CELLS[3]
    elif change == "source_label": docs["distance"]["rows"][0]["prediction_semantics"] = "unknown"
    with pytest.raises(ValueError): independent.verify(result, docs)


def bound_export(tmp_path, monkeypatch):
    docs = legacy_fixture()
    source = current.legacy.record(__file__)
    refs = {}
    for key in ("fixed", "distance", "profile", "snapshot", "fixed_reader", "distance_reader", "profile_reader"):
        docs[key]["source"] = source
        if key == "fixed_reader": docs[key]["report"] = refs["fixed"]
        if key == "distance_reader": docs[key]["report"] = refs["distance"]
        if key == "profile_reader": docs[key]["audit"] = refs["profile"]
        path = tmp_path / ("invented_" + key + ".json")
        path.write_text(json.dumps(docs[key], sort_keys=True, allow_nan=False))
        refs[key] = current.legacy.record(path)
    protocol = current.legacy.record(ROOT / "benchmark_tools/results" / current.legacy.PROTOCOL)
    checked = [*refs.values(), source, protocol]
    monkeypatch.setattr(current.legacy, "prepare", lambda repo:(docs, refs, checked, protocol))
    monkeypatch.setattr(independent.rational, "PINS", {k:r["sha256"] for k,r in refs.items()})
    output = tmp_path / "invented_v2"
    current.export(ROOT, output)
    return output, current.legacy.record(output / "report.json")


def test_bound_round_trip_checks_truthful_schema_and_tables(tmp_path, monkeypatch):
    output, report = bound_export(tmp_path, monkeypatch)
    receipt = independent.readback(output / "report.json", report["sha256"], tmp_path / "reader.json")
    assert receipt["schema"] == "native_qfo_four_cell_strata_rational_readback_v2"
    assert receipt["report_schema_checked"] == "native_qfo_four_cell_strata_v2"
    assert receipt["internal_arithmetic_contract"] == "native_qfo_four_cell_strata_v1"
    assert receipt["compatibility_rows_checked"] == 3 and receipt["human_rows_checked"] == 23
    assert receipt["inherited_score_rows_reproduced"] == 69 and receipt["inherited_difference_rows_reproduced"] == 46
    assert receipt["publication_ready"] is False and len(receipt["components"]) == 1


@pytest.mark.parametrize("change", ["report_sha", "input", "table", "amendment", "components", "reader_component"])
def test_bound_provenance_refusal(tmp_path, monkeypatch, change):
    output, report = bound_export(tmp_path, monkeypatch)
    path = output / "report.json"
    result = json.loads(path.read_text())
    if change == "report_sha": report["sha256"] = "0"*64
    elif change == "input":
        p = Path(result["inputs"]["distance"]["path"])
        p.write_text(p.read_text() + " ")
    elif change == "table":
        p = output / "TABLE.md"
        p.write_text(p.read_text().replace("all", "tampered", 1))
    elif change == "reader_component": monkeypatch.setattr(independent, "READER_COMPONENT_SHA", "0"*64)
    else:
        result["amendment" if change == "amendment" else "components"] = (
            dict(result["amendment"], sha256="0"*64) if change == "amendment" else [])
        path.write_text(json.dumps(result, sort_keys=True, allow_nan=False))
        report = current.legacy.record(path)
    receipt_path = tmp_path / "refused.json"
    with pytest.raises(ValueError): independent.readback(path, report["sha256"], receipt_path)
    assert not receipt_path.exists()


def test_existing_output_rejected_before_loading(tmp_path, monkeypatch):
    monkeypatch.setattr(current.legacy, "prepare", lambda repo:pytest.fail("must not load"))
    with pytest.raises(FileExistsError): current.export(ROOT, tmp_path)
    with pytest.raises(FileExistsError): independent.readback(tmp_path / "absent.json", "0"*64, tmp_path)


def test_real_compatibility_metadata_only_without_projection():
    docs, _, _, _ = current.legacy.prepare(ROOT)
    before = copy.deepcopy(docs)
    view, mapping = current.compatibility_view(docs)
    assert docs == before and len(mapping["rows"]) == 3
    assert all(r["retained"] == "native_pair" and r["normalized"] == "resolved_native_pairs" for r in mapping["rows"])
    assert view["distance"]["family_rows"] == docs["distance"]["family_rows"]
