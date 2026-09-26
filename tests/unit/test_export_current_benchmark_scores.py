import copy
import csv
import json
from pathlib import Path
import shutil

import pytest

from benchmark_tools.export_current_benchmark_scores import (
    SOURCES, OB_AUDITS, METRICS, assemble, export, supplement_orthobench,
)
from benchmark_tools.audit_failed_recovery_refinement import record

RESULTS = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


def reports():
    return [json.loads((RESULTS / name).read_text()) for name in SOURCES]


def test_current_sources_and_exact_units(tmp_path):
    ob, qfo, kingdoms = reports()
    result = export(RESULTS, tmp_path / "scores")
    assert len(result["rows"]) == 8
    assert len(result["inputs"]) == 5
    assert sum("orthobench_supplemental_readback" in row for row in result["rows"]) == 5
    assert (tmp_path / "scores/scores.tsv").read_bytes() == (
        RESULTS / "current_benchmark_scores_20260926/scores.tsv").read_bytes()
    for row, old, corrected in zip(result["rows"], ob["methods"], qfo["methods"]):
        assert row["scores"]["OrthoBench"] == old["orthobench"]["f_score_percent"] / 100
        assert all(row["scores"][metric] == corrected["scores"][metric] for metric in METRICS)
    sonic = next(r for r in result["rows"] if r["key"] == "sonicparanoid_2_0_9")
    assert sonic["three_kingdoms"]["run"] == "contemporary matched input"
    selected = next(r for r in kingdoms["rows"] if r["key"] == sonic["key"])
    assert sonic["scores"]["ThreeKingdoms"] == selected["counts"]["f_score"]
    with (tmp_path / "scores/scores.tsv").open() as stream:
        assert len(list(csv.DictReader(stream, delimiter="\t"))) == 8
    with pytest.raises(FileExistsError):
        export(RESULTS, tmp_path / "scores")


@pytest.mark.parametrize("change", ["original", "duplicate", "missing", "admission", "mean", "nan", "diagnostic"])
def test_rejects_mixed_or_invalid_inputs(change):
    ob, qfo, kingdoms = reports()
    if change == "original":
        qfo["status"] = "historical"
    elif change == "duplicate":
        qfo["methods"].append(copy.deepcopy(qfo["methods"][0]))
    elif change == "missing":
        ob["methods"].pop()
    elif change == "admission":
        qfo["methods"][0]["status"] = "missing"
    elif change == "mean":
        qfo["methods"][0]["secondary_mean"] = 0
    elif change == "nan":
        qfo["methods"][0]["scores"]["GO"] = float("nan")
    else:
        kingdoms["rows"][-1]["use"] = "comparison"
    with pytest.raises(ValueError):
        assemble(ob, qfo, kingdoms)


def test_wrong_source_hash_fails_before_output(tmp_path):
    source = tmp_path / "inputs"
    source.mkdir()
    (source / next(iter(SOURCES))).write_text("{}")
    output = tmp_path / "output"
    with pytest.raises(ValueError, match="Changed retained"):
        export(source, output)
    assert not output.exists()


def audit_inputs():
    paths = [RESULTS / name for name in OB_AUDITS]
    return [json.loads(path.read_text()) for path in paths], record(RESULTS / next(iter(SOURCES))), record(paths[0])


@pytest.mark.parametrize("change", [
    "scope", "agreement", "ready", "source_chain", "duplicate", "missing", "prediction",
    "unchecked", "retained", "exact", "local_score", "upstream_score", "nan", "difference", "metric_count",
])
def test_rejects_inconsistent_supplement_without_partial_attachment(change):
    rows = assemble(*reports())
    (readback, upstream), comparison_record, readback_record = audit_inputs()
    if change == "scope":
        upstream["status"] = "historical"
    elif change == "agreement":
        readback["all_scores_agree"] = False
    elif change == "ready":
        readback["publication_ready"] = True
    elif change == "source_chain":
        upstream["checked_records"] = []
    elif change == "duplicate":
        readback["rows"].append(copy.deepcopy(readback["rows"][0]))
    elif change == "missing":
        upstream["rows"].pop()
    elif change == "prediction":
        upstream["rows"][0]["prediction"]["sha256"] = "0" * 64
    elif change == "unchecked":
        target = upstream["rows"][0]["prediction"]
        upstream["checked_records"] = [r for r in upstream["checked_records"] if r != target]
    elif change == "retained":
        readback["rows"][0]["retained_score"]["coverage_percent"] = -1
    elif change == "exact":
        readback["rows"][0]["score"]["exact_refogs"] = -1
    elif change == "local_score":
        readback["rows"][0]["score"]["precision"] += 1
    elif change == "upstream_score":
        upstream["rows"][-1]["upstream_scores"][0] += 1
    elif change == "nan":
        upstream["rows"][0]["upstream_scores"][0] = float("nan")
    elif change == "difference":
        upstream["rows"][0]["differences"]["precision"] = 1
    else:
        upstream["rows"][0]["upstream_scores"].append(0)
    with pytest.raises(ValueError):
        supplement_orthobench(rows, readback, upstream, comparison_record, readback_record)
    assert all("orthobench_supplemental_readback" not in row for row in rows)


def test_relocated_reports_need_no_raw_predictions_and_reject_changed_audit(tmp_path):
    relocated = tmp_path / "reports"
    for name in {**SOURCES, **OB_AUDITS}:
        destination = relocated / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(RESULTS / name, destination)
    report = export(relocated, tmp_path / "ok")
    for row in report["rows"]:
        if "orthobench_supplemental_readback" in row:
            assert row["orthobench_supplemental_readback"]["historical_consumption_established"] is False
    (relocated / next(iter(OB_AUDITS))).write_text("{}")
    with pytest.raises(ValueError, match="Changed retained score source"):
        export(relocated, tmp_path / "bad")
    assert not (tmp_path / "bad").exists()
