"""Identifier diagnosis retains failed export without approximate join/retry."""

import json
from pathlib import Path

import pytest

from benchmark_tools import diagnose_native_qfo_candidate_group_ids as diagnosis
from benchmark_tools import trace_native_qfo_candidate_groups as primary
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_native_qfo_candidate_groups import retained


def failed(build, output):
    with pytest.raises(ValueError, match="unmapped changed pair"):
        primary.execute(*build("identifier"), output)
    return record(output / "failure.json")


def test_independent_failure_diagnosis_checks_complete_groups(retained, tmp_path):
    _, build = retained
    original = failed(build, tmp_path / "original")
    result = diagnosis.execute(original, tmp_path / "diagnosis")
    assert result["whole_candidate_partition_reconstructed"]
    assert (result["genes"], result["baseline_groups"], result["candidate_groups"], result["accepted_merges"]) == (6, 5, 2, 3)
    assert result["unmapped_accessions"] == ["z"]
    assert result["changed_pairs"] == 4 and result["pairs_with_unmapped_accessions"] == 1
    assert result["affected_pairs_by_candidate_state"] == {"FP": 1}
    assert result["primary_export_retried"] is result["identifiers_substituted"] is result["pair_localization_admitted"] is False
    assert record(tmp_path / "original/failure.json") == original
    assert not (tmp_path / "original/report.json").exists()
    assert Path(result["unmapped_pair_table"]["path"]).read_text().splitlines()[1] == "t\tz\tnot_scored\tFP\tz"
    with pytest.raises(ValueError, match="fresh direct"):
        diagnosis.execute(original, tmp_path / "diagnosis")


@pytest.mark.parametrize("field,value", [("schema", "wrong"), ("error", "different error"),
    ("error_type", "KeyError"), ("automatic_retry", True), ("accuracy_rescored", True),
    ("native_inference_reexecuted", True), ("publication_ready", True)])
def test_other_failure_or_claims_rejected(retained, tmp_path, field, value):
    write, build = retained
    original = failed(build, tmp_path / "original")
    failure = json.loads(Path(original["path"]).read_text()); failure[field] = value
    ref = write("changed/failure.json", failure)
    with pytest.raises(ValueError):
        diagnosis.execute(ref, tmp_path / "diagnosis")
    result = json.loads((tmp_path / "diagnosis/failure.json").read_text())
    assert result["automatic_retry"] is result["pair_localization_admitted"] is False


def test_partial_primary_report_rejected(retained, tmp_path):
    _, build = retained
    original = failed(build, tmp_path / "original")
    (tmp_path / "original/report.json").write_text("{}")
    with pytest.raises(ValueError, match="unexpected report"):
        diagnosis.collect(original)


def test_stale_input_rejected(retained, tmp_path):
    _, build = retained
    original = failed(build, tmp_path / "original")
    Path(json.loads(Path(original["path"]).read_text())["decomposition"]["path"]).write_text("{}")
    with pytest.raises(ValueError, match="Changed evidence"):
        diagnosis.collect(original)


def test_mapping_must_explain_failure(retained, tmp_path):
    write, build = retained
    original = failed(build, tmp_path / "original")
    failure = json.loads(Path(original["path"]).read_text())
    failure["decomposition"], failure["readback"] = build()
    ref = write("changed/failure.json", failure)
    with pytest.raises(ValueError, match="do not explain"):
        diagnosis.collect(ref)


def test_actual_retained_diagnosis_contract_without_scientific_reexecution():
    root = Path(__file__).resolve().parents[2]
    path = root / "benchmark_tools/results/native_qfo_candidate_group_id_diagnosis_20261006_v1/report.json"
    if not path.is_file():
        pytest.skip("Selected retained diagnosis not installed")
    assert record(path)["sha256"] == "3bc6646e1dc7cb6488ceb6aaaf763e3a4e3cdcbd41268f2066ef2d8008f5be74"
    result = json.loads(path.read_text())
    assert (result["genes"], result["baseline_groups"], result["candidate_groups"], result["accepted_merges"]) == (
        984137, 394328, 353638, 40690)
    assert result["baseline_groups"] - result["candidate_groups"] == result["accepted_merges"]
    assert (result["union_scored_pairs"], result["changed_pairs"], result["pairs_with_unmapped_accessions"],
            result["changed_pairs_without_missing_accessions"]) == (42080, 2295, 8, 2287)
    assert result["unmapped_accessions"] == ["Q17QN5_BOVIN", "Q1RMT5_BOVIN"]
    assert result["affected_pairs_by_candidate_state"] == {"FP": 8}
    assert record(result["failure"]["path"])["sha256"] == "bfad1eefd8f084adc6ab66d5da385c9ce8a1f17440948c59330cf1a761efb695"
    assert record(result["unmapped_pair_table"]["path"])["sha256"] == "6ec0ed60e8004f99ed8822d09ad87da665ddc25883bcf9ffd1b8c073400073b7"
    assert len(Path(result["unmapped_pair_table"]["path"]).read_text().splitlines()) == 9
    for ref in result["checked_records"]:
        assert record(ref["path"]) == ref
    for key in ("pair_localization_admitted", "primary_export_retried", "identifiers_substituted", "accuracy_rescored",
                "native_inference_reexecuted", "uncertainty_admitted", "publication_ready"):
        assert result[key] is False


def test_manuscript_claims_and_result_preserve_join_limitation():
    directory = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
    manuscript = (directory / "PUBLICATION_MANUSCRIPT_DRAFT_20260916.md").read_text()
    claims = (directory / "PUBLICATION_CLAIMS_20260916.md").read_text()
    result = (directory / "NATIVE_QFO_CANDIDATE_GROUP_TRACE_RESULT_20261006.md").read_text()
    assert "Complete pair-path localization\nremains unadmitted" in manuscript
    assert "353,638 candidate groups" in manuscript and "984,137 genes" in manuscript
    assert "No guessed replacement or partial localization is admitted" in claims
    assert "eight changed FP pairs" in claims
    assert "The2,287 other changed pairs are NOT admitted as a complete localization" in result
    assert "Q17QN5_BOVIN" in result and "Q1RMT5_BOVIN" in result
