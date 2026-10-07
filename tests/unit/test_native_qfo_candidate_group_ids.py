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
