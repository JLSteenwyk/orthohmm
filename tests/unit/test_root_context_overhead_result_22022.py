import gzip
import json
from pathlib import Path
import tarfile

import pytest

from benchmark_tools.audit_root_context_overhead_session import validate
from benchmark_tools.summarize_root_context_overhead import summarize
from benchmark_tools.verify_root_context_overhead_provenance import scheduler_identity

RESULTS = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
SHA = "ef3c4d0e083e31273a098cb342fa89b797911a84f7fe643ca256b27829c40b82"


@pytest.fixture(scope="module")
def audit():
    with gzip.open(RESULTS / "root_context_overhead_audit_22022_20260920.json.gz", "rt") as stream:
        return json.load(stream)


def test_retained_panel_recomputes_all_pairs_and_preserves_flags(audit):
    plan = json.loads((RESULTS / "dgx_root_context_overhead_plan_20260919.json").read_text())
    assert audit["validated_tasks"] == 18 and audit["issues"] == []
    assert summarize(plan, audit["runs"], audit["issues"]) == audit["comparison"]
    assert audit["comparison"]["engineering_budget_passed"] is True
    assert [m["valid_pairs"] for m in audit["comparison"]["methods"]] == [3, 3, 3]
    assert [m["median_signed_ratio"] for m in audit["comparison"]["methods"]] == [
        -0.009296314819665974, 0.013428412156370806, -0.0034416163143453637]
    assert [len(r["original_flagged_intervals"]) for r in audit["runs"]] == [
        1, 0, 49, 57, 1, 2, 59, 53, 2, 2, 0, 0, 1, 4, 1, 0, 60, 52]
    assert [len(r["narrow_flagged_intervals"]) for r in audit["runs"]] == [
        1, 0, 7, 8, 0, 1, 13, 10, 0, 1, 0, 0, 0, 2, 1, 0, 8, 8]
    assert sum(r["observation_intervals"] for r in audit["runs"]) == 11846
    for row in audit["runs"]:
        assert row["status"] == "validated" and row["output_equivalent"] is True
        assert row["work_identity"] == row["prior_work_identity"]
        assert row["native_counts"]["input_genes"] == 73266
        assert (row["root_context"] is None) == (row["arm"] == "lineage")
        assert row["scientific_timings_admitted"] is False
    for key in ("scientific_timings_admitted", "environmental_validity_established", "publication_ready"):
        assert audit[key] is False


def test_retained_session_binds_actual_terminal_scheduler(audit):
    with tarfile.open(RESULTS / "root_context_overhead_receipts_22022.tar.gz") as archive:
        receipts = [json.load(archive.extractfile("root_context_overhead_submission_v1/" + name))
                    for name in ("queue.json", "launch.json", "result.json")]
        raw = archive.extractfile("root_context_overhead_scheduler_22022/scheduler_22022.txt").read().decode()
        capture = json.load(archive.extractfile("root_context_overhead_scheduler_22022/capture.json"))
    allocation = scheduler_identity(raw, 22022)
    recipe = json.loads((RESULTS / "dgx_root_context_overhead_recipe_20260919.json").read_text())
    result = validate(*receipts, allocation, recipe, SHA, 22022, "America/New_York")
    assert result["status"] == "overhead_bounded_session_verified"
    assert all(audit["waiting_session"][k] == v for k, v in result.items())
    assert capture["status"] == "complete" and capture["observation_errors"] == 0
    assert allocation["RunTime"] == "03:21:00"


def test_relocated_audit_matches_all_result_fields():
    result = json.loads((RESULTS / "root_context_overhead_relocation_22022_20260920.json").read_text())
    assert result["same_results"] is True
    assert result["original_validated_tasks"] == result["relocated_validated_tasks"] == 18
    assert result["original_issues"] == result["relocated_issues"] == []
    assert result["excluded_run_keys"] == ["evidence"]
    assert result["scientific_timings_admitted"] is False
