import gzip
import json
from pathlib import Path
import re
import tarfile

from benchmark_tools.audit_root_context_native_session import validate
from benchmark_tools.verify_root_context_native_provenance import scheduler_identity

RESULTS = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


def test_retained_native_audit_and_description_keep_all_results_and_flags():
    audit = json.loads(gzip.decompress((RESULTS / "root_context_native_audit_22021_20260919.json.gz").read_bytes()))
    description = json.loads(gzip.decompress((RESULTS / "root_context_native_description_22021_20260919.json.gz").read_bytes()))
    assert audit["validated_tasks"] == 3 and audit["all_tasks_validated"] is True
    assert audit["issues"] == [] and audit["all_outputs_equivalent"] is True
    assert audit["scientific_timings_admitted"] is False and audit["publication_ready"] is False
    assert [r["native_wall_s"] for r in audit["runs"]] == [551.830510409, 809.677239777, 613.159839346]
    assert [len(r["original_flagged_intervals"]) for r in audit["runs"]] == [0, 51, 1]
    assert [len(r["narrow_flagged_intervals"]) for r in audit["runs"]] == [0, 10, 1]
    assert [r["native_counts"]["input_genes"] for r in audit["runs"]] == [73266]*3
    assert [r["subsets"]["all"]["intervals"] for r in description["runs"]] == [552, 810, 614]
    for source, described in zip(audit["runs"], description["runs"]):
        assert source["work_identity"] == source["prior_work_identity"]
        assert described["subsets"]["all"]["original_flags"] == len(source["original_flagged_intervals"])
        assert described["subsets"]["all"]["narrow_flags"] == len(source["narrow_flagged_intervals"])
        assert described["subsets"]["all"]["distributions"]["root_minus_three_named_children_cpu_usec"]["minimum"] < 0
    assert description["runs"][0]["subsets"]["narrow_flagged"]["distributions"] is None


def test_retained_wait_receipts_match_real_terminal_scheduler():
    with tarfile.open(RESULTS / "root_context_native_receipts_22021.tar.gz") as archive:
        receipts = [json.load(archive.extractfile("root_context_native_submission_v1/" + name))
                    for name in ("queue.json", "launch.json", "result.json")]
        raw = archive.extractfile("root_context_native_scheduler_22021/scheduler_22021.txt").read().decode()
        journal = archive.extractfile("root_context_native_manager_journal_22021.txt").read().decode()
    scheduler_identity(raw, 22021)
    allocation = dict(re.findall(r"(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)", raw))
    recipe = json.loads((RESULTS / "dgx_root_context_native_recipe_20260919.json").read_text())
    result = validate(*receipts, allocation, recipe,
        "2190269735646d5996832f841c4439d38e0b7f05fc5b3402af000f8b2bfe73f4", 22021, "America/New_York")
    assert result["status"] == "native_bounded_session_verified"
    restarts = [line for line in journal.splitlines() if "Scheduled restart job" in line
        and "2026-09-19T16:53:58-07:00" <= line.split()[0] <= "2026-09-19T17:27:32-07:00"]
    assert len(restarts) == 383


def test_relocated_native_audit_results_match():
    result = json.loads((RESULTS / "root_context_native_relocation_22021_20260919.json").read_text())
    assert result["same_results"] is True
    assert result["original_validated_tasks"] == result["relocated_validated_tasks"] == 3
    assert result["scientific_timings_admitted"] is False
