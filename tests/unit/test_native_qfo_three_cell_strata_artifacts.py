"""Bindings and numeric claims of the actual retained three-cell projection."""

import json
from pathlib import Path
import subprocess

import pytest

from benchmark_tools import export_native_qfo_three_cell_strata as exporter
from benchmark_tools import readback_native_qfo_three_cell_strata_v2 as reader

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "benchmark_tools/results"
REPORT = RESULTS / "native_qfo_three_cell_strata_20261007_v1/report.json"
READBACK = RESULTS / "native_qfo_three_cell_strata_readback_20261007_v2.json"


def test_actual_report_reader_and_every_direct_byte_binding():
    assert exporter.record(REPORT)["sha256"] == "55ed72ce771c756424877cfa412da329ff104efb7dca4b35113b0d80703a59e5"
    assert exporter.record(READBACK)["sha256"] == "6f3eef63c0e96e1d85805fcfb70ed9b66fbb1db76a08ca2a105483422a0f6969"
    report, readback = (json.loads(p.read_text()) for p in (REPORT, READBACK))
    assert readback["report"] == exporter.record(REPORT)
    assert readback["source"] == exporter.record(reader.__file__)
    assert report["source"] == exporter.record(exporter.__file__)
    assert len(readback["checked_inputs"]) == 28
    for ref in readback["checked_inputs"]:
        assert exporter.record(ref["path"]) == ref
    assert len(report["rows"]) == readback["score_rows_checked"] == 60
    assert len(report["differences"]) == readback["differences_checked"] == 40
    assert len(report["family_rows"]) == readback["family_rows_checked"] == 54
    assert readback["inherited_score_rows_reproduced"] == 40
    assert readback["proteins_checked"] == 563
    assert report["cells"][1]["timing_admitted"] is report["cells"][1]["timing_eligible"] is False
    assert all(report[k] is False and readback[k] is False for k in reader.FLAGS)
    for field in ("rows", "differences"):
        assert all(r[k] is None for r in report[field] if not r["families"] for k in reader.METRICS)


def test_actual_tsv_and_human_table_match_verified_json_without_scoring_replay():
    report = json.loads(REPORT.read_text())
    outputs = {Path(r["path"]).name: r for r in report["outputs"]}
    reader.tsv(outputs["scores.tsv"]["path"], report["rows"], reader.SCORE_FIELDS)
    reader.tsv(outputs["differences.tsv"]["path"], report["differences"], reader.DIFF_FIELDS)
    rows = [{"cell": r["cell"], "family": r["family"], **r["counts_without_prior"],
             **{m: r[m] for m in reader.METRICS}} for r in report["family_rows"]]
    reader.tsv(outputs["family_counts.tsv"]["path"], rows, reader.FAMILY_FIELDS)
    reader.human_table(outputs["TABLE.md"]["path"], report)


def test_execution_receipt_matches_observed_sources_commits_and_outputs():
    receipt = json.loads((RESULTS / "native_qfo_three_cell_strata_execution_20261007_v1.json").read_text())
    for name in ("export", "v2_readback"):
        item = receipt[name]
        current = exporter.record(ROOT / item["source_path"])
        assert current["sha256"] == item["source_sha256"]
        assert current["bytes"] == item["source_bytes"]
        assert item["exit_code"] == 0 and item["elapsed_seconds"] == .05 and item["swaps"] == 0
    for name, sha in receipt["outputs"].items():
        assert exporter.record(RESULTS / name)["sha256"] == sha
    assert receipt["export_reexecuted"] is False
    failure = json.loads((ROOT / receipt["failed_v1_readback"]["receipt_path"]).read_text())
    assert failure["report"] == exporter.record(REPORT)
    assert exporter.record(failure["source"]["path"]) == failure["source"]
    assert failure["execution_exit_code"] == 1 and failure["output_written"] is False
    assert not (ROOT / failure["intended_output"]).exists()
    for commit, path in ((receipt["source_commit_before_export"], receipt["export"]["source_path"]),
                         (receipt["source_commit_before_v2_readback"], receipt["v2_readback"]["source_path"])):
        result = subprocess.run(["git", "show", f"{commit}:{path}"], cwd=ROOT, capture_output=True, check=True)
        assert result.stdout == (ROOT / path).read_bytes()


@pytest.mark.parametrize("suite,bin_name,rounded_change", [
    ("sequence", "higher_entropy", "-0.736"), ("sequence", "lower_entropy", "+0.599"),
    ("sequence", "short_relative", "-0.703"), ("sequence", "not_short_relative", "-0.236"),
    ("domain", "median_pfam_types_at_least_two", "+0.643"), ("domain", "median_pfam_types_below_two", "-0.672"),
    ("domain", "repeated_type_fraction_at_least_quarter", "+0.987"),
    ("domain", "repeated_type_fraction_below_quarter", "-0.591"),
    ("duplication", "lower_duplication_fraction", "-0.467"),
    ("duplication", "upper_duplication_fraction", "-0.273"),
])
def test_every_reported_distinct_bin_change_is_bound(suite, bin_name, rounded_change):
    report = json.loads(REPORT.read_text())
    row = [r for r in report["differences"] if r["contrast"] == "C_at_P0_R0"
           and r["suite"] == suite and r["stratum"] == bin_name]
    assert len(row) == 1
    assert format(row[0]["F1"] * 100, "+.3f") == rounded_change
    assert rounded_change in (RESULTS / "NATIVE_QFO_THREE_CELL_STRATA_RESULT_20261007.md").read_text()


def test_summary_tradeoffs_match_all_nonempty_bins():
    report = json.loads(REPORT.read_text())
    candidate = [r for r in report["differences"] if r["contrast"] == "C_at_P0_R0" and r["families"]]
    assert len(candidate) == 15
    assert all(r["PPV"] < 0 and r["TPR"] > 0 for r in candidate)
    overall = next(r for r in candidate if r["suite"] == "sequence" and r["stratum"] == "all")
    assert format(-overall["PPV"] * 100, ".3f") == "3.523"
    assert format(overall["TPR"] * 100, ".3f") == "4.375"
    assert format(overall["F1"] * 100, "+.3f") == "-0.347"


def test_old_main_and_rc5_index_preserved():
    assert exporter.record(RESULTS / "PUBLICATION_MAIN_TEXT_20261006_v2.md")["sha256"] == (
        "2804448cbd99ffa3164fcb9e997640f966103a827c9da761bfb7c3732d79a9c1"
    )
    assert exporter.record(RESULTS / "publication_package_rc5_index_20261006.json")["sha256"] == (
        "66bc02d7821d1097808a870c38198bded9c88c379c24052b1828867af8692922"
    )
