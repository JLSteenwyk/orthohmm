"""Check retained diagnostic output bindings without replaying original analyses."""

from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "benchmark_tools/results"
OUT = RESULTS / "native_qfo_candidate_support_20261006_v1"


def data(path):
    return json.loads(path.read_text())


def table(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def test_all_actual_execution_output_anchors():
    receipt = data(RESULTS / "native_qfo_candidate_support_execution_20261006_v1.json")
    assert receipt["status"] == "actual_export_and_independent_readback_succeeded"
    assert receipt["automatic_retry"] is False and receipt["original_analysis_replayed"] is False
    assert [o["exit_code"] for o in receipt["observations"]] == [0, 0]
    for ref in receipt["outputs"]:
        raw = (ROOT / ref["path"]).read_bytes()
        assert len(raw) == ref["bytes"] and hashlib.sha256(raw).hexdigest() == ref["sha256"]


def test_actual_sources_were_committed_before_selected_execution():
    report = data(OUT / "report.json")
    receipt = data(RESULTS / "native_qfo_candidate_support_execution_20261006_v1.json")
    for key in ("source", "reader_source", "protocol"):
        ref = report[key]
        path = Path(ref["path"])
        committed = subprocess.run(["git", "show", receipt["tested_source_pushed_before_selected_analysis"] + ":"
            + str(path.relative_to(ROOT))], cwd=ROOT, check=True, capture_output=True).stdout
        assert len(committed) == ref["bytes"] and hashlib.sha256(committed).hexdigest() == ref["sha256"]
        assert path.read_bytes() == committed


def test_actual_readback_scope_counts_and_summaries():
    report = data(OUT / "report.json")
    reader = data(RESULTS / "native_qfo_candidate_support_readback_20261006_v1.json")
    assert reader["report"]["sha256"] == hashlib.sha256((OUT / "report.json").read_bytes()).hexdigest()
    assert reader["totals"] == report["totals"]
    assert reader["localized_summary"] == report["localized_summary"]
    assert reader["cohort_counts"] == [{k: s[k] for k in ("iteration", "cohort", "events")} for s in report["cohort_summaries"]]
    assert reader["rational_summaries_verified"] is True and reader["exporter_or_union_kernel_imported"] is False
    for field in ("defaults_changed", "uncertainty_admitted", "calibrated_confidence", "causal_mechanism_established", "publication_ready"):
        assert report[field] is False and reader[field] is False


def test_full_original_pair_identity_and_transitive_blank_features():
    report = data(OUT / "report.json")
    original = table(Path(report["ledger"]["path"]))
    annotated = table(OUT / "annotated_pairs.tsv")
    assert len(original) == len(annotated) == 2295
    assert [{k: row[k] for k in original[0]} for row in annotated] == original
    transitive = [row for row in annotated if row["connection_path"] == "transitive_union"]
    assert len(transitive) == 112
    for row in transitive:
        assert row["first_direct_event"] == "" and row["event_cohort"] == "transitive_no_direct_event"
        assert all(row[k] == "" for k in row if k not in original[0] and k != "event_cohort")


def test_implicated_events_unique_and_actual_cohorts():
    rows = table(OUT / "implicated_events.tsv")
    assert len(rows) == len({r["event_index"] for r in rows}) == 449
    assert Counter(row["cohort"] for row in rows) == {"TP-only": 51, "FP-only": 352, "mixed": 46}
    assert sum(int(r["recovered_tp_pairs"]) for r in rows) == 150
    assert sum(int(r["added_fp_pairs"]) for r in rows) == 2033
    assert Counter((int(r["iteration"]), r["cohort"]) for r in rows) == {
        (0, "TP-only"): 49, (1, "TP-only"): 2, (0, "FP-only"): 322, (1, "FP-only"): 30, (0, "mixed"): 37, (1, "mixed"): 9}


def test_all_event_once_summaries_and_nonfinite_policy():
    report = data(OUT / "report.json")
    summaries = report["cohort_summaries"]
    assert [s["events"] for s in summaries[:4]] == [51, 352, 46, 40241]
    assert sum(s["events"] for s in summaries[:4]) == report["totals"]["accepted_events"] == 40690
    for s in summaries:
        for feature, numeric in s["features"].items():
            assert numeric["missing"] == 0 and numeric["finite"] + numeric["positive_infinity"] == s["events"]
            if feature != "margin":
                assert numeric["positive_infinity"] == 0
    assert [s["features"]["margin"]["positive_infinity"] for s in summaries[:4]] == [12, 55, 10, 20423]
