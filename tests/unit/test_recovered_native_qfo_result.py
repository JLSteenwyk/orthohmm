"""Portable checks on retained publication artifacts, not new raw admissions."""

import csv
import hashlib
import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
DIRECTORY = ROOT / "native_qfo_scientific_scores_20261006_v1"


def read(name):
    return json.loads((ROOT / name).read_text())


def test_published_admission_exact_bytes_and_failed_timing_preserved():
    data = (ROOT / "recovered_native_qfo_assessment_admission_22449.json").read_bytes()
    assert len(data) == 781476
    assert hashlib.sha256(data).hexdigest() == "4c3a17a76eff5043c8b40d0c9f1e8ead6c7ed32d4e248dbd485988c33537ae3b"
    report = json.loads(data)
    assert report["accuracy_admitted"] is True and report["resources"] is None
    assert report["native_job_id"] == 22437 and report["native_index"] == 7
    for key in ("original_native_scheduler_success", "scientific_timings_admitted",
                "eligible_for_timing_comparison", "native_inference_reexecuted", "publication_ready"):
        assert report[key] is False


def test_retained_fresh_tasks_and_native_records_complete():
    report = read("recovered_native_qfo_assessment_admission_22449.json")
    assert len(report["native_tasks"]) == 15
    assert all(row["status"] == "COMPLETED" for row in report["native_tasks"])
    assert len(report["assessment"]["native_assessments"]) == 48
    assert len(report["checked_records"]) == 1702


def test_portable_table_fields_match_snapshot_without_raw_replay():
    snapshot = read("native_qfo_scientific_scores_20261006_v1/report.json")
    with (DIRECTORY / "scores.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert snapshot["supplied_admissions"] == 2 and snapshot["supplied_recovered_admissions"] == 1
    assert len(rows) == len(snapshot["rows"]) == 7 and len(rows[0]) == 16
    fields = {"Index": "index", "Cell": "cell", "Accuracy status": "status",
        "Measurement status": "measurement_status", "Secondary mean": "secondary_mean",
        "Submitted pairs": "submitted_pairs", "Input accessions": "input_accessions",
        "Relation accessions": "relation_accessions", "Relation coverage": "relation_coverage",
        "Prediction semantics": "prediction_semantics"}
    metrics = {"VGNC F1": "VGNC", "SwissTrees F1": "SwissTrees", "TreeFam-A F1": "TreeFam-A",
        "GO similarity": "GO", "EC similarity": "EC", "FAS": "FAS"}
    for actual, expected in zip(rows, snapshot["rows"]):
        for label, key in fields.items():
            assert actual[label] == ("" if expected[key] is None else str(expected[key]))
        for label, key in metrics.items():
            value = expected["scores"][key]
            assert actual[label] == ("" if value is None else str(value))
    assert all(value is None for row in snapshot["rows"][2:] for value in row["scores"].values())


def test_summary_scores_and_coverage_are_bound_to_machine_results():
    snapshot = read("native_qfo_scientific_scores_20261006_v1/report.json")
    text = (ROOT / "RECOVERED_NATIVE_QFO_SCORE_RESULT_22449.md").read_text()
    for name, label in (("VGNC", "VGNC F1"), ("SwissTrees", "SwissTrees F1"),
        ("TreeFam-A", "TreeFam-A F1"), ("GO", "GO similarity"), ("EC", "EC similarity"), ("FAS", "FAS")):
        values = [row["scores"][name] for row in snapshot["rows"][:2]]
        assert f"| {label} | {values[0]:.6f} | {values[1]:.6f} |" in text
    a, b = snapshot["rows"][:2]
    assert f'| Secondary six-metric mean | {a["secondary_mean"]:.6f} | {b["secondary_mean"]:.6f} |' in text
    assert f'| Submitted pairs | {a["submitted_pairs"]:,} | {b["submitted_pairs"]:,} |' in text
    assert f'| Inputs with any relation | {100*a["relation_coverage"]:.4f}% | {100*b["relation_coverage"]:.4f}% |' in text
    assert "not selected tool defaults" in text


def test_exact_family_records_permit_only_one_retained_contrast():
    audit = read("recovered_native_qfo_swiss_counts_22449_20261006.json")
    retained = read("qfo_corrected_factorial_complete_20260919/swiss_counts.json")
    native = audit["cells"][0]
    original = next(row for row in retained["cells"] if row["cell"] == native["cell"])
    assert native["families"] == original["families"] and native["aggregate"] == original["aggregate"]
    assert len(native["families"]) == 18 and audit["reference_relation_count"] == 10765
    binding = read("native_qfo_swiss_uncertainty_binding_22449_20261006.json")
    assert len(binding["bound_cells"]) == 2 and len(binding["contrasts"]) == 14
    assert binding["new_bootstrap_draws"] == 0 and binding["multiplicity_endpoints"] == 42
    assert [r["name"] for r in binding["contrasts"] if r["status"] == "native_records_matched"] == ["R_at_P0_C0"]
    assert all(r["metrics"] is None for r in binding["contrasts"] if r["status"] != "native_records_matched")


def test_summary_intervals_preserve_zero_crossing_and_family_signs():
    binding = read("native_qfo_swiss_uncertainty_binding_22449_20261006.json")
    effect = next(row for row in binding["contrasts"] if row["name"] == "R_at_P0_C0")
    text = (ROOT / "RECOVERED_NATIVE_QFO_SCORE_RESULT_22449.md").read_text()
    for name, label in (("F1", "F1"), ("PPV", "Precision"), ("TPR", "Recall")):
        metric = effect["metrics"][name]
        lo, hi = metric["bonferroni_percentile_ci"]
        signs = f'{metric["family_wins"]}/{metric["family_ties"]}/{metric["family_losses"]}'
        assert f'| {label} | {100*metric["difference"]:+.4f} | [{100*lo:.4f}, {100*hi:.4f}] | {signs} |' in text
    lo, hi = effect["metrics"]["F1"]["bonferroni_percentile_ci"]
    assert lo < 0 < hi and "not a clear adjusted" in text


def test_independent_readbacks_keep_scopes_and_failures_explicit():
    score = read("recovered_native_qfo_score_readback_22449_20261006.json")
    swiss = read("recovered_native_qfo_swiss_readback_22449_20261006.json")
    assert score["table_fields_checked"] == 112 and score["scientific_timings_admitted"] is False
    assert score["new_fas_sampling"] is False and score["independent_biological_validation"] is False
    assert swiss["raw_pairs_checked"] == 10765 and swiss["raw_families_checked"] == 18
    assert swiss["ordinary_raw_recounted"] is False and swiss["new_bootstrap_draws"] == 0
    assert swiss["publication_ready"] is False


def test_result_links_exist_and_timing_caveat_is_present():
    text = (ROOT / "RECOVERED_NATIVE_QFO_SCORE_RESULT_22449.md").read_text()
    links = re.findall(r"\[[^]]+\]\(([^)]+)\)", text)
    assert links and all((ROOT / link).is_file() for link in links)
    assert "unknown and potentially tool-dependent impact" in text
    assert "not estimates of isolated performance" in text
    assert "Native pair-IID SEM" in text and "not a paired-method or family CI" in text
