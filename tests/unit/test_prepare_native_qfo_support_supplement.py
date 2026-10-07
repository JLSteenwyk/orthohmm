"""Keep supplementary descriptions separate from tuning and uncertainty claims."""

import copy
import hashlib
from pathlib import Path
import subprocess
import sys

from PIL import Image
import pytest

from benchmark_tools import export_native_qfo_candidate_support as exporter
from benchmark_tools import prepare_native_qfo_support_supplement as presentation


@pytest.fixture
def documents():
    def event(source, target, iteration):
        row = {k: 1 for k in exporter.FEATURES}
        row.update(source_genes=source, target_genes=target, source_size=len(source), target_size=len(target),
                   iteration=iteration, source_cluster=0, target_cluster=1, support=.1, margin=float("inf"))
        return row
    trace = [event(["a", "b"], ["c"], 0), event(["d"], ["e"], 0), event(["f"], ["g"], 1), event(["h"], ["i"], 1)]
    def pair(a, b, state, index, iteration):
        return [a, b, a, b, "FN" if state == "TP" else "not_scored", state, a, b, "new", str(iteration), str(index), "direct_cross_endpoint"]
    rows = [pair("a", "c", "TP", 0, 0), pair("b", "c", "FP", 0, 0), pair("d", "e", "TP", 1, 0), pair("f", "g", "FP", 2, 1)]
    summary = exporter.analyze(trace, rows)[2]
    ref = {"path": "fixture", "bytes": 0, "sha256": "fixture"}
    common = {k: False for k in ("defaults_changed", "uncertainty_admitted", "calibrated_confidence", "causal_mechanism_established", "publication_ready")}
    report = dict(common, schema="native_qfo_candidate_accepted_support_v1", status="descriptive_accepted_event_support_exported",
                  failed_r1_timing_remains_ineligible=True, **summary)
    reader = dict(common, schema="native_qfo_candidate_accepted_support_readback_v1", status="accepted_event_features_independently_verified",
        report=ref, totals=copy.deepcopy(summary["totals"]), localized_summary=copy.deepcopy(summary["localized_summary"]),
        rational_summaries_verified=True, exporter_or_union_kernel_imported=False,
        numeric_absolute_tolerance=1e-12, numeric_relative_tolerance=1e-12,
        cohort_counts=[{k: s[k] for k in ("iteration", "cohort", "events")} for s in summary["cohort_summaries"]])
    return report, reader, ref


def test_valid_descriptive_cohort_round_and_feature_inventory(documents):
    report, reader, ref = documents
    groups = presentation.snapshot(report, reader, ref)
    assert len(groups) == 12 and set(groups) == {(r, c) for r in (None, 0, 1) for c in presentation.COHORTS}
    assert sum(groups[None, c]["events"] for c in presentation.COHORTS) == 4
    assert groups[1, "TP-only"]["events"] == 0
    assert groups[None, "TP-only"]["features"]["margin"]["median"] is None


@pytest.mark.parametrize("field", ["defaults_changed", "uncertainty_admitted", "calibrated_confidence", "causal_mechanism_established", "publication_ready"])
@pytest.mark.parametrize("side", [0, 1])
def test_scientific_claims_are_not_upgraded(documents, side, field):
    documents[side][field] = True
    with pytest.raises(ValueError, match="scope"):
        presentation.snapshot(*documents)


@pytest.mark.parametrize("field,value", [("rational_summaries_verified", False), ("exporter_or_union_kernel_imported", True),
    ("numeric_absolute_tolerance", 1e-6), ("numeric_relative_tolerance", 1e-6), ("report", {})])
def test_independent_verification_cannot_be_reinterpreted(documents, field, value):
    documents[1][field] = value
    with pytest.raises(ValueError):
        presentation.snapshot(*documents)


@pytest.mark.parametrize("change", ["duplicate", "missing", "count", "finite", "sentinel", "nan", "missing_feature", "round_bool"])
def test_partial_or_invalid_summary_rejected(documents, change):
    report, reader, ref = documents
    group = report["cohort_summaries"][0]
    if change == "duplicate":
        report["cohort_summaries"][-1] = copy.deepcopy(group)
    elif change == "missing":
        report["cohort_summaries"].pop()
    elif change == "count":
        group["events"] = True
    elif change == "finite":
        group["features"]["support"]["finite"] = 2
    elif change == "sentinel":
        group["features"]["support"].update(finite=0, positive_infinity=1)
    elif change == "nan":
        group["features"]["support"]["median"] = float("nan")
    elif change == "missing_feature":
        group["features"].pop("margin")
    else:
        group["iteration"] = True
    reader["cohort_counts"] = [{k: s[k] for k in ("iteration", "cohort", "events")} for s in report["cohort_summaries"]]
    with pytest.raises(ValueError):
        presentation.snapshot(report, reader, ref)


def test_invalid_path_counts_rejected_even_if_reader_echoes_them(documents):
    report, reader, ref = documents
    report["localized_summary"][0]["pairs"] += 1
    reader["localized_summary"] = copy.deepcopy(report["localized_summary"])
    with pytest.raises(ValueError, match="Path counts"):
        presentation.snapshot(report, reader, ref)


def test_manuscript_caption_and_claim_to_evidence_scope(documents):
    report, reader, ref = documents
    text = presentation.manuscript(report, presentation.snapshot(report, reader, ref))
    for phrase in ("Every accepted event contributes once", "unlabeled, not true negatives", "not confidence intervals",
        "dependent pair counts", "does not inspect rejected alternatives", "neither changes defaults nor demonstrates superiority",
        "Failed R-on timing remains ineligible", "without reexecuting inference", "Claim-To-Evidence Checklist"):
        assert phrase in text
    assert "No cutoff fitting" in text and "positive-infinity margin" in text
    assert "not available" in text


def test_fixture_figure_contains_pixels_and_handles_empty_finite_cohorts(documents, tmp_path):
    groups = presentation.snapshot(*documents)
    presentation.render(groups, tmp_path)
    with Image.open(tmp_path / "accepted_event_support.png") as image:
        assert image.size == (2560, 1760)
        assert image.convert("RGB").getextrema() != ((255, 255), (255, 255), (255, 255))
    assert (tmp_path / "accepted_event_support.pdf").stat().st_size > 1000
    svg = (tmp_path / "accepted_event_support.svg").read_text()
    assert "These are not confidence intervals" in svg and "Unlabeled does not mean correct" in svg


def test_existing_namespace_refused_before_input_read(tmp_path):
    with pytest.raises(ValueError, match="Fresh"):
        presentation.generate(tmp_path / "missing", tmp_path / "missing2", tmp_path)


def test_unanchored_input_rejected_before_output_creation(tmp_path):
    a, b, output = tmp_path / "a.json", tmp_path / "b.json", tmp_path / "output"
    a.write_text("{}"); b.write_text("{}")
    with pytest.raises(ValueError, match="anchor"):
        presentation.generate(a, b, output)
    assert not output.exists()


def test_cli_isolated_from_project_import_paths():
    result = subprocess.run([sys.executable, "-I", "-B", presentation.__file__, "--help"], capture_output=True, text=True)
    assert result.returncode == 0 and "--output" in result.stdout


def test_existing_main_and_rc5_index_bytes_are_preserved():
    pinned = {"PUBLICATION_MAIN_TEXT_20261006_v2.md": "2804448cbd99ffa3164fcb9e997640f966103a827c9da761bfb7c3732d79a9c1",
        "publication_package_rc5_index_20261006.json": "66bc02d7821d1097808a870c38198bded9c88c379c24052b1828867af8692922"}
    for name, digest in pinned.items():
        assert hashlib.sha256((presentation.RESULTS / name).read_bytes()).hexdigest() == digest
