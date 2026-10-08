"""Prospective manuscript integration with explicitly synthetic final evidence."""

from copy import deepcopy
import csv
import json
from pathlib import Path

import pytest

from benchmark_tools import prepare_composed_native_main_text as main
from benchmark_tools import review_composed_native_qfo_figure as reviewer
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_composed_native_qfo_figure import ROOT, data, rendered
from tests.unit.test_export_native_qfo_factorial_scores import write


@pytest.fixture(scope="module")
def evidence(tmp_path_factory, rendered):
    root = tmp_path_factory.mktemp("synthetic_composed_manuscript")
    plot = json.loads((rendered / "manifest.json").read_text())
    snapshot = json.loads(Path(plot["snapshot"]["path"]).read_text())
    failed = ROOT / "benchmark_tools/results/native11_qfo_scoring_failure_addendum_20261008_v1/report.json"
    snapshot["failure_addendum"] = record(failed)
    snapshot_ref = write(root / "snapshot.json", snapshot)
    uncertainty = json.loads(Path(plot["native_intervals"]["path"]).read_text())
    uncertainty["snapshot"] = snapshot_ref
    intervals_ref = write(root / "intervals.json", uncertainty)
    plot.update(snapshot=snapshot_ref, native_intervals=intervals_ref)
    plot_ref = write(rendered / "manifest.json", plot)
    reader = reviewer.review(rendered, root / "preview.png")
    reader_ref = write(root / "reader.json", reader)
    parent = ROOT / "benchmark_tools/results/PUBLICATION_MAIN_TEXT_20261007_v4.md"
    refs = dict(parent=record(parent), snapshot=snapshot_ref, intervals=intervals_ref,
        figure=plot_ref, reader=reader_ref, failure=record(failed))
    specifications = {key: (ref["path"], ref["sha256"]) for key, ref in refs.items()}
    docs, refs, _ = main.inputs(specifications)
    return root, docs, refs, specifications


def test_generated_native_scope_and_unrelated_sections_preserved(evidence):
    root, docs, refs, _ = evidence
    text = main.manuscript(docs, refs, root, root / "claims.tsv")
    parent = docs["parent"]
    assert parent[parent.index("## Results\n"):parent.index(main.START)] in text
    assert parent[parent.index(main.MECHANISM):parent.index(main.END)] in text
    assert parent[parent.index(main.END):parent.index("## Reproducibility And Availability\n")] in text
    assert parent[parent.index("## References\n"):] in text
    assert "Five of seven fresh corrected-QfO cells have admitted accuracy" in text
    assert "4 of 14 planned SwissTrees contrasts" in text
    assert "All 42 planned metric endpoints" in text
    assert "Five completed endpoint tasks supply zero admitted scores" in text
    assert "OUT_OF_MEMORY" in text and "59.4612%" in text
    assert main.ABSTRACT_OLD not in text
    assert "All three adjusted F1 intervals include zero." not in text
    assert "not a new mechanistic diagnosis of native12" in text
    assert "new manuscript\nrender/archive receipt" in text


def test_interval_interpretation_comes_from_data_not_parent_prose(evidence):
    root, original, refs, _ = evidence
    docs = deepcopy(original)
    effect = next(row for row in docs["intervals"]["comparisons"] if row["name"] == "R_at_P0_C0")
    effect["metrics"]["F1"].update(paired_percentile_ci=[.03, .17], bonferroni_percentile_ci=[.02, .2])
    text = main.manuscript(docs, refs, root, root / "claims.tsv")
    assert "1 strictly above zero" in text
    assert "3 include zero" in text
    summary = main.interpretation([effect])
    assert summary["positive"] == ["R_at_P0_C0"]
    effect["metrics"]["F1"]["bonferroni_percentile_ci"] = [-.2, -.02]
    assert main.interpretation([effect])["negative"] == ["R_at_P0_C0"]
    effect["metrics"]["F1"]["bonferroni_percentile_ci"] = [0, .02]
    assert main.interpretation([effect])["includes_zero"] == ["R_at_P0_C0"]


@pytest.mark.parametrize("change", ["source", "figure_binding", "reader_binding", "failure_binding", "partial_readback",
    "mean", "failure_score", "failure_peak", "failure_coverage", "failure_allocation", "figure_draws", "reader_schema"])
def test_mixed_or_invalid_evidence_refused(evidence, change):
    _, original, refs, _ = evidence
    docs = deepcopy(original)
    if change == "source": docs["snapshot"]["source"] = {}
    elif change == "figure_binding": docs["figure"]["snapshot"] = {}
    elif change == "reader_binding": docs["reader"]["manifest"] = {}
    elif change == "failure_binding": docs["snapshot"]["failure_addendum"] = {}
    elif change == "partial_readback": docs["reader"]["score_rows"] = 30
    elif change == "mean": docs["snapshot"]["rows"][6]["secondary_mean"] = 0
    elif change == "failure_score": docs["failure"]["row"]["admitted_endpoint_score_count"] = 1
    elif change == "failure_peak": docs["failure"]["row"]["scoring_peak_memory_bytes"] = 0
    elif change == "failure_coverage": docs["failure"]["row"]["submitted_pairs"] = 0
    elif change == "failure_allocation": docs["failure"]["row"]["scoring_memory_limit_bytes"] = 128*1024**3
    elif change == "figure_draws": docs["figure"]["new_bootstrap_draws"] = 100000
    else: docs["reader"]["schema"] = "ordinary"
    with pytest.raises(ValueError): main.validate(docs, refs)


@pytest.mark.parametrize("anchor", [main.START, main.MECHANISM, main.END, main.ABSTRACT_OLD])
def test_missing_or_duplicated_manuscript_anchor_refused(evidence, anchor):
    root, original, refs, _ = evidence
    docs = deepcopy(original)
    for parent in (docs["parent"].replace(anchor, "", 1), docs["parent"]+anchor):
        docs["parent"] = parent
        with pytest.raises(ValueError, match="anchors"): main.manuscript(docs, refs, root, root / "claims.tsv")


def test_actual_frozen_parent_and_failure_pins_required(evidence, tmp_path):
    _, _, _, specifications = evidence
    changed = dict(specifications)
    parent = tmp_path / "changed.md"
    parent.write_text("Changed source")
    ref = record(parent)
    changed["parent"] = (ref["path"], ref["sha256"])
    with pytest.raises(ValueError, match="frozen v4"): main.inputs(changed)
    missing = dict(specifications)
    missing.pop("failure")
    with pytest.raises(ValueError, match="six manuscript inputs"): main.inputs(missing)


def test_generate_fresh_fixture_outputs_and_claims_not_a_production_manuscript(evidence, tmp_path):
    _, _, _, specifications = evidence
    output, receipt, claims = (tmp_path / name for name in ("synthetic_v5.md", "receipt.json", "claims.tsv"))
    result = main.generate(tmp_path, specifications, output, receipt, claims)
    assert result["admitted_endpoint_count"] == 30 and result["unavailable_score_cells"] == 2
    assert result["conditional_interval_rows"] == 12
    assert result["new_bootstrap_draws"] == 0 and result["consumed_native_bootstrap_draws"] == 100000
    assert not result["publication_ready"] and not result["render_review_established"]
    assert not result["historical_source_modified"] and not result["independent_confirmation"]
    with claims.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == 5
    assert all(row["evidence"] and row["excluded"] for row in rows)
    with pytest.raises(ValueError, match="fresh outputs"):
        main.generate(tmp_path, specifications, output, receipt, claims)


def test_existing_aliased_or_outside_outputs_refused_before_inputs(tmp_path, monkeypatch):
    monkeypatch.setattr(main, "inputs", lambda *args: pytest.fail("unexpected input reads"))
    output = tmp_path / "output"
    for receipt, claims in ((output, tmp_path / "claims"),
        (tmp_path / "receipt", tmp_path.parent / "outside")):
        with pytest.raises(ValueError): main.generate(tmp_path, {}, output, receipt, claims)
    output.write_text("retain me")
    with pytest.raises(ValueError): main.generate(tmp_path, {}, output, tmp_path / "receipt", tmp_path / "claims")
    assert output.read_text() == "retain me"


@pytest.mark.parametrize("value", [True, None, float("nan"), "0.5"])
def test_secondary_mean_has_numeric_not_boolean_or_missing_scope(evidence, value):
    _, original, refs, _ = evidence
    docs = deepcopy(original)
    docs["snapshot"]["rows"][6]["secondary_mean"] = value
    with pytest.raises(ValueError): main.validate(docs, refs)
