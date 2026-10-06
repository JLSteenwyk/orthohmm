"""Check the dated native-evidence main text against actual admitted summaries."""

import hashlib
import json
from pathlib import Path
import subprocess

import pytest

from benchmark_tools.render_manuscript_review import citation_ids, local_assets


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "benchmark_tools/results"
MAIN = BASE / "PUBLICATION_MAIN_TEXT_20261006.md"


def read(name):
    return json.loads((BASE / name).read_text())


def text():
    return " ".join(MAIN.read_text().split())


def section(document, heading, next_heading):
    return document.split(heading, 1)[1].split(next_heading, 1)[0]


def test_original_main_and_frozen_goal_remain_byte_identical():
    for path, sha in (
        (BASE / "PUBLICATION_MAIN_TEXT_20261004_v3.md",
         "18baad343773fbf3c142a64b21a5d0c4245c1f37b9ddec7610af3999c20bd5ec"),
        (ROOT / "benchmark_tools/PUBLICATION_GOAL_20261003.txt",
         "7d99ecb39a740b689101e885ca9a8e8d337aaa51d2aa78efb5295e4de27acde0"),
    ):
        assert hashlib.sha256(path.read_bytes()).hexdigest() == sha


def test_new_native_evidence_does_not_remove_original_scientific_sections():
    old = (BASE / "PUBLICATION_MAIN_TEXT_20261004_v3.md").read_text()
    new = MAIN.read_text()
    start = "### Synthetic Null Tails Depend On Composition"
    end = "## Reproducibility And Availability"
    assert section(old, start, end) == section(new, start, end)
    start = "### Accuracy Depends On The Endpoint"
    end = "### Synthetic Null Tails Depend On Composition"
    old_results = section(old, start, end)
    assert new.split(start, 1)[1].startswith(old_results)
    # Preserve the broad scientific methods before adding the native protocol.
    assert section(old, "## Methods", "### Shared-Host Resource Measurement") == section(
        new, "## Methods", "### Fresh Native Ablation Protocol")


@pytest.mark.parametrize("cell", (
    "p0_c0_r0", "p0_c0_r1", "p0_c1_r0", "p0_c1_r1", "p1_c0_r1", "p1_c1_r0",
))
def test_all_six_native_orthobench_rows_match_snapshot(cell):
    rows = read("native_factorial_progress_20261005_v6/report.json")["rows"]
    row = next(r for r in rows if r["dataset"] == "orthobench" and r["cell"] == cell)
    label = cell.upper().replace("_", "/")
    values = [f'{100 * row[k]:.4f}' for k in ("f1", "precision", "recall")]
    assert "| " + label + " | " + " | ".join(values) + " |" in MAIN.read_text()


def test_native_ob_contrasts_retain_adjustment_and_unavailable_scope():
    data = read("native_factorial_uncertainty_binding_20261005_v4.json")
    assert len(data["bound_cells"]) == 6 and data["planned_contrasts"] == 12
    assert data["multiplicity_endpoints"] == 36 and data["replicates_reused"] == 20000
    assert data["new_bootstrap_draws"] == 0
    matched = [c for c in data["contrasts"] if c["metrics"] is not None]
    assert len(matched) == 6
    main = text()
    for row in matched:
        metric = row["metrics"]["f_score"]
        value = metric["difference_percentage_points"]
        assert f'{value:+.4f}' in main
        low, high = metric["bonferroni_percentile_ci"]
        assert f'[{low:.4f}, {high:.4f}]' in main
        if row["factor"] == "profile_expansion":
            assert low <= 0 <= high
            assert '/'.join(str(row[k]) for k in (
                "family_f1_wins", "family_f1_ties", "family_f1_losses")) in main
    assert "Six of 12 planned simple effects" in main
    assert "neither complete a fresh eight-cell factorial" in main
    assert "not the total HMM contribution" in main


@pytest.mark.parametrize("endpoint", ("VGNC", "SwissTrees", "TreeFam-A", "GO", "EC", "FAS"))
def test_both_native_qfo_rows_match_actual_admitted_snapshot(endpoint):
    rows = read("native_qfo_scientific_scores_20261006_v1/report.json")["rows"]
    rows = [r for r in rows if r["accuracy_admitted"]]
    assert [r["cell"] for r in rows] == ["p0_c0_r0", "p0_c0_r1"]
    label = "F1" if endpoint in ("VGNC", "SwissTrees", "TreeFam-A") else (
        "Sample mean" if endpoint == "FAS" else "Similarity")
    expected = f'| {endpoint} | {label} | ' + " | ".join(
        f'{r["scores"][endpoint]:.8f}' for r in rows) + " |"
    assert expected in MAIN.read_text()
    for row in rows:
        assert f'{row["secondary_mean"]:.8f}' in text()
        assert f'{100 * row["relation_coverage"]:.4f}%' in text()


def test_native_swiss_effects_are_not_rounded_endpoint_subtraction_or_new_draws():
    data = read("native_qfo_swiss_uncertainty_binding_22449_20261006.json")
    assert len(data["bound_cells"]) == 2 and data["new_bootstrap_draws"] == 0
    assert data["multiplicity_endpoints"] == 42
    matched = [r for r in data["contrasts"] if r["metrics"] is not None]
    assert len(matched) == 1 and matched[0]["name"] == "R_at_P0_C0"
    main = text()
    for name, metric in matched[0]["metrics"].items():
        assert f'{100 * metric["difference"]:+.4f}' in main
        low, high = metric["bonferroni_percentile_ci"]
        assert f'[{100 * low:.4f}, {100 * high:.4f}]' in main
        if name == "F1":
            assert low <= 0 <= high
            assert '/'.join(str(metric[k]) for k in (
                "family_wins", "family_ties", "family_losses")) in main
    for phrase in ("other 13 planned contrasts remain unavailable",
                   "count-based arithmetic differs slightly from native serialized endpoints",
                   "F1 interval includes zero", "initial HMM search on"):
        assert phrase in main


def test_functional_pair_composition_preserves_original_means_and_dependence():
    data = read("native_qfo_functional_pair_composition_20261006_v1.json")
    main = text()
    for row in data["comparisons"]:
        result = row["result"]
        if row["metric"] in ("GO", "EC"):
            assert result["left_only_pairs"] == 0
            assert result["left_pairs"] == result["shared_pairs"]
            assert result["shared_pairs_with_different_serialized_scores"] == 0
            for key in ("left_pairs", "right_pairs", "right_only_pairs"):
                assert f'{result[key]:,}' in main
            excluded_mean = (result["right_score_sum_millionths"] -
                             result["shared_right_sum_millionths"]) / (
                                 1e6 * result["right_only_pairs"])
            assert f'{excluded_mean:.10f}' in main
        else:
            assert result["shared_sample_pairs"] == 1007
            for key in ("left_sample_pairs", "right_sample_pairs", "shared_sample_pairs"):
                assert f'{result[key]:,}' in main
            for key in ("shared_fraction_of_left", "shared_fraction_of_right"):
                assert f'{100 * result[key]:.4f}%' in main
            assert f'{result["original_sample_mean_difference"]:+.10f}' in main
    for phrase in ("Intersection-only means would change the endpoint",
                   "not changed similarities on common pairs",
                   "Unseeded method-specific sample mixtures",
                   "paired uncertainty remains unfinished"):
        assert phrase in main
    assert data["new_bootstrap_draws"] == 0 and data["uncertainty_admitted"] is False


def test_failure_scope_and_current_archive_limitations_are_explicit():
    data = read("native_qfo_scientific_scores_20261006_v1/report.json")
    recovered = next(r for r in data["rows"] if r["cell"] == "p0_c0_r1")
    assert recovered["resources"] is None
    assert recovered["timing_admitted"] is recovered["timing_eligible"] is False
    for phrase in ("Not submission-ready", "five remain unavailable", "failed 22437 timing",
                   "not successful timing", "not the selected-default all-tool comparison",
                   "unknown and potentially tool-dependent impact",
                   "not estimates of isolated performance", "safe capacity and valid accounting",
                   "neither archive includes this 6 October native-evidence revision",
                   "No submission-ready release or archival DOI is claimed"):
        assert phrase in text()


def test_all_local_links_and_original_eighteen_citations_resolve():
    parsed = subprocess.check_output(["pandoc", "--from=markdown", "--to=json", str(MAIN)], text=True)
    document = json.loads(parsed)
    _, targets = local_assets(document, MAIN, ROOT)
    assert len(citation_ids(document)) == 18
    old = json.loads(subprocess.check_output([
        "pandoc", "--from=markdown", "--to=json",
        str(BASE / "PUBLICATION_MAIN_TEXT_20261004_v3.md")], text=True))
    assert citation_ids(document) == citation_ids(old)
    _, old_targets = local_assets(old, BASE / "PUBLICATION_MAIN_TEXT_20261004_v3.md", ROOT)
    assert old_targets.keys() <= targets.keys()
    for name in (
        "native_factorial_receipt_amendment_20261004/plan.json",
        "native_factorial_progress_20261005_v6/report.json",
        "native_factorial_uncertainty_binding_20261005_v4.json",
        "native_factorial_uncertainty_readback_20261005_v2.json",
        "native_qfo_scientific_scores_20261006_v1/report.json",
        "native_qfo_swiss_uncertainty_binding_22449_20261006.json",
        "recovered_native_qfo_swiss_readback_22449_20261006.json",
        "native_qfo_p0c0_figure_20261006_v1/native_qfo_p0c0.pdf",
        "native_qfo_functional_pair_sql_readback_20261006.json",
    ):
        assert "benchmark_tools/results/" + name in targets
    bibliography = read("publication_bibliography_20260920_v5.csl.json")
    assert citation_ids(document) <= {r["id"] for r in bibliography}
