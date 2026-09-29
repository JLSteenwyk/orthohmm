import json
from pathlib import Path


BASE = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


def test_current_manuscript_intervals_match_complete_result():
    result = json.loads((BASE / "ob_complete_uncertainty_20260928.json").read_text())
    main = (BASE / "PUBLICATION_MAIN_TEXT_20260927.md").read_text()
    extended = (BASE / "PUBLICATION_MANUSCRIPT_DRAFT_20260916.md").read_text()
    claims = (BASE / "PUBLICATION_CLAIMS_20260916.md").read_text()
    for method, documents in (
            ("orthohmm_phylogeny_satellite_v2", [main, extended]),
            ("orthohmm_high_sensitivity", [extended])):
        for metric, row in result["comparisons"][method]["metrics"].items():
            if method.endswith("high_sensitivity") and metric == "f_score":
                continue
            lo, hi = row["bonferroni_percentile_ci"]
            interval = f"[{lo:.3f}, {hi:.3f}]"
            assert all(interval in document for document in documents)
    lo, hi = result["comparisons"]["orthohmm_phylogeny_satellite_v2"]["metrics"]["f_score"]["bonferroni_percentile_ci"]
    assert f"[{lo:.3f}, {hi:.3f}]" in claims
    assert "Those narrower historical intervals must not be substituted" in extended
    assert "100,000 paired RefOG" in main
    assert "family exchangeability" in main


def test_vgnc_sensitivity_is_not_presented_as_confidence_interval():
    result = json.loads((BASE / "corrected_vgnc_influence_20260928.json").read_text())
    contrast = result["comparisons"]["orthohmm_phylogeny_satellite_v2"]
    interval = f'[{100 * contrast["minimum_deleted_difference"]:.4f}, {100 * contrast["maximum_deleted_difference"]:.4f}]'
    for name in ("PUBLICATION_MAIN_TEXT_20260927.md", "PUBLICATION_MANUSCRIPT_DRAFT_20260916.md",
                 "PUBLICATION_CLAIMS_20260916.md"):
        text = (BASE / name).read_text()
        assert interval in text
        assert "CORRECTED_VGNC_INFLUENCE_RESULT_20260928.md" in text
        assert "not confidence intervals" in text.lower()
    assert result["rows"] == 8 * 16844
    assert all(c["negative"] == 16844 and c["positive"] == c["zero"] == 0
               for c in result["comparisons"].values())


def test_complete_interval_figure_is_linked_in_both_manuscripts():
    target = "figures_ob_complete_uncertainty_20260928/ob_complete_uncertainty.pdf"
    assert (BASE / target).is_file()
    for name in ("PUBLICATION_MAIN_TEXT_20260927.md", "PUBLICATION_MANUSCRIPT_DRAFT_20260916.md"):
        assert target in (BASE / name).read_text()
