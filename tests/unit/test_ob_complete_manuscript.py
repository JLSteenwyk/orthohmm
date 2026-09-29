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
