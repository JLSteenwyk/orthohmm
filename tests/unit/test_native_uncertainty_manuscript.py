"""Bind the native profile interpretation to its actual family-level evidence."""

import json
from pathlib import Path


BASE = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


def binding():
    return json.loads((BASE / "native_factorial_uncertainty_binding_20261005_v3.json").read_text())


def paragraph():
    text = (BASE / "PUBLICATION_MANUSCRIPT_DRAFT_20260916.md").read_text()
    return text.split("A separate [full-native uncertainty binding]", 1)[1].split(
        "The original-release QfO factorial", 1)[0]


def test_profile_numbers_and_interpretation_match_bound_contrast():
    data = binding()
    contrast = next(row for row in data["contrasts"]
                    if row["on"] == "p1_c0_r1" and row["off"] == "p0_c0_r1")
    metric = contrast["metrics"]["f_score"]
    text = paragraph()
    assert contrast["status"] == "native_records_matched"
    assert f'{metric["difference_percentage_points"]:+.6f}' in text
    lo, hi = metric["bonferroni_percentile_ci"]
    assert f"[{lo:.6f},{hi:.6f}]" in text
    assert lo < 0 < hi
    assert (contrast["family_f1_wins"], contrast["family_f1_ties"],
            contrast["family_f1_losses"]) == (4, 61, 5)
    assert "four family wins,61 ties and five losses" in text
    assert "interval includes zero" in text


def test_native_scope_and_missing_contrasts_remain_explicit():
    data, text = binding(), paragraph()
    assert len(data["bound_cells"]) == 5
    assert sum(row["status"] == "native_records_matched" for row in data["contrasts"]) == 5
    assert sum(row["metrics"] is None for row in data["contrasts"]) == 7
    assert data["new_bootstrap_draws"] == 0
    assert data["replicates_reused"] == 20000
    assert data["multiplicity_endpoints"] == 36
    for literal in ("All350", "five supported contrasts", "Seven native contrasts",
                    "not imputed", "initial HMM search", "conditional downstream",
                    "resolved-pair truth", "independent", "failed timing",
                    "does not replace the complete retained eight-cell analysis"):
        assert literal in text


def test_claim_checklist_retains_neutral_effect_and_source_links():
    claims = (BASE / "PUBLICATION_CLAIMS_20260916.md").read_text()
    row = next(line for line in claims.splitlines()
               if "Fresh native OrthoBench profile expansion" in line)
    assert "Not established" in row and "spanning zero" in row
    assert "not independent confirmation" in row and "initial HMM search" in row
    for filename in ("NATIVE_FACTORIAL_UNCERTAINTY_BINDING_20261005.md",
                     "native_factorial_uncertainty_binding_20261005_v3.json"):
        assert filename in row and (BASE / filename).is_file()


def test_independent_readback_matches_new_bound_result():
    result = json.loads((BASE / "native_factorial_uncertainty_readback_20261005.json").read_text())
    assert result["status"] == "passed" and result["checked_helpers"] == 920
    assert len(result["cells"]) == 5 and len(result["contrasts"]) == 5
    assert all(row["family_records_exact"] for row in result["cells"])
    assert result["binding"]["path"].endswith("native_factorial_uncertainty_binding_20261005_v3.json")
