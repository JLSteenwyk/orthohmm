import json
from pathlib import Path


BASE = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


def text():
    return " ".join((BASE / "PUBLICATION_MAIN_TEXT_20260927.md").read_text().split())


def test_parameter_prose_matches_complete_admitted_summary():
    result = json.loads((BASE / "qfo_private_cpm_parameter_result_22395.json").read_text())
    main = text()
    assert result["complete_panel"] is True
    assert result["endpoint_count"] == result["multiplicity_endpoints"] == 18
    assert result["replicates"] == 100000 and result["seed"] == 20260925
    assert len(result["point_estimates"]) == 7 and len(result["families"]) == 18
    row = next(c for c in result["comparisons"] if c["candidate"] == "cpm_high")["metrics"]["F1"]
    for arm in ("control", "cpm_high"):
        assert f'{100 * result["point_estimates"][arm]["F1"]:.6f}%' in main
    assert f'{100 * row["difference"]:.6f}' in main
    lo, hi = row["bonferroni_percentile_ci"]
    assert f"[{100 * lo:.6f}, {100 * hi:.6f}]" in main
    endpoints = [m for c in result["comparisons"] for m in c["metrics"].values()]
    assert len(endpoints) == 18
    assert all(m["bonferroni_percentile_ci"][0] <= 0 <= m["bonferroni_percentile_ci"][1]
               for m in endpoints)
    assert "All 18 adjusted intervals include zero" in main
    assert "100,000 shared family-bootstrap draws" in main and "20260925" in main
    assert "not averaged across family F1 values" in main
    assert "No defaults were changed" in main


def test_recovered_result_does_not_reclassify_failed_attempts():
    main = text()
    assert "high-CPM experiments remain unavailable" not in main
    assert "Historical SIGSEGV, allocator and admission failures remain failed" in main
    assert "validated content-equivalent private control" in main
    assert "nor proves their cause or memory safety" in main
    assert "no missing score was imputed" in main
    assert "not an OrthoFinder superiority test" in main
    for target in ("QFO_PARAMETER_NEIGHBORHOOD_PROTOCOL_20260919.md",
                   "QFO_COMPLETE_PARAMETER_UNCERTAINTY_20261001.md",
                   "QFO_PRIVATE_CPM_SCORE_RESULT_22394.md",
                   "qfo_parameter_complete_export_20261001/qfo_parameter_neighborhood.pdf"):
        assert target in main and (BASE / target).is_file()


def test_fas_diagnostic_is_not_a_corrected_endpoint():
    result = json.loads((BASE / "qfo_fas_stratum_weights_20261001.json").read_text())
    main = text()
    rows = {r["method"]: r for r in result["methods"]}
    assert len(rows) == 8
    assert result["benchmark_scores_changed"] is result["uncertainty_admitted"] is False
    for method in ("orthomcl_1_4", "orthohmm_phylogeny_satellite_v2"):
        assert f'+{rows[method]["native_minus_diagnostic"]:.9f}' in main
    assert "not corrected benchmark scores or estimates of population bias" in main
    assert "Omissions can also change the stratum means" in main
    assert "No FAS confidence interval or new method ranking follows" in main
    assert "QFO_FAS_STRATUM_WEIGHT_RESULT_20261001.md" in main


def test_numerical_restoration_scope_matches_actual_receipt():
    result = json.loads((BASE / "qfo_parameter_component_relocation_20261001.json").read_text())
    main = text()
    numerical = result["numerical_reproduction"]["report"]
    guard = result["original_path_guard"]
    assert numerical["endpoints"] == 18 and numerical["absolute_tolerance"] == 1e-12
    assert guard["canary_rejected"] is True
    assert guard["original_path_events_after_canary"] == 0 and not guard["loaded_project_modules"]
    assert "all 18 endpoints within 1e-12 using unchanged arithmetic" in main
    assert "not OS containment" in main
    assert "raw-reference recount, native inference restoration" in main
    assert "QFO_PARAMETER_NUMERICAL_COMPONENT_RESULT_20261001.md" in main
    assert "Not submission-ready" in main
