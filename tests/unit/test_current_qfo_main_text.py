import gzip
import json
from pathlib import Path


BASE = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


def text():
    return " ".join((BASE / "PUBLICATION_MAIN_TEXT_20260927.md").read_text().split())


def test_null_diagnostic_prose_matches_retained_counts_and_design():
    report = json.loads(gzip.decompress((BASE / "frozen_null_score_observations_20261002.json.gz").read_bytes()))
    main = text()
    assert report["independent_pairs"] == 90000 and report["tail_endpoints"] == 90
    cells = {(s["regime"], s["length"]): s["bands"]["64"][-1] for s in report["summaries"]}
    assert [cells["blosum_background", length]["hits"] for length in (50, 150, 400)] == [0, 1, 0]
    for length in (50, 150, 400):
        assert f'{100 * cells["half_glutamine", length]["fraction"]:.2f}%' in main or cells["half_glutamine", length]["fraction"] == 1
    low, high = cells["half_glutamine", 50]["bonferroni_clopper_pearson"]
    assert f"[{100 * low:.4f}%, {100 * high:.4f}%]" in main
    assert "90,000 independent pairs, not 180,000 independent observations" in main
    assert "Five fixed cutoffs yielded 90 tail endpoints" in main
    assert "0/1/0 of 10,000" in main and "89.83%, 100% and 100%" in main
    assert "180 sparse reference-Python scores, not all native scores" in main


def test_null_diagnostic_interpretation_stays_bounded_in_both_manuscripts():
    main = text()
    extended = " ".join((BASE / "PUBLICATION_MANUSCRIPT_DRAFT_20260916.md").read_text().split())
    for document in (main, extended):
        assert "not real-data orthology false-positive rates" in document or "does not measure real-data orthology false-positive rates" in document
        assert "do not" in document and "rare-tail calibration" in document
        assert "not a measured real-proteome prevalence" in document or "not a measured proteome prevalence" in document
        assert "No coefficients, thresholds or defaults were fitted or promoted" in document or "No new constants, bands, thresholds or defaults are fitted or promoted" in document
        assert "not all native scores" in document or "not all-native-score equivalence" in document
    assert "not a predicted ortholog or an observed pipeline false positive" in main
    assert "No prefilter or biological inference ran" in main
    assert "not reference-family or orthology uncertainty" in main
    for name in ("FROZEN_NULL_SCORE_PROTOCOL_20261002.md", "FROZEN_NULL_SCORE_RESULT_20261002.md",
                 "figures_frozen_null_scores_20261002_v2/frozen_null_scores.pdf",
                 "frozen_null_score_observations_20261002.json.gz"):
        assert name in main and (BASE / name).is_file()


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


def test_descriptive_restoration_prose_matches_retained_component():
    result = json.loads((BASE / "swiss_descriptive_component_20261002.json").read_text())
    main = text()
    restored = result["restoration"]
    assert restored["returncode"] == 0 and restored["python_version"] == "3.12.3"
    assert restored["rows"] == 208 and restored["logical_cells"] == 984
    assert restored["canary_rejected"] and restored["original_path_events_after_canary"] == 0
    assert result["raw_source_admission"] is False
    assert result["bootstrap_intervals_recomputed"] is False
    assert "all 208 rows and 984 score/difference cells" in main
    assert "not raw-source admission, new benchmark scoring or bootstrap uncertainty" in main
    target = "SWISS_DESCRIPTIVE_COMPONENT_20261002.md"
    assert target in main and (BASE / target).is_file()


def test_raw_archive_prose_retains_exact_restoration_scope():
    result = json.loads((BASE / "swiss_raw_archive_restoration_20261002.json").read_text())
    main = text()
    assert result["original_record_occurrences"] == dict(duplication=11, fragment=1774)
    assert [result["archives"][kind]["members"] for kind in ("duplication", "fragment")] == [12, 1774]
    runs = result["copied_verification"]["runs"]
    assert runs["restore"]["python"].startswith("3.12.3")
    assert runs["pytest"]["junit"] == dict(tests=27, errors=0, failures=0, skipped=0)
    for run in runs.values():
        assert run["canaries_blocked"] == 3 and run["later_forbidden_opens"] == 0
        assert run["child_subprocesses_forbidden"] is True
    assert "all 11 and 1,774 original record occurrences" in main
    assert "All 27 affected regression/options cases passed" in main
    assert "including four raw-export regressions with unchanged assertions" in main
    assert "Independently pinned archive and binding digests" in main
    assert "Python-event guards are not OS containment" in main
    assert "Native annotation extraction/admission was not rerun" in main
    assert "private archives remain unuploaded with redistribution uncleared" in main
    assert result["native_annotation_admission_rerun"] is False
    assert result["raw_data_committed_or_uploaded"] is False
    assert result["redistribution_authorized"] is False
    assert result["publication_ready"] is False
    target = "SWISS_RAW_ARCHIVE_RESTORATION_20261002.md"
    assert target in main and (BASE / target).is_file()


def test_completed_components_do_not_claim_complete_release():
    main = text()
    assert "Neither component establishes complete executable study restoration" in main
    assert "all-method cross-host portability" in main
    assert "No submission-ready release or archival DOI is claimed" in main
    assert "older versions are not rendered copies of revised Markdown" in main
    assert "replacement controlled Threadripper timing panel has not run" in main


def test_simulation_completion_and_unavailable_contrasts_remain_separate():
    fixed = json.loads((BASE / "simulation_fixed_native_results_20260916.json").read_bytes())
    variable = json.loads((BASE / "simulation_variable_native_results_20260916.json").read_bytes())
    main = text()
    for report, expected in ((fixed, [70, 64, 0, 0]), (variable, [70, 67, 65, 65])):
        methods = ("orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full", "orthofinder_sequence_only")
        assert [sum(r["method"] == m and r["status"] == "complete" for r in report["records"]) for m in methods] == expected
    contrasts = [c for b in fixed["conditions"].values() for c in b["contrasts"].values()]
    assert len(contrasts) == 14 and all(c["status"] == "no_complete_pairs" for c in contrasts)
    assert "70 and 64 admitted datasets out of 70" in main
    assert "all 14 planned comparisons are unavailable, not wins for OrthoHMM" in main
    assert "admitted 70 and 67 OrthoHMM datasets and 65 full OrthoFinder datasets" in main
    assert "success-conditioned comparisons, not failure-adjusted population estimates" in main


def test_simulation_paired_effects_and_adjustment_match_retained_report():
    report = json.loads((BASE / "simulation_variable_native_results_20260916.json").read_bytes())
    main = text()
    assert report["bootstrap"]["replicates"] == 20000
    assert report["bootstrap"]["f1_multiplicity_count"] == 14
    for method, negatives in (("orthohmm_high_sensitivity", 7), ("orthohmm_satellite_v2", 4)):
        effects = [b["contrasts"][method]["metrics"]["f1"] for b in report["conditions"].values()]
        assert all(e["difference_percentage_points"] < 0 for e in effects)
        assert sum(e["bonferroni_14_ci"][1] < 0 for e in effects) == negatives
    for condition, seeds in (("divergent", 5), ("divergent_turnover", 8)):
        contrast = report["conditions"][condition]["contrasts"]["orthohmm_satellite_v2"]
        assert len(contrast["included_seeds"]) == seeds
        assert f'{contrast["metrics"]["f1"]["difference_percentage_points"]:.4f}' in main
    assert "all seven high-sensitivity contrasts and four phylogenetic contrasts" in main
    assert "available-case means from different seed sets are not paired effects" in main


def test_tree_outcomes_and_adjusted_endpoint_classification_match_summary():
    report = json.loads((BASE / "simulation_tree_robustness_summary_20260917.json").read_bytes())
    main = text()
    assert len(report["records"]) == 560
    assert sum(r["status"] == "complete" for r in report["records"]) == 537
    assert sum(r["status"] == "failed" for r in report["records"]) == 23
    assert report["bootstrap"]["replicates"] == 20000
    assert report["bootstrap"]["seed"] == 20260918 and report["bootstrap"]["multiplicity"] == 126
    oracle = [metric for c in report["contrasts"] if c["target"] == "generating" for metric in c["metrics"].values()]
    assert all(m["bonferroni_126_ci"][0] <= 0 <= m["bonferroni_126_ci"][1] for m in oracle)
    negative = [(c["method"], c["condition"], name) for c in report["contrasts"] for name, metric in c["metrics"].items()
                if metric["bonferroni_126_ci"][1] < 0]
    assert len(negative) == 12 and {name for _, _, name in negative} == {"f1", "recall"}
    assert len({condition for method, condition, _ in negative if method == "orthohmm_satellite_v2"}) == 4
    assert len({condition for method, condition, _ in negative if method == "orthofinder_full"}) == 2
    assert "560 arm outcomes: 537 scored and 23 failed" in main
    assert "No generating-minus-inferred adjusted interval excluded zero" in main
    assert "four OrthoHMM conditions and two OrthoFinder conditions: 12 endpoints" in main
    assert "no precision interval excluded zero" in main


def test_tree_upstream_artifact_boundary_and_oracle_scope_retained():
    audit = json.loads((BASE / "simulation_tree_artifacts_verified_20260917.json").read_bytes())
    main = text()
    assert [sum(c["status"] == s for c in audit["contrasts"]) for s in
            ("retained_upstream_equivalent", "retained_upstream_different", "unavailable")] == [400, 2, 18]
    assert "400 comparisons, differed in two and were unavailable in 18" in main
    assert "prevent a strict tree-only causal attribution there" in main
    assert "oracle diagnostics, not achievable end-to-end inference" in main
    assert "not establish robustness to arbitrary trees or empirical posterior uncertainty" in main


def test_simulation_links_and_attribution_resolve_existing_bibliography():
    main = text()
    entries = json.loads((BASE / "publication_bibliography_20260920_v5.csl.json").read_bytes())
    ids = {entry["id"] for entry in entries}
    assert {"zombi2019online", "pyvolve2015"} <= ids
    assert "[@zombi2019online; @pyvolve2015]" in main
    for target in ("SIMULATION_VARIABLE_NATIVE_INTERPRETATION_20260916.md",
                   "SIMULATION_FIXED_NATIVE_INTERPRETATION_20260916.md",
                   "SIMULATION_TREE_CONTROL_PROTOCOL_20260917.md",
                   "SIMULATION_TREE_ROBUSTNESS_RESULTS_20260917.md", "../SIMULATION_ARITHMETIC_REPLAY.md"):
        assert target in main and (BASE / target).is_file()
    assert "native admission and tree-control inference are not rerun" in main
