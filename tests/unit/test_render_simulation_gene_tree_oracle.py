import pytest

from benchmark_tools.render_simulation_gene_tree_oracle import render


def test_renderer_rejects_incomplete_results():
    with pytest.raises(ValueError, match="Incomplete"):
        render({"status": "unverified"})


def test_tables_are_computed_from_counts_and_means():
    row = {"condition": "divergent", "scored_cells": 10,
        "mean_metrics": {a: {"f1": v} for a, v in (("inferred", .70), ("generating_root", .71), ("generating_rerooted", .705))},
        "oracle_total_fn": 1000, "true_pairs_across_candidates": 998,
        "candidate_counts": {"oracle_eligible": 3},
        "residual_by_arm": {"generating_root": {"oracle_eligible": {"fn": 2}}}}
    report = {"status": "independent_count_and_candidate_decomposition_verified", "cells": [{}] * 70,
        "summary": [row] * 7, "topology": {"eligible": 21, "unrooted_disagreement": 2,
            "root_only_disagreement": 1, "rooted_disagreement": 3},
        "input_identities_rechecked": 45, "detailed_report": {"path": "local.json", "sha256": "abc"}}
    text = render(report)
    assert "70.000 | 71.000 | 70.500 | +1.000" in text
    assert "1,000 | 998 | 99.800% | 2" in text
    assert "local.json" in text and "`abc`" in text
