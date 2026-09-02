import json

import pytest

from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.resume_phylogeny_benchmark import (
    CANDIDATE_RELATIVE_PATH,
    MEMBERSHIP_RELATIVE_PATH,
    combine_run_metrics,
    load_membership_constraints,
    verify_candidate_checkpoint,
)


def test_verified_candidate_checkpoint_requires_manifest_match(tmp_path):
    candidate = tmp_path / "clusters.txt"
    candidate.write_text("a b\n", encoding="utf-8")
    expected = file_provenance(candidate)
    expected["path"] = CANDIDATE_RELATIVE_PATH
    failed = {
        "status": "failed",
        "harness": {"output_manifest": [expected]},
    }

    actual = verify_candidate_checkpoint(failed, candidate)
    assert actual["sha256"] == expected["sha256"]

    candidate.write_text("a c\n", encoding="utf-8")
    with pytest.raises(ValueError, match="does not match"):
        verify_candidate_checkpoint(failed, candidate)


def test_combine_run_metrics_sums_wall_and_preserves_failed_stage():
    failed = {
        "wall_s": 100.0,
        "peak_process_tree_rss_bytes": 300,
        "stages": {"search": {"wall_s": 80}, "phylogeny": {"wall_s": 2}},
        "counts": {"genes": 10},
    }
    resumed = {
        "wall_s": 25.5,
        "peak_process_tree_rss_bytes": 200,
        "stages": {"phylogeny": {"wall_s": 20}},
    }

    combined = combine_run_metrics(failed, resumed, {"root_hogs": 4})

    assert combined["wall_s"] == 125.5
    assert combined["peak_process_tree_rss_bytes"] == 300
    assert combined["stages"]["phylogeny_initial_failure"]["wall_s"] == 2
    assert combined["stages"]["phylogeny_restart"]["wall_s"] == 20
    assert combined["counts"] == {"genes": 10, "root_hogs": 4}


def test_membership_constraints_require_verified_nonempty_groups(tmp_path):
    path = tmp_path / "merges.json"
    path.write_text(json.dumps([
        {"source_genes": ["a"], "target_genes": ["b", "c"]}
    ]))
    expected = file_provenance(path)
    expected["path"] = MEMBERSHIP_RELATIVE_PATH
    failed = {
        "status": "failed",
        "harness": {"output_manifest": [expected]},
    }

    constraints, actual = load_membership_constraints(failed, path)

    assert constraints[0]["source_genes"] == ["a"]
    assert actual["sha256"] == expected["sha256"]

    path.write_text("[]")
    with pytest.raises(ValueError, match="does not match"):
        load_membership_constraints(failed, path)
