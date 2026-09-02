import pytest

from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.resume_phylogeny_benchmark import (
    CANDIDATE_RELATIVE_PATH,
    combine_run_metrics,
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
