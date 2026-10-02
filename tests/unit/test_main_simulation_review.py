import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "benchmark_tools/results"


def test_closed_main_simulation_review_artifacts_remain_exact():
    review = json.loads((BASE / "publication_main_simulation_visual_review_20261002.json").read_bytes())
    artifacts = review["closed_review_artifacts"]
    assert len(artifacts) == 14 and len({a["path"] for a in artifacts}) == 14
    for artifact in artifacts:
        path = ROOT / artifact["path"]
        assert path.resolve().is_relative_to(ROOT)
        raw = path.read_bytes()
        assert len(raw) == artifact["bytes"]
        assert hashlib.sha256(raw).hexdigest() == artifact["sha256"]
    pdf = next(a for a in artifacts if a["path"].endswith("document.pdf"))
    assert pdf["bytes"] == 160993
    assert pdf["sha256"] == "65e271daee8bd8047543d614d937d26f285a212ab8e3823eea06a5a6fce87b53"


def test_main_simulation_review_is_a_dated_visual_snapshot_not_release():
    review = json.loads((BASE / "publication_main_simulation_visual_review_20261002.json").read_bytes())
    assert review["status"] == "main_simulation_review_visually_checked"
    assert review["page_count"] == 9 and review["all_nine_pages_inspected"] is True
    assert review["bounds_violations"] == [] and review["observed_clipping_or_overlap"] is False
    assert review["direct_input_hashes_and_git_blobs_checked"] is True
    assert len(review["checked_committed_direct_inputs"]) == 50
    assert all(r["revision"] == review["render_time_source_commit"] for r in review["checked_committed_direct_inputs"])
    assert len(review["citation_ids"]) == 18
    assert review["local_occurrences"] == 48 and review["unique_targets"] == 46
    assert review["untracked_targets"] == []
    for key in ("native_inference_or_scoring_rerun", "tree_bootstrap_rerun", "benchmark_scores_or_defaults_changed",
                "controlled_timing_executed", "new_study_archive_built", "public_release_or_deposition_executed", "publication_ready"):
        assert review[key] is False
    assert review["figures_linked_not_embedded"] is True
    assert review["other_linked_figures_newly_revalidated"] is False
