import json
from pathlib import Path

import pytest

from benchmark_tools import deferred_qfo_inferred as qfo
from benchmark_tools import report_deferred_cells as reporting
from benchmark_tools import run_deferred_inferred_supplement as supplement
from benchmark_tools.run_membership_hillclimb import score_qfo
from benchmark_tools.run_qfo_recovered_assessment import require_darwin_work_path
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def test_short_scoring_namespace_fits_actual_repo_where_old_identifier_fails():
    repo = Path(qfo.__file__).resolve().parent.parent
    output = repo / "benchmarks/results/deferred_qfo_inferred_r1_v2_20261013"
    ids = [qfo.scoring_id(repo, output, i) for i in (0, 1)]
    assert ids[0] != ids[1]
    assert qfo.scoring_id(repo, output, 0) == ids[0]
    assert qfo.scoring_id(repo, output.with_name("other_attempt"), 0) != ids[0]
    for name in ids:
        require_darwin_work_path(repo / "qfo_benchmark/w" / ("hc1_" + name))
    with pytest.raises(ValueError, match="length limit"):
        require_darwin_work_path(repo / "qfo_benchmark/w/hc1_deferred_own_0123456789abcdef_0")


@pytest.mark.parametrize("index", [-1, 2, True, 0., "0", None])
def test_scoring_arm_index_is_strict(index):
    with pytest.raises(ValueError, match="indices"):
        qfo.scoring_id(Path("/repo"), Path("/output"), index)


def test_namespace_rejects_unreasonably_long_root_before_pair_conversion(tmp_path):
    with pytest.raises(ValueError, match="length limit"):
        score_qfo(tmp_path / ("x" * 160), tmp_path / "absent_native", tmp_path, "short")
    assert list(tmp_path.iterdir()) == []


def test_darwin_path_boundary_unchanged():
    suffix = "/xx/xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx/predictions.db"
    allowed = 159 - len(suffix)
    require_darwin_work_path("/" + "x" * (allowed - 1))
    with pytest.raises(ValueError, match="length limit"):
        require_darwin_work_path("/" + "x" * allowed)


def test_orthobench_control_uses_matching_fresh_native_products_and_candidates(tmp_path):
    fresh = tmp_path / "fresh"
    original = ["--output-directory", str(tmp_path / "retained"), "--candidate-clusters", "old_candidate",
                "--membership-constraints", "old_constraints"]
    phylo, candidate, constraints, checkpoint = supplement.native_control_paths("OrthoBench", fresh, original)
    assert phylo == fresh / "orthohmm_phylogeny"
    assert candidate == fresh / "orthohmm_working_res/phylogeny_candidate_superfamilies.txt"
    assert constraints == fresh / "orthohmm_working_res/phylogeny_candidate_merges.json"
    assert checkpoint == fresh
    tk = supplement.native_control_paths("ThreeKingdoms", fresh, original)
    assert tk == (tmp_path / "retained/orthohmm_phylogeny", Path("old_candidate"), Path("old_constraints"), tmp_path / "retained")


@pytest.mark.parametrize("names", [[], ["QfO"], ["OrthoBench", "OrthoBench"], ["ThreeKingdoms", "OrthoBench"]])
def test_invalid_retry_dataset_inventory_is_rejected(names):
    with pytest.raises(ValueError, match="subset"):
        supplement.selected_datasets(names)


def test_orthobench_only_retry_is_prespecified_without_rerunning_three_kingdoms():
    assert supplement.selected_datasets(["OrthoBench"]) == ("OrthoBench",)
    assert supplement.selected_datasets(supplement.DATASETS) == supplement.DATASETS
    plan = {"cases": [{"dataset": "OrthoBench", "arms": [{"label": label} for label in reporting.LABELS]}]}
    assert reporting.planned_case(plan, "OrthoBench", "control")[0]["label"] == "control"


def test_terminal_preparation_failure_is_not_mislabeled_pending(tmp_path):
    root, snapshots = tmp_path / "attempt", tmp_path / "snapshots"
    snapshots.mkdir()
    root.mkdir()
    plan = root / "plan.json"
    plan.write_text(json.dumps({"output": str(root), "source_revision": "fixture"}))
    (root / "results.json").write_text(json.dumps({"execution_plan": record(plan), "source_revision": "fixture"}))
    (root / "OrthoBench").mkdir()
    (root / "OrthoBench/results.json").write_text(json.dumps({"dataset": "OrthoBench", "status": "failed_preparation_or_gate",
        "arms": [], "error_type": "ValueError", "error": "exact control mismatch"}))
    cells, _ = reporting.collect_deferred_cells(tmp_path, "deferred_inferred", "OrthoBench", root, snapshots)
    assert all(c["status"] == "not_run_terminal_producer_gap" and c["scores"] == {} for c in cells)
    assert all(c["error"] == "exact control mismatch" for c in cells)


def test_failed_scoring_reason_is_retained_without_any_claimed_score(tmp_path):
    root, snapshots = tmp_path / "attempt", tmp_path / "snapshots"
    snapshots.mkdir()
    root.mkdir()
    plan = root / "plan.json"
    plan.write_text(json.dumps({"output": str(root), "source_revision": "fixture"}))
    (root / "results.json").write_text(json.dumps({"execution_plan": record(plan), "source_revision": "fixture",
        "status": "finished_with_explicit_gaps", "arms": [{"label": "control", "status": "failed",
        "scores": {"VGNC F": .9}, "error": "Darwin work path exceeds validated length limit", "error_type": "ValueError"}]}))
    cells, _ = reporting.collect_deferred_cells(tmp_path, "deferred_inferred", "QfO", root, snapshots)
    assert cells[0]["scores"] == {} and "Darwin" in cells[0]["error"]
    assert cells[1]["status"] == "not_run_terminal_producer_gap"
