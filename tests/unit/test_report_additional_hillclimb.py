from copy import deepcopy
import json

import pytest

from benchmark_tools import report_additional_hillclimb as reporter
from benchmark_tools.metric_frontier import ENDPOINTS
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_membership_hillclimb import PRODUCTS, write_json


def verified_cell(dataset="QfO"):
    return {"dataset": dataset, "status": "independent_content_readback_verified",
            "scores": {k: .8 for k in ENDPOINTS if k.startswith(dataset + ":")},
            "products": {name: {"path": name, "bytes": 1, "sha256": "same"} for name in PRODUCTS}}


def panel(tmp_path, monkeypatch, status="fresh_original_cap400_running"):
    root = tmp_path / "continuation"
    root.mkdir()
    plan = root / "execution_plan.json"
    write_json(plan, {"repo": str(tmp_path), "output": str(root), "inputs": [], "sources": [],
                      "runtime": {"files": []}})
    write_json(root / "results.json", {"status": status, "execution_plan": record(plan), "job_id": "live",
                                      "scores": {"VGNC F": 1}})
    parent = tmp_path / "parent.json"
    table = tmp_path / "parent.tsv"
    table.write_text("old table\n")
    failed = {"dataset": "QfO", "status": "not_run_failed_control", "scores": {}}
    unrelated = {"key": "other", "label": "other", "scores": {}, "cells": []}
    write_json(parent, {"score_table": record(table), "methods": [unrelated,
        {"key": "full_search:cap400", "label": "cap400", "scores": {}, "cells": [failed]}],
        "round_screening": {}, "limitations": [], "reporter_revision": "earlier_reporter"})
    monkeypatch.setattr(reporter, "runtime", lambda: {"files": []})
    return root, parent, unrelated, failed


def test_pending_producer_scores_are_not_admitted(tmp_path, monkeypatch):
    root, parent, unrelated, _ = panel(tmp_path, monkeypatch)
    monkeypatch.setattr(reporter, "continuation_cell", lambda *args: pytest.fail("Pending cell was read back"))
    result = reporter.build(tmp_path, parent, tmp_path / "out", continuation_root=root)
    assert result["methods"][0] == unrelated
    assert result["methods"][1]["scores"] == {}
    assert "NA" in (tmp_path / "out/scores.tsv").read_text()
    assert "reporter_revision" not in result
    assert result["reporter_files"][0] == record(reporter.__file__)


def test_completed_continuation_keeps_original_unrun_cell(tmp_path, monkeypatch):
    root, parent, unrelated, failed = panel(tmp_path, monkeypatch, "original_cap400_native_and_all_six_scores_verified")
    good = verified_cell()
    monkeypatch.setattr(reporter, "continuation_cell", lambda *args: deepcopy(good))
    result = reporter.build(tmp_path, parent, tmp_path / "out", continuation_root=root)
    assert result["methods"][0] == unrelated
    updated = result["methods"][1]
    assert updated["scores"] == good["scores"]
    assert updated["cells"][0] == failed
    assert updated["cells"][1]["status"] == good["status"]
    assert result["defaults_promoted"] is False and result["goal_completed"] is False
    assert "full_search:cap400" not in result["round_screening"]["shared_hit"]["comparisons"]


def test_native_conflict_is_not_hidden_by_identical_scores(tmp_path, monkeypatch):
    root, parent, _, _ = panel(tmp_path, monkeypatch, "original_cap400_native_and_all_six_scores_verified")
    original = json.loads(parent.read_text())
    good = verified_cell()
    original["methods"][1].update(cells=[good], scores=good["scores"])
    write_json(parent, original)
    changed = deepcopy(good)
    changed["products"][PRODUCTS[0]]["sha256"] = "different"
    monkeypatch.setattr(reporter, "continuation_cell", lambda *args: changed)
    with pytest.raises(ValueError, match="native products"):
        reporter.build(tmp_path, parent, tmp_path / "out", continuation_root=root)


def test_missing_producer_is_not_called_running_or_successful(tmp_path, monkeypatch):
    _, parent, _, _ = panel(tmp_path, monkeypatch)
    result = reporter.build(tmp_path, parent, tmp_path / "out", shared_root=tmp_path / "missing")
    assert result["additional_trial_observations"][0]["status"] == "report_missing_execution_state_not_inferred"
    assert len(result["methods"]) == 2


def recovery_panel(tmp_path, monkeypatch, status):
    from benchmark_tools.recover_canonical_qfo_scoring import REVISION, PATH_ERROR
    root, parent, unrelated, _ = panel(tmp_path, monkeypatch, status)
    row = {"dataset": "QfO", "label": "canonical", "status": "failed", "scores": {},
           "error": PATH_ERROR, "error_type": "ValueError", "outputs": [{"path": "native"}],
           "graph_execution": {"returncode": 0}, "native_execution": {"returncode": 0},
           "native_summary": {"genes": 984137}}
    snapshot = root / "failed_parent_snapshot.json"
    write_json(snapshot, {"source_revision": REVISION, "job_id": "24240", "arms": [row]})
    path = root / "execution_plan.json"
    plan = json.loads(path.read_text())
    plan.update(parent_snapshot=record(snapshot), original_failed_row=row)
    write_json(path, plan)
    report = json.loads((root / "results.json").read_text())
    report["execution_plan"] = record(path)
    write_json(root / "results.json", report)
    return root, parent, unrelated, row


def test_pending_canonical_recovery_keeps_failure_and_no_scores(tmp_path, monkeypatch):
    root, parent, unrelated, row = recovery_panel(tmp_path, monkeypatch, "scoring_only_recovery_running")
    monkeypatch.setattr(reporter, "canonical_recovery_cell", lambda *args: pytest.fail("Pending scoring admitted"))
    result = reporter.build(tmp_path, parent, tmp_path / "out", canonical_recovery_root=root)
    assert result["methods"][0] == unrelated
    method = next(m for m in result["methods"] if m["key"] == "canonical_collection:canonical")
    assert method["scores"] == {}
    assert any(c.get("error") == row["error"] and c["status"] == "failed" for c in method["cells"])
    assert "NA" in (tmp_path / "out/scores.tsv").read_text()


def test_completed_canonical_recovery_retains_original_failure(tmp_path, monkeypatch):
    root, parent, unrelated, row = recovery_panel(tmp_path, monkeypatch,
        "canonical_qfo_scoring_recovered_independently_verified")
    good = verified_cell()
    monkeypatch.setattr(reporter, "canonical_recovery_cell", lambda *args: deepcopy(good))
    result = reporter.build(tmp_path, parent, tmp_path / "out", canonical_recovery_root=root)
    assert result["methods"][0] == unrelated
    method = next(m for m in result["methods"] if m["key"] == "canonical_collection:canonical")
    assert method["scores"] == good["scores"]
    assert method["cells"][0]["status"] == "failed"
    assert method["cells"][0]["error"] == row["error"]
    assert method["cells"][1]["status"] == "independent_content_readback_verified"
    assert result["defaults_promoted"] is False and result["goal_completed"] is False


def test_recovery_cannot_change_archived_failure(tmp_path, monkeypatch):
    root, parent, _, _ = recovery_panel(tmp_path, monkeypatch, "scoring_only_recovery_running")
    path = root / "execution_plan.json"
    plan = json.loads(path.read_text())
    plan["original_failed_row"]["error"] = "different"
    write_json(path, plan)
    report = json.loads((root / "results.json").read_text())
    report["execution_plan"] = record(path)
    write_json(root / "results.json", report)
    with pytest.raises(ValueError, match="original failed row"):
        reporter.build(tmp_path, parent, tmp_path / "out", canonical_recovery_root=root)


def test_shared_result_joins_only_its_verified_canonical_control(tmp_path, monkeypatch):
    root, parent, unrelated, failed = panel(tmp_path, monkeypatch)
    report = json.loads((root / "results.json").read_text())
    report.update(status="running", arms=[{"dataset": "OrthoBench", "label": "shared_hit1",
                                          "status": "independent_content_and_score_readback_verified"}])
    write_json(root / "results.json", report)
    good = verified_cell("OrthoBench")
    control = verified_cell("OrthoBench")
    good["matched_control"] = {"cell": control}
    monkeypatch.setattr(reporter, "shared_cell", lambda *args: deepcopy(good))
    result = reporter.build(tmp_path, parent, tmp_path / "out", shared_root=root)
    assert result["methods"][0] == unrelated
    assert result["methods"][1]["cells"] == [failed]
    new = {m["key"]: m for m in result["methods"][2:]}
    assert set(new) == {"canonical_collection:canonical", "shared_hit:shared_hit1"}
    assert new["canonical_collection:canonical"]["scores"] == control["scores"]
    assert new["shared_hit:shared_hit1"]["scores"] == good["scores"]
    assert all(not m["complete_eight_endpoints"] for m in new.values())


def test_duplicate_shared_dataset_is_rejected_before_readback(tmp_path, monkeypatch):
    root, parent, _, _ = panel(tmp_path, monkeypatch)
    report = json.loads((root / "results.json").read_text())
    report["arms"] = [{"dataset": "OrthoBench"}, {"dataset": "OrthoBench"}]
    write_json(root / "results.json", report)
    with pytest.raises(ValueError, match="dataset inventory"):
        reporter.build(tmp_path, parent, tmp_path / "out", shared_root=root)


@pytest.mark.parametrize("change", ["repo", "output", "runtime"])
def test_changed_producer_namespace_or_runtime_is_rejected(tmp_path, monkeypatch, change):
    root, parent, _, _ = panel(tmp_path, monkeypatch)
    path = root / "execution_plan.json"
    plan = json.loads(path.read_text())
    plan[change] = {"files": [], "changed": True} if change == "runtime" else str(tmp_path / "wrong")
    write_json(path, plan)
    report = json.loads((root / "results.json").read_text())
    report["execution_plan"] = record(path)
    write_json(root / "results.json", report)
    with pytest.raises(ValueError):
        reporter.build(tmp_path, parent, tmp_path / "out", continuation_root=root)


def test_continuation_parameter_mismatch_precedes_native_readback(tmp_path, monkeypatch):
    monkeypatch.setattr(reporter, "verified_native", lambda *args: pytest.fail("Invalid case was audited"))
    plan = {"case": {"fresh_case": {"dataset": "QfO", "label": "cap400", "search_max_candidates_per_query": 100}}}
    with pytest.raises(ValueError, match="cap400 continuation"):
        reporter.continuation_cell(tmp_path, plan, {"execution": {"returncode": 0}})


def shared_fixture(tmp_path, monkeypatch):
    source = tmp_path / "source"
    (source / "benchmark_tools").mkdir(parents=True)
    (source / "orthohmm").mkdir()
    driver = source / "benchmark_tools/benchmark_production.py"
    driver.write_text("# fixture\n")
    control_root = tmp_path / "control"
    control_root.mkdir()
    control_plan = control_root / "execution_plan.json"
    write_json(control_plan, {"source_revision": "control_revision"})
    control_path = control_root / "results.json"
    control = {"dataset": "OrthoBench", "label": "canonical", "status": "independent_content_readback_verified",
               "scores": {"OrthoBench F1": .8}}
    write_json(control_path, {"arms": [control], "execution_plan": record(control_plan),
                              "source_revision": "control_revision"})
    folder = tmp_path / "variant"
    checkpoint = folder / "output/orthohmm_working_res/high_sensitivity_checkpoint"
    checkpoint.mkdir(parents=True)
    phylo = folder / "output/orthohmm_phylogeny"
    phylo.mkdir()
    for path in (folder / "metrics.json", checkpoint / "manifest.json",
                 phylo / "provenance_manifest.json", phylo / "reconciliation_summary.json"):
        write_json(path, {})
    pin = record(driver)
    native = {"genes": 2, "species": 2, "groups": 1, "pairs": 1, "products": {},
              "metrics": {"metadata": {"search_result_order": "canonical", "search_min_total_hits": 1},
                          "harness": {"git_commit": "revision", "exit_code": 0,
                                      "source_manifest": [{**pin, "path": "benchmark_tools/benchmark_production.py"}],
                                      "input_manifest": []}}}
    case = {"dataset": "OrthoBench", "label": "shared_hit1", "output": str(folder),
            "search_min_total_hits": 1, "search_max_candidates_per_query": 100, "inputs": []}
    plan = {"cases": [case], "source_revision": "revision", "source_root": str(source),
            "control_root": str(control_root), "inputs": [record(control_plan)]}
    report = {"source_revision": "revision", "control_report": record(control_path)}
    row = {"dataset": "OrthoBench", "execution": {"returncode": 0},
           "control": {"initial_min_total_hits": 4, "scores": {"OrthoBench F1": .8},
                       "source_revision": "control_revision"},
           "scores": {"OrthoBench F1": .82}}
    monkeypatch.setattr(reporter, "verify_input_inventory", lambda *args: None)
    monkeypatch.setattr(reporter, "validate_native", lambda *args: deepcopy(native))
    monkeypatch.setattr(reporter, "canonical_cell", lambda *args: {
        "scores": {"OrthoBench: OrthoBench F1": .8}, "coverage": {"genes": 2, "species": 2}})
    monkeypatch.setattr(reporter, "readback_scores", lambda *args: (row["scores"], {}, []))
    return plan, report, row, native


def test_shared_hit_admission_rechecks_actual_mechanism_and_matched_control(tmp_path, monkeypatch):
    plan, report, row, _ = shared_fixture(tmp_path, monkeypatch)
    cell = reporter.shared_cell(tmp_path, tmp_path, plan, report, row)
    assert cell["scores"] == {"OrthoBench: OrthoBench F1": .82}
    assert cell["matched_control"]["cell"]["scores"] == {"OrthoBench: OrthoBench F1": .8}
    assert cell["execution_kind"] == "fresh_canonical_initial_shared_hit1_cap100"


@pytest.mark.parametrize("field,value", [("search_result_order", "arrival"), ("search_min_total_hits", 4)])
def test_wrong_actual_shared_hit_mechanism_is_rejected(tmp_path, monkeypatch, field, value):
    plan, report, row, native = shared_fixture(tmp_path, monkeypatch)
    native["metrics"]["metadata"][field] = value
    with pytest.raises(ValueError, match="mechanism differs"):
        reporter.shared_cell(tmp_path, tmp_path, plan, report, row)


def test_changed_canonical_control_score_is_rejected(tmp_path, monkeypatch):
    plan, report, row, _ = shared_fixture(tmp_path, monkeypatch)
    row["control"]["scores"]["OrthoBench F1"] = .9
    with pytest.raises(ValueError, match="canonical control scores"):
        reporter.shared_cell(tmp_path, tmp_path, plan, report, row)


def test_changed_native_universe_is_rejected(tmp_path, monkeypatch):
    plan, report, row, native = shared_fixture(tmp_path, monkeypatch)
    native["genes"] = 3
    with pytest.raises(ValueError, match="universe differ"):
        reporter.shared_cell(tmp_path, tmp_path, plan, report, row)


def test_wrong_native_source_manifest_is_rejected(tmp_path, monkeypatch):
    plan, report, row, native = shared_fixture(tmp_path, monkeypatch)
    native["metrics"]["harness"]["source_manifest"][0]["sha256"] = "wrong"
    with pytest.raises(ValueError, match="source/input manifests"):
        reporter.shared_cell(tmp_path, tmp_path, plan, report, row)


def test_unpinned_canonical_control_is_rejected(tmp_path, monkeypatch):
    plan, report, row, _ = shared_fixture(tmp_path, monkeypatch)
    plan["inputs"] = []
    with pytest.raises(ValueError, match="bound to the frozen trial plan"):
        reporter.shared_cell(tmp_path, tmp_path, plan, report, row)


def test_plain_frozen_copy_git_label_cannot_override_matching_source_hashes(tmp_path, monkeypatch):
    plan, report, row, native = shared_fixture(tmp_path, monkeypatch)
    native["metrics"]["harness"]["git_commit"] = "newer_ancestor_checkout"
    cell = reporter.shared_cell(tmp_path, tmp_path, plan, report, row)
    assert cell["source_revision"] == "revision"
    assert cell["harness_git_metadata"]["git_commit"] == "newer_ancestor_checkout"


def test_plain_copy_different_git_label_does_not_excuse_changed_source_hashes(tmp_path, monkeypatch):
    plan, report, row, native = shared_fixture(tmp_path, monkeypatch)
    native["metrics"]["harness"]["git_commit"] = "newer_ancestor_checkout"
    native["metrics"]["harness"]["source_manifest"][0]["sha256"] = "changed_source"
    with pytest.raises(ValueError, match="source/input manifests"):
        reporter.shared_cell(tmp_path, tmp_path, plan, report, row)


def test_real_frozen_worktree_revision_mismatch_is_still_rejected(tmp_path, monkeypatch):
    plan, report, row, native = shared_fixture(tmp_path, monkeypatch)
    (tmp_path / "source/.git").write_text("gitdir: /fixture\n")
    native["metrics"]["harness"]["git_commit"] = "wrong_worktree_revision"
    with pytest.raises(ValueError, match="frozen worktree revision"):
        reporter.shared_cell(tmp_path, tmp_path, plan, report, row)
