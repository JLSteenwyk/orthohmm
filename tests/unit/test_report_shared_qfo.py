from copy import deepcopy
import hashlib
import json

import pytest

from benchmark_tools import continue_shared_qfo as continuation
from benchmark_tools import report_additional_hillclimb as reporter
from benchmark_tools.metric_frontier import ENDPOINTS
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_full_search_hillclimb import native_command
from benchmark_tools.run_membership_hillclimb import DATASETS, PRODUCTS, write_json


@pytest.fixture
def trial(tmp_path, monkeypatch):
    root, parent_root, recovery_root = (tmp_path / name for name in ("trial", "original", "recovery"))
    for folder in (root, parent_root, recovery_root):
        folder.mkdir()
    science = tmp_path / "science"
    old_folder = tmp_path / "original_native"
    fasta = tmp_path / "fasta"
    command = native_command(science, fasta, old_folder, ["mafft", "FastTree"], 100, 16)
    command.extend(["--search-result-order", "canonical", "--search-min-total-hits", "1"])
    qfo = {"dataset": "QfO", "label": "shared_hit1", "search_min_total_hits": 1,
           "search_max_candidates_per_query": 100, "command": command, "output": str(old_folder),
           "fasta": str(fasta), "inputs": [{"path": f"fasta{i}"} for i in range(78)]}
    parent = {"source_revision": continuation.REVISION, "source_root": str(science), "cpu": 16,
              "cases": [qfo if dataset == "QfO" else {"dataset": dataset} for dataset in DATASETS]}
    parent_path = parent_root / "execution_plan.json"
    write_json(parent_path, parent)
    archived = {"status": "finished_with_explicit_gaps", "job_id": "24242",
                "source_revision": continuation.REVISION, "execution_plan": record(parent_path),
                "arms": [{"dataset": dataset, "label": "shared_hit1", "scores": {},
                          "status": "not_run_unverified_control" if dataset == "QfO" else "other",
                          "error_type": "ValueError",
                          "error": "Canonical control unavailable or not independently verified"} for dataset in DATASETS]}
    write_json(root / "original_skip_snapshot.json", archived)
    recovery_plan = recovery_root / "execution_plan.json"
    write_json(recovery_plan, {"repo": str(tmp_path), "output": str(recovery_root), "inputs": [],
                               "sources": [], "runtime": {"files": []}})
    write_json(recovery_root / "results.json", {"execution_plan": record(recovery_plan), "job_id": "24246"})
    case = continuation.relocated_case(parent, root)
    folder = root / "native"
    phylo = folder / "output/orthohmm_phylogeny"
    phylo.mkdir(parents=True)
    (folder / "output/orthohmm_working_res/high_sensitivity_checkpoint").mkdir(parents=True)
    for path in (folder / "metrics.json", phylo / "provenance_manifest.json", phylo / "reconciliation_summary.json"):
        write_json(path, {})
    artifact = root / "product"
    artifact.write_text("native fixture\n")
    products = {name: record(artifact) for name in PRODUCTS}
    coverage = {"genes": 984137, "species": 78, "groups": 42, "pairs": 123}
    control = {"status": "independent_content_readback_verified", "dataset": "QfO",
               "scores": {key: .8 for key in ENDPOINTS if key.startswith("QfO:")},
               "products": products, "coverage": coverage, "source_revision": "control-frozen", "checked": []}
    plan = {"repo": str(tmp_path), "output": str(root), "scientific_revision": continuation.REVISION,
            "cpu": 16, "dependency_jobs": ["24242", "24246"], "parent_root": str(parent_root),
            "parent_plan": record(parent_path), "scientific_source": str(science), "case": case,
            "recovery_root": str(recovery_root), "inputs": [record(parent_path), record(recovery_plan)],
            "sources": [], "runtime": {"files": []},
            "scoring_id": "sh_" + hashlib.sha256(str(root).encode()).hexdigest()[:12]}
    write_json(root / "execution_plan.json", plan)
    report = {"status": "original_shared_qfo_native_and_all_six_scores_verified", "job_id": "24247",
              "execution_plan": record(root / "execution_plan.json"),
              "dependency_accounting": "24242|FAILED|1:0\n24246|COMPLETED|0:0\n",
              "original_skip": record(root / "original_skip_snapshot.json"),
              "execution": {"command": deepcopy(case["command"]), "returncode": 0},
              "control": {"report": record(recovery_root / "results.json"), "cell": deepcopy(control)},
              "scores": {key.split(": ", 1)[1]: .81 for key in control["scores"]}, "checked": []}
    native = {**coverage, "products": products, "metrics": {}}
    monkeypatch.setattr(reporter, "runtime", lambda: {"files": []})
    monkeypatch.setattr(reporter, "canonical_recovery_cell", lambda *args: deepcopy(control))
    monkeypatch.setattr(reporter, "shared_native", lambda *args: deepcopy(native))
    monkeypatch.setattr(reporter, "readback_scores", lambda repo, dataset, row, *args: (row["scores"], {}, []))
    return root, plan, report, control, native


def test_shared_qfo_admission_verifies_control_and_retains_skip(tmp_path, trial):
    root, plan, report, control, _ = trial
    cell = reporter.shared_qfo_cell(tmp_path, root, plan, report)
    assert len(cell["scores"]) == 6 and all(value == .81 for value in cell["scores"].values())
    assert cell["matched_control"]["cell"] == control
    assert cell["source_revision"] == continuation.REVISION
    assert cell["original_unrun_attempt"]["status"] == "not_run_unverified_control"
    assert cell["original_unrun_attempt"]["scores"] == {}
    assert cell["resources"] == report["execution"]
    assert cell["execution_kind"] == "fresh_original_canonical_shared_hit1_cap100_continuation"


@pytest.mark.parametrize("key,value", [("cpu", 32), ("scientific_revision", "new"),
                                       ("dependency_jobs", ["24242", "other"])])
def test_shared_qfo_wrong_science_or_resources_rejected(tmp_path, trial, key, value):
    root, plan, report, _, _ = trial
    plan[key] = value
    with pytest.raises(ValueError, match="original-science"):
        reporter.shared_qfo_cell(tmp_path, root, plan, report)


@pytest.mark.parametrize("change", [{"returncode": 1}, {"command": ["other"]}])
def test_shared_qfo_actual_execution_must_match(tmp_path, trial, change):
    root, plan, report, _, _ = trial
    report["execution"].update(change)
    with pytest.raises(ValueError, match="actual command"):
        reporter.shared_qfo_cell(tmp_path, root, plan, report)


def test_shared_qfo_wrong_scoring_namespace_rejected(tmp_path, trial):
    root, plan, report, _, _ = trial
    plan["scoring_id"] = "another"
    with pytest.raises(ValueError, match="scoring namespace"):
        reporter.shared_qfo_cell(tmp_path, root, plan, report)


@pytest.mark.parametrize("key,value", [("scientific_source", "another"), ("parent_root", "another")])
def test_shared_qfo_wrong_original_source_or_parent_rejected(tmp_path, trial, key, value):
    root, plan, report, _, _ = trial
    plan[key] = value
    with pytest.raises(ValueError, match="original scientific"):
        reporter.shared_qfo_cell(tmp_path, root, plan, report)


def test_shared_qfo_does_not_admit_changed_case(tmp_path, trial):
    root, plan, report, _, _ = trial
    plan["case"]["search_max_candidates_per_query"] = 400
    with pytest.raises(ValueError, match="original scientific"):
        reporter.shared_qfo_cell(tmp_path, root, plan, report)


@pytest.mark.parametrize("key", ["scores", "products", "coverage", "source_revision"])
def test_shared_qfo_claimed_control_must_match_independent_readback(tmp_path, trial, key):
    root, plan, report, _, _ = trial
    report["control"]["cell"][key] = "changed"
    with pytest.raises(ValueError, match="control content"):
        reporter.shared_qfo_cell(tmp_path, root, plan, report)


@pytest.mark.parametrize("key,value", [("genes", 984136), ("species", 77)])
def test_shared_qfo_full_universe_required(tmp_path, trial, key, value):
    root, plan, report, _, native = trial
    native[key] = value
    with pytest.raises(ValueError, match="gene universe"):
        reporter.shared_qfo_cell(tmp_path, root, plan, report)


def test_shared_qfo_original_inference_must_have_been_unrun(tmp_path, trial):
    root, plan, report, _, _ = trial
    path = root / "original_skip_snapshot.json"
    archived = json.loads(path.read_text())
    archived["arms"][1]["status"] = "failed"
    write_json(path, archived)
    report["original_skip"] = record(path)
    with pytest.raises(ValueError, match="do not duplicate"):
        reporter.shared_qfo_cell(tmp_path, root, plan, report)


def test_shared_qfo_live_dependency_cannot_support_admission(tmp_path, trial):
    root, plan, report, _, _ = trial
    report["dependency_accounting"] = "24242|RUNNING|0:0\n24246|COMPLETED|0:0\n"
    with pytest.raises(ValueError, match="terminal original"):
        reporter.shared_qfo_cell(tmp_path, root, plan, report)


def test_shared_qfo_wrong_control_job_rejected(tmp_path, trial):
    root, plan, report, _, _ = trial
    path = root.parent / "recovery/results.json"
    recovery = json.loads(path.read_text())
    recovery["job_id"] = "another"
    write_json(path, recovery)
    report["control"]["report"] = record(path)
    with pytest.raises(ValueError, match="planned dependency"):
        reporter.shared_qfo_cell(tmp_path, root, plan, report)


@pytest.mark.parametrize("key", ["parent_plan", "recovery_plan"])
def test_shared_qfo_controls_must_be_pinned_in_frozen_plan(tmp_path, trial, key):
    root, plan, report, _, _ = trial
    index = 0 if key == "parent_plan" else 1
    plan["inputs"].pop(index)
    with pytest.raises(ValueError, match="original scientific|planned dependency"):
        reporter.shared_qfo_cell(tmp_path, root, plan, report)


def test_shared_qfo_changed_archived_skip_hash_rejected(tmp_path, trial):
    root, plan, report, _, _ = trial
    (root / "original_skip_snapshot.json").write_text("changed")
    with pytest.raises(ValueError):
        reporter.shared_qfo_cell(tmp_path, root, plan, report)


def test_shared_qfo_full_verification_inventory_checked_but_not_duplicated(tmp_path, trial):
    root, plan, report, _, _ = trial
    artifact = root / "extra_inventory_artifact"
    artifact.write_text("verified")
    report["checked"] = [record(artifact)]
    cell = reporter.shared_qfo_cell(tmp_path, root, plan, report)
    assert record(artifact) not in cell["checked"]
    assert cell["verification_record_count"] > 0
    artifact.write_text("changed")
    with pytest.raises(ValueError):
        reporter.shared_qfo_cell(tmp_path, root, plan, report)


@pytest.mark.parametrize("status", ["fresh_original_shared_qfo_running", "native_verified_scoring_pending",
                                    "failed_or_unverified_shared_qfo_continuation",
                                    "original_shared_qfo_native_and_all_six_scores_verified"])
def test_shared_qfo_build_admits_only_verified_scores(tmp_path, monkeypatch, trial, status):
    root, _, report, control, _ = trial
    report["status"] = status
    write_json(root / "results.json", report)
    table = tmp_path / "parent.tsv"
    table.write_text("old table\n")
    unrelated = {"key": "other", "label": "other", "scores": {}, "cells": []}
    parent = tmp_path / "parent.json"
    write_json(parent, {"score_table": record(table), "methods": [unrelated],
                        "round_screening": {}, "limitations": []})
    if status != "original_shared_qfo_native_and_all_six_scores_verified":
        monkeypatch.setattr(reporter, "shared_qfo_cell", lambda *args: pytest.fail("Unverified scores admitted"))
    result = reporter.build(tmp_path, parent, tmp_path / "out", shared_qfo_root=root)
    assert result["methods"][0] == unrelated
    method = next(row for row in result["methods"] if row["key"] == "shared_hit:shared_hit1")
    if status == "original_shared_qfo_native_and_all_six_scores_verified":
        assert len(method["scores"]) == 6
        assert method["cells"][0]["status"] == "not_run_unverified_control"
        matched = next(row for row in result["methods"] if row["key"] == "canonical_collection:canonical")
        assert matched["scores"] == control["scores"]
    else:
        assert method["scores"] == {}
    assert not method["complete_eight_endpoints"]
    assert result["defaults_promoted"] is False and result["goal_completed"] is False
