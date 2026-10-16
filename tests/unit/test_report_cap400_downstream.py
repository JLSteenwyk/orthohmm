from copy import deepcopy
import json
from pathlib import Path
import sys

import pytest

from benchmark_tools import report_additional_hillclimb as reporter
from benchmark_tools import report_cap400_downstream as recovery
from benchmark_tools.continue_cap400_downstream import REVISION, phylogeny_command
from benchmark_tools.metric_frontier import ENDPOINTS
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.recover_full_qfo_search import short_id
from benchmark_tools.run_membership_hillclimb import PRODUCTS, write_json


ACCOUNTING = "JobIDRaw|State|ExitCode|AllocCPUS|NodeList\n24350|COMPLETED|0:0|32|bizon\n"


def test_terminal_scheduler_requires_exact_success_and_resources(monkeypatch):
    monkeypatch.setattr(recovery.subprocess, "check_output", lambda *a, **kw: ACCOUNTING)
    assert recovery.terminal({"job_id": "24350"})["State"] == "COMPLETED"
    with pytest.raises(ValueError):
        recovery.terminal({"job_id": "other"})


@pytest.mark.parametrize("text", [ACCOUNTING.replace("COMPLETED", "RUNNING"),
    ACCOUNTING.replace("0:0", "0:11"), ACCOUNTING.replace("|32|", "|1|"),
    ACCOUNTING.replace("bizon", "other"), ACCOUNTING + "24350|COMPLETED|0:0|32|bizon\n"])
def test_nonterminal_or_changed_scheduler_rejected(monkeypatch, text):
    monkeypatch.setattr(recovery.subprocess, "check_output", lambda *a, **kw: text)
    with pytest.raises(ValueError):
        recovery.terminal({"job_id": "24350"})


def panel(tmp_path, monkeypatch, status):
    root = tmp_path / "producer"
    root.mkdir()
    write_json(root / "results.json", {"status": status, "job_id": "24350", "scores": {"VGNC F": 1}})
    old_table = tmp_path / "old.tsv"
    old_table.write_text("old table\n")
    other = {"key": "other", "label": "other", "scores": {ENDPOINTS[0]: .5}, "cells": []}
    old_skip = {"dataset": "QfO", "status": "not_run_failed_control", "scores": {}}
    method = {"key": "full_search:cap400", "label": "cap400", "scores": {ENDPOINTS[0]: .72},
              "cells": [old_skip]}
    parent = tmp_path / "parent.json"
    write_json(parent, {"methods": [other, method], "score_table": record(old_table),
                        "round_screening": {}, "limitations": []})
    failed = {"dataset": "QfO", "status": "failed", "scores": {}, "producer_job_id": "24245",
              "error": "retained original final Leiden crash"}
    monkeypatch.setattr(recovery, "load_plan", lambda *args: ({"runtime": {"files": []}}, [], deepcopy(failed)))
    monkeypatch.setattr(reporter, "runtime", lambda: {"files": []})
    monkeypatch.setattr(recovery.subprocess, "check_output", lambda *a, **kw: ACCOUNTING)
    return root, parent, other, old_skip, failed


@pytest.mark.parametrize("status", ["original_inferred_phylogeny_running",
    "failed_cap400_downstream_recovery", "native_verified_scoring_pending"])
def test_pending_or_failed_scores_never_admitted_and_original_failure_kept(tmp_path, monkeypatch, status):
    root, parent, other, old_skip, failed = panel(tmp_path, monkeypatch, status)
    monkeypatch.setattr(recovery, "cell", lambda *args: pytest.fail("Unverified producer scores admitted"))
    result = reporter.build(tmp_path, parent, tmp_path / "out", cap400_recovery_root=root)
    assert result["methods"][0] == other
    method = result["methods"][1]
    assert method["scores"] == {ENDPOINTS[0]: .72}
    assert old_skip in method["cells"] and failed in method["cells"]
    assert "NA" in (tmp_path / "out/scores.tsv").read_text()
    assert result["defaults_promoted"] is result["goal_completed"] is False


def test_completed_recovery_updates_only_cap400_qfo_and_keeps_failure_history(tmp_path, monkeypatch):
    root, parent, other, old_skip, failed = panel(tmp_path, monkeypatch, recovery.SUCCESS)
    scores = {k: .81 for k in ENDPOINTS if k.startswith("QfO:")}
    good = {"dataset": "QfO", "status": "independent_content_readback_verified", "scores": scores,
            "products": {name: {"bytes": 1, "sha256": "same", "path": name} for name in PRODUCTS}}
    monkeypatch.setattr(recovery, "cell", lambda *args: deepcopy(good))
    result = reporter.build(tmp_path, parent, tmp_path / "out", cap400_recovery_root=root)
    assert result["methods"][0] == other
    method = result["methods"][1]
    assert method["scores"] == {ENDPOINTS[0]: .72, **scores}
    assert old_skip in method["cells"] and failed in method["cells"]
    assert (tmp_path / "out/snapshots/cap400_recovery_reference.json").exists()
    assert not (tmp_path / "out/snapshots/cap400_recovery.json").exists()


@pytest.fixture
def completed(tmp_path, monkeypatch):
    root, source, science = (tmp_path / n for n in ("producer", "source", "science"))
    root.mkdir()
    (science / "benchmark_tools").mkdir(parents=True)
    replay_source = science / "benchmark_tools/replay_phylogeny.py"
    replay_source.write_text("# original fixture\n")
    artifact = tmp_path / "artifact"
    artifact.write_text("native fixture\n")
    pin = record(artifact)
    control = {"report": pin, "scores": {"VGNC F": .8}}
    admission = {"complete_gene_order": {"path": "/checkpoint/gene_names.txt"}, "partition": pin,
                 "producer": pin, "execution_plan": pin}
    graph_plan = {"output": "/graph"}
    plan = {"repo": "/r", "output": str(root), "source_root": str(source),
        "scientific_source": str(science), "scientific_revision": REVISION, "fasta": "/fasta",
        "aligner": "mafft", "tree_builder": "FastTree", "runtime": {}, "tool_path": "/bin",
        "control": control, "checkpoint": "/checkpoint", "partition": pin, "payload": "/graph/payload",
        "scoring_id": short_id(Path("/r"), root, 1), "checked": [], "sources": [], "admission": pin}
    plan["command"] = phylogeny_command(plan)
    stage = {"refined_groups": 2, "candidate_groups": 1, "products": [], "modules": []}
    stage_path = root / "graph_report.json"
    write_json(stage_path, stage)
    metrics = {"status": "complete", "source": record(replay_source),
               "git": {"commit": REVISION, "dirty": False},
               "parameters": {"species_tree_mode": "infer", "species_tree": None,
                   "aligner": "mafft", "tree_builder": "FastTree", "root_duplication_rule": "species_overlap",
                   "pair_orthology_rule": "positive_paralogy", "species_tree_rooting": "min_variance", "cpu": 32,
                   "satellite_membership_policy": "high_confidence_pair", "explicit_unconstrained_ablation": False,
                   "checkpoint_source": None}}
    working = root / "graph/orthohmm_working_res"
    working.mkdir(parents=True)
    candidate, constraints = working / "phylogeny_candidate_superfamilies.txt", working / "phylogeny_candidate_merges.json"
    candidate.write_text("native fixture\n")
    constraints.write_text("[]\n")
    metrics["input"] = {"candidate_clusters": record(candidate), "membership_constraints": record(constraints),
                        "fasta_directory": "/fasta"}
    write_json(root / "native_metrics.json", metrics)
    raw_scores = {key.split(": ", 1)[1]: .81 for key in ENDPOINTS if key.startswith("QfO:")}
    control["scores"] = {key: .8 for key in raw_scores}
    write_json(root / "execution_plan.json", plan)
    report = {"status": recovery.SUCCESS, "job_id": "24350", "postflight_pins_verified": True,
        "initial_hmm_search_rerun": False, "profile_expansion_rerun": False, "graph_clustering_rerun": False,
        "graph_execution": {"returncode": 0, "command": [sys.executable,
            str(source / "benchmark_tools/continue_cap400_downstream_worker.py"), "--execution-plan",
            str(root / "execution_plan.json")]},
        "execution": {"returncode": 0, "command": plan["command"]},
        "execution_plan": record(root / "execution_plan.json"),
        "graph": {"report": record(stage_path), "refined_groups": 2, "candidate_groups": 1},
        "scores": raw_scores, "deltas": {key: .81 - .8 for key in raw_scores},
        "control": control, "metric_files": []}
    original_path = root / "original.json"
    write_json(original_path, {"control": control, "execution": {"failed_fresh": True}})
    graph_result = root / "graph_recovery.json"
    write_json(graph_result, {"execution": {"cached_graph": True}})
    admission["producer"] = record(graph_result)
    parent = {"scientific_source": str(science), "runtime": {}, "tool_path": "/bin",
              "case": {"fresh_case": {"fasta": "/fasta", "command": ["--aligner", "mafft", "--tree-builder", "FastTree"]}}}
    failure = {"parent_result": record(original_path)}
    write_json(root / "results.json", report)
    native = {"genes": 984137, "species": 78, "groups": 10, "pairs": 20,
              "products": {name: pin for name in PRODUCTS}, "checked": [pin]}
    monkeypatch.setattr(recovery, "verify_graph", lambda *a: (parent, failure, admission, graph_plan))
    monkeypatch.setattr(recovery, "validate_graph_stage", lambda *a: None)
    monkeypatch.setattr(recovery, "replay_native", lambda *a: deepcopy(native))
    monkeypatch.setattr(recovery, "readback_scores", lambda *a: (deepcopy(raw_scores), {}, []))
    monkeypatch.setattr(recovery, "control_readback", lambda *a: ({"report": pin}, []))
    monkeypatch.setattr(recovery.subprocess, "check_output", lambda *a, **kw: ACCOUNTING)
    control["checked"] = []
    write_json(original_path, {"control": control, "execution": {"failed_fresh": True}})
    failure["parent_result"] = record(original_path)
    return root, plan, report, native, metrics


def test_completed_cell_has_separate_cost_scopes_and_all_scores(tmp_path, completed):
    root, plan, report, _, _ = completed
    cell = recovery.cell(tmp_path, root, plan, report)
    assert len(cell["scores"]) == 6
    assert cell["resources"]["failed_fresh_execution"] == {"failed_fresh": True}
    assert cell["resources"]["cached_final_graph_recovery"] == {"cached_graph": True}
    assert cell["execution_kind"] == "cached_final_graph_recovery_of_original_arrival_cap400"


@pytest.mark.parametrize("field,value", [("status", "running"), ("postflight_pins_verified", False),
    ("initial_hmm_search_rerun", True), ("profile_expansion_rerun", True), ("graph_clustering_rerun", True)])
def test_incomplete_or_changed_execution_refused(tmp_path, completed, field, value):
    root, plan, report, _, _ = completed
    report[field] = value
    with pytest.raises(ValueError):
        recovery.cell(tmp_path, root, plan, report)


@pytest.mark.parametrize("field,value", [("genes", 984136), ("species", 77)])
def test_incomplete_native_universe_refused(tmp_path, completed, field, value):
    root, plan, report, native, _ = completed
    native[field] = value
    with pytest.raises(ValueError, match="gene universe"):
        recovery.cell(tmp_path, root, plan, report)


def test_wrong_producer_delta_refused(tmp_path, completed):
    root, plan, report, _, _ = completed
    report["deltas"]["VGNC F"] = 1
    with pytest.raises(ValueError, match="delta"):
        recovery.cell(tmp_path, root, plan, report)


@pytest.mark.parametrize("change", [{"cpu": 16}, {"checkpoint_source": "/unplanned"}])
def test_wrong_native_scope_refused(tmp_path, completed, change):
    root, plan, report, _, metrics = completed
    metrics["parameters"].update(change)
    write_json(root / "native_metrics.json", metrics)
    with pytest.raises(ValueError, match="metrics"):
        recovery.cell(tmp_path, root, plan, report)


@pytest.mark.parametrize("stage", ["graph_execution", "execution"])
def test_different_actual_command_is_rejected(tmp_path, completed, stage):
    root, plan, report, _, _ = completed
    report[stage]["command"] = ["changed"]
    with pytest.raises(ValueError, match="Actual"):
        recovery.cell(tmp_path, root, plan, report)


def test_changed_native_candidate_input_is_rejected(tmp_path, completed):
    root, plan, report, _, metrics = completed
    metrics["input"]["candidate_clusters"]["sha256"] = "changed"
    write_json(root / "native_metrics.json", metrics)
    with pytest.raises(ValueError, match="inputs"):
        recovery.cell(tmp_path, root, plan, report)


@pytest.fixture
def control_fixture(tmp_path, monkeypatch):
    root = tmp_path / "control"
    root.mkdir()
    raw = {key.split(": ", 1)[1]: .8 for key in ENDPOINTS if key.startswith("QfO:")}
    write_json(root / "execution_plan.json", {"command": ["original-control"]})
    write_json(root / "scoring_preflight.json", {"stage": {"participant": "matched-original-control"}})
    write_json(root / "results.json", {"status": "cached_native_control_and_all_six_scores_verified",
        "execution_plan": record(root / "execution_plan.json"), "scores": raw})
    claimed = {"report": record(root / "results.json"), "scores": deepcopy(raw)}
    plan, report, original = ({"control": deepcopy(claimed), "scientific_revision": REVISION},
                              {"control": deepcopy(claimed)}, {"control": deepcopy(claimed)})
    native = {"genes": 984137, "species": 78, "groups": 1, "pairs": 1,
              "products": {}, "checked": []}
    monkeypatch.setattr(recovery, "replay_native", lambda *a: native)
    def raw_reader(repo, dataset, row, candidate, folder, participant):
        assert participant == "matched-original-control" and folder == root
        return deepcopy(raw), {}, []
    monkeypatch.setattr(recovery, "readback_scores", raw_reader)
    return plan, report, original, raw, native


def test_matched_control_is_independently_read_back(tmp_path, control_fixture):
    plan, report, original, _, _ = control_fixture
    cell, pins = recovery.control_readback(tmp_path, plan, report, original)
    assert len(cell["scores"]) == 6 and all(value == .8 for value in cell["scores"].values())
    assert pins and cell["coverage"]["genes"] == 984137


@pytest.mark.parametrize("which", [0, 1])
def test_changed_claimed_control_is_rejected(tmp_path, control_fixture, which):
    args = list(control_fixture[:3])
    args[which]["control"]["scores"]["VGNC F"] = 1
    with pytest.raises(ValueError, match="matched original control"):
        recovery.control_readback(tmp_path, *args)


def test_changed_raw_control_score_is_rejected(tmp_path, control_fixture):
    plan, report, original, raw, _ = control_fixture
    raw["VGNC F"] = 1
    with pytest.raises(ValueError, match="raw scores"):
        recovery.control_readback(tmp_path, plan, report, original)


def test_incomplete_control_universe_is_rejected(tmp_path, control_fixture):
    plan, report, original, _, native = control_fixture
    native["genes"] = 1
    with pytest.raises(ValueError, match="content"):
        recovery.control_readback(tmp_path, plan, report, original)
