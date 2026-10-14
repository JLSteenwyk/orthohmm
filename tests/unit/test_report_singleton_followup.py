from copy import deepcopy
import json

import pytest

from benchmark_tools import report_additional_hillclimb as reporter
from benchmark_tools import run_singleton_tie_followup as driver
from benchmark_tools.metric_frontier import ENDPOINTS
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.run_membership_hillclimb import write_json


def fixture(tmp_path, monkeypatch, dataset):
    root = tmp_path / dataset
    root.mkdir()
    artifact = tmp_path / "source_fixture"
    artifact.write_text("pinned")
    pin = record(artifact)
    raw = {key.split(": ", 1)[1]: .8 for key in ENDPOINTS if key.startswith(dataset + ":")}
    case = {"dataset": dataset, "case": {}, "checkpoint": "checkpoint", "checkpoint_sha": "hash",
            "candidate": {"candidate_partition": pin, "membership_constraints": pin}, "golden": {},
            "control": {"scores": {dataset + ": " + key: value for key, value in raw.items()}}}
    prepared = {"status": "prespecified_followup_inputs_prepared_not_submitted", "scientific_revision": driver.REVISION,
        "cpu": 16, "pilot_gate": {"eligible": True, "control_f1": .74, "variant_f1": .736, "delta": .736 - .74},
        "cases": [case if name == dataset else {"dataset": name} for name in driver.UNIVERSES],
        "scientific_source": str(tmp_path / "science"), "runtime": {"files": []}, "tool_path": "fixed",
        "preparer_source": pin, "inputs": [pin]}
    prepared_path = root / "prepared.json"
    write_json(prepared_path, prepared)
    plan = {"repo": str(tmp_path), "output": str(root), "dataset": dataset, "case": case,
        "scientific_source": prepared["scientific_source"], "scientific_revision": driver.REVISION,
        "driver_revision": "driver", "driver_root": str(tmp_path / "driver"),
        "prepared_inputs": record(prepared_path), "inputs": [pin], "driver_sources": [],
        "frozen_reporting_evidence": [], "runtime": prepared["runtime"], "tool_path": "fixed", "cpu": 16,
        "scoring_id": driver.scoring_id(root)}
    rows = []
    for policy in driver.POLICIES:
        folder = root / policy
        folder.mkdir()
        graph, native = ["graph", "--singleton-tie-policy", policy], ["native"]
        write_json(folder / "preflight.json", {"graph": graph, "native": native})
        write_json(folder / "graph_metrics.json", {"parameters": {"singleton_tie_policy": policy,
            "cpm_resolution": .1, "leiden_seed": 4, "profile_iterations": 1, "profile_min_species": 1, "matrix": "BLOSUM62"}})
        rows.append({"dataset": dataset, "label": policy, "status": "native_and_score_independently_verified",
                     "scores": deepcopy(raw), "candidate": case["candidate"],
                     "graph_execution": {"command": graph, "returncode": 0},
                     "native_execution": {"command": native, "returncode": 0}})
    if dataset == "QfO":
        recovery = tmp_path / "benchmarks/results/canonical_qfo_scoring_recovery_r1_v1_20261014"
        recovery.mkdir(parents=True)
        write_json(recovery / "results.json", {})
        write_json(recovery / "execution_plan.json", {"scoring_id": "recovered"})
        rows[0]["control_scoring_reuse"] = record(recovery / "results.json")
        plan["inputs"].extend([rows[0]["control_scoring_reuse"], record(recovery / "execution_plan.json")])
    plan_path = root / "execution_plan.json"
    write_json(plan_path, plan)
    report = {"status": "both_singleton_policies_scored", "source_revision": driver.REVISION,
        "driver_revision": "driver", "execution_plan": record(plan_path), "job_id": "live", "dataset": dataset, "arms": rows}
    write_json(root / "results.json", report)
    import benchmark_tools.run_canonical_downstream as producer
    import benchmark_tools.run_membership_hillclimb as membership
    monkeypatch.setattr(producer, "commands", lambda *a: (["graph"], ["native"]))
    monkeypatch.setattr(membership, "compare_golden", lambda *a: {})
    native = {"genes": driver.UNIVERSES[dataset][0], "species": driver.UNIVERSES[dataset][1],
              "groups": 1, "pairs": 1, "phylogeny": str(root), "checked": [], "products": {}}
    monkeypatch.setattr(reporter, "replay_native", lambda *a: deepcopy(native))
    calls = []
    def readback(repo, name, row, native, folder, participant):
        calls.append((name, folder, participant))
        return row["scores"], {}, []
    monkeypatch.setattr(reporter, "readback_scores", readback)
    monkeypatch.setattr(reporter, "runtime", lambda: {"files": []})
    return root, plan, report, native, calls


@pytest.mark.parametrize("dataset", list(driver.UNIVERSES))
def test_followup_plan_independently_verifies_original_prepared_inputs(tmp_path, monkeypatch, dataset):
    root, plan, report, _, _ = fixture(tmp_path, monkeypatch, dataset)
    actual, pins = reporter.singleton_followup_plan(tmp_path, root, report)
    assert actual == plan and record(root / "prepared.json") in pins
    plan["inputs"] = []
    write_json(root / "execution_plan.json", plan)
    report["execution_plan"] = record(root / "execution_plan.json")
    with pytest.raises(ValueError, match="inventory is incomplete"):
        reporter.singleton_followup_plan(tmp_path, root, report)


@pytest.mark.parametrize("dataset", list(driver.UNIVERSES))
def test_verified_variant_rechecks_its_matched_control(tmp_path, monkeypatch, dataset):
    root, plan, report, _, calls = fixture(tmp_path, monkeypatch, dataset)
    cell = reporter.singleton_followup_cell(tmp_path, root, plan, report, report["arms"][1])
    assert cell["scores"] == plan["case"]["control"]["scores"]
    assert cell["matched_control"]["singleton_tie_policy"] == "first_hit"
    assert cell["source_revision"] == driver.REVISION
    if dataset == "QfO":
        assert calls == [(dataset, root / "defer_ambiguous", "ohmm_hillclimb_r1_" + plan["scoring_id"]),
                         (dataset, tmp_path / "benchmarks/results/canonical_qfo_scoring_recovery_r1_v1_20261014",
                          "ohmm_hillclimb_r1_recovered")]


@pytest.mark.parametrize("key,value", [("singleton_tie_policy", "first_hit"), ("cpm_resolution", .2),
    ("leiden_seed", 1), ("profile_iterations", 2), ("profile_min_species", 2), ("matrix", "BLOSUM45")])
def test_changed_actual_graph_settings_rejected(tmp_path, monkeypatch, key, value):
    root, plan, report, _, _ = fixture(tmp_path, monkeypatch, "ThreeKingdoms")
    path = root / "defer_ambiguous/graph_metrics.json"
    graph = json.loads(path.read_text())
    graph["parameters"][key] = value
    write_json(path, graph)
    with pytest.raises(ValueError, match="graph settings differ"):
        reporter.singleton_followup_cell(tmp_path, root, plan, report, report["arms"][1])


@pytest.mark.parametrize("status", ["graph_running", "native_running", "failed", "not_run_failed_control_gate"])
def test_unverified_control_never_admits_variant(tmp_path, monkeypatch, status):
    root, plan, report, _, _ = fixture(tmp_path, monkeypatch, "ThreeKingdoms")
    report["arms"][0]["status"] = status
    with pytest.raises(ValueError, match="verified singleton follow-up"):
        reporter.singleton_followup_cell(tmp_path, root, plan, report, report["arms"][1])


def test_changed_native_gene_universe_rejected(tmp_path, monkeypatch):
    root, plan, report, native, _ = fixture(tmp_path, monkeypatch, "QfO")
    native["genes"] -= 1
    with pytest.raises(ValueError, match="gene universe differs"):
        reporter.singleton_followup_cell(tmp_path, root, plan, report, report["arms"][1])


def test_qfo_reused_assessment_must_remain_pinned(tmp_path, monkeypatch):
    root, plan, report, _, _ = fixture(tmp_path, monkeypatch, "QfO")
    plan["inputs"].pop()
    with pytest.raises(ValueError, match="bound to the follow-up"):
        reporter.singleton_followup_cell(tmp_path, root, plan, report, report["arms"][1])


def test_changed_actual_command_rejected(tmp_path, monkeypatch):
    root, plan, report, _, _ = fixture(tmp_path, monkeypatch, "ThreeKingdoms")
    report["arms"][1]["native_execution"]["command"] = ["other"]
    with pytest.raises(ValueError, match="actual execution differs"):
        reporter.singleton_followup_cell(tmp_path, root, plan, report, report["arms"][1])


@pytest.mark.parametrize("invalid", [None, "path", "hash"])
def test_frozen_reporting_evidence_preserves_admitted_code_identity(tmp_path, monkeypatch, invalid):
    root, plan, report, _, _ = fixture(tmp_path, monkeypatch, "ThreeKingdoms")
    original = tmp_path / "benchmark_tools/reporter.py"
    original.parent.mkdir()
    original.write_text("admitted source")
    pin = record(original)
    captured = tmp_path / "driver/evidence/benchmark_tools/reporter.py"
    captured.parent.mkdir(parents=True)
    captured.write_text("admitted source")
    frozen = record(captured)
    prepared_path = root / "prepared.json"
    prepared = json.loads(prepared_path.read_text())
    prepared["inputs"].append(pin)
    write_json(prepared_path, prepared)
    plan["prepared_inputs"] = record(prepared_path)
    plan["frozen_reporting_evidence"] = [{"original": pin, "frozen": deepcopy(frozen)}]
    plan["inputs"].append(frozen)
    original.write_text("later reporting edit, not a scientific change")
    if invalid == "path":
        plan["frozen_reporting_evidence"][0]["frozen"]["path"] = str(tmp_path / "outside.py")
    if invalid == "hash":
        plan["frozen_reporting_evidence"][0]["frozen"]["sha256"] = "different"
    write_json(root / "execution_plan.json", plan)
    report["execution_plan"] = record(root / "execution_plan.json")
    if invalid:
        with pytest.raises(ValueError, match="evidence substitution"):
            reporter.singleton_followup_plan(tmp_path, root, report)
    else:
        actual, _ = reporter.singleton_followup_plan(tmp_path, root, report)
        assert actual["frozen_reporting_evidence"][0]["original"] == pin


def test_failed_producer_observation_never_erases_previously_admitted_scores(tmp_path, monkeypatch):
    root, plan, report, _, _ = fixture(tmp_path, monkeypatch, "ThreeKingdoms")
    good = reporter.singleton_followup_cell(tmp_path, root, plan, report, report["arms"][1])
    report["arms"][1].update(status="failed", error="later failed observation")
    write_json(root / "results.json", report)
    table = tmp_path / "old.tsv"
    table.write_text("old")
    parent = tmp_path / "old.json"
    write_json(parent, {"score_table": record(table), "methods": [{"key": "singleton_tie:defer_ambiguous",
        "label": "singleton", "scores": good["scores"], "cells": [good]}], "round_screening": {}, "limitations": []})
    result = reporter.build(tmp_path, parent, tmp_path / "out", singleton_tk_root=root)
    assert result["methods"][0]["scores"] == good["scores"]
    assert any(row.get("error") == "later failed observation" for row in result["methods"][0]["cells"])


@pytest.mark.parametrize("dataset", list(driver.UNIVERSES))
@pytest.mark.parametrize("status", ["native_and_score_independently_verified", "graph_running", "failed", None])
def test_dataset_admission_preserves_ob_and_unrelated_methods(tmp_path, monkeypatch, dataset, status):
    root, plan, report, _, _ = fixture(tmp_path, monkeypatch, dataset)
    if status is None:
        report["arms"].pop()
    else:
        report["arms"][1]["status"] = status
    write_json(root / "results.json", report)
    table = tmp_path / "old.tsv"
    table.write_text("old table")
    unrelated = {"key": "other", "label": "other", "scores": {}, "cells": []}
    ob = {"dataset": "OrthoBench", "status": "independent_content_readback_verified",
          "scores": {"OrthoBench: OrthoBench F1": .7365896292967041}}
    method = {"key": "singleton_tie:defer_ambiguous", "label": "singleton", "scores": ob["scores"], "cells": [ob]}
    parent = tmp_path / "old.json"
    write_json(parent, {"score_table": record(table), "methods": [unrelated, method],
                        "round_screening": {}, "limitations": []})
    kwargs = {"singleton_qfo_root" if dataset == "QfO" else "singleton_tk_root": root}
    result = reporter.build(tmp_path, parent, tmp_path / "out", **kwargs)
    assert result["methods"][0] == unrelated
    updated = result["methods"][1]
    assert updated["scores"]["OrthoBench: OrthoBench F1"] == .7365896292967041
    expected = 1 + (6 if dataset == "QfO" else 1) if status == "native_and_score_independently_verified" else 1
    assert len(updated["scores"]) == expected
    assert result["defaults_promoted"] is False and result["goal_completed"] is False
    assert not updated["complete_eight_endpoints"]


def test_both_verified_followups_complete_the_same_eight_endpoint_vector(tmp_path, monkeypatch):
    qfo_root, _, _, _, _ = fixture(tmp_path, monkeypatch, "QfO")
    tk_root, _, _, _, _ = fixture(tmp_path, monkeypatch, "ThreeKingdoms")
    def native(row, *args):
        genes, species = driver.UNIVERSES[row["dataset"]]
        return {"genes": genes, "species": species, "groups": 1, "pairs": 1,
                "phylogeny": str(tmp_path), "checked": [], "products": {}}
    monkeypatch.setattr(reporter, "replay_native", native)
    table = tmp_path / "old.tsv"
    table.write_text("old")
    ob = {"dataset": "OrthoBench", "status": "independent_content_readback_verified",
          "scores": {"OrthoBench: OrthoBench F1": .7365896292967041}}
    parent = tmp_path / "old.json"
    write_json(parent, {"score_table": record(table), "methods": [{"key": "singleton_tie:defer_ambiguous",
        "label": "singleton", "scores": ob["scores"], "cells": [ob]}], "round_screening": {}, "limitations": []})
    result = reporter.build(tmp_path, parent, tmp_path / "out", singleton_qfo_root=qfo_root, singleton_tk_root=tk_root)
    method = result["methods"][0]
    assert set(method["scores"]) == set(ENDPOINTS) and method["complete_eight_endpoints"] is True
    assert {cell["dataset"] for cell in method["cells"]} == {"OrthoBench", "QfO", "ThreeKingdoms"}
    assert result["defaults_promoted"] is False
