from copy import deepcopy
import json
import sys

import pytest

from benchmark_tools import audit_cap400_scoring_parent as audit
from benchmark_tools.prepare_ob_candidate_neighborhood import record


@pytest.fixture
def failure():
    return {"job_id": "24350", "status": "failed_cap400_downstream_recovery",
        "error_type": "ModuleNotFoundError", "error": "No module named 'qfo_benchmark'",
        "scores": {}, "postflight_pins_verified": True, "initial_hmm_search_rerun": False,
        "profile_expansion_rerun": False, "graph_clustering_rerun": False}


ROW = {"JobIDRaw": "24350", "State": "FAILED", "ExitCode": "1:0",
       "AllocCPUS": "32", "NodeList": "bizon"}


def test_specific_terminal_import_failure_accepted(failure):
    assert audit.require_failure(failure, [ROW]) == ROW


@pytest.mark.parametrize("change", [{"job_id": "other"}, {"status": "running"},
    {"error": "different failure"}, {"error_type": "ValueError"}, {"scores": {"GO": .5}},
    {"postflight_pins_verified": False}, {"initial_hmm_search_rerun": True},
    {"profile_expansion_rerun": True}, {"graph_clustering_rerun": True}, {"postflight_error": "changed"}])
def test_different_failure_or_scope_refused(failure, change):
    failure.update(change)
    with pytest.raises(ValueError):
        audit.require_failure(failure, [ROW])


@pytest.mark.parametrize("rows", [[], [ROW, ROW], [{**ROW, "State": "RUNNING"}],
    [{**ROW, "ExitCode": "0:0"}], [{**ROW, "AllocCPUS": "8"}], [{**ROW, "NodeList": "other"}]])
def test_nonterminal_or_wrong_scheduler_refused(failure, rows):
    with pytest.raises(ValueError):
        audit.require_failure(failure, rows)


@pytest.fixture
def native(tmp_path):
    source = tmp_path / "science/benchmark_tools/replay_phylogeny.py"
    source.parent.mkdir(parents=True)
    source.write_text("# frozen scientific fixture\n")
    candidate, constraints = tmp_path / "candidates", tmp_path / "constraints"
    candidate.write_text("a b c\n")
    constraints.write_text("[]\n")
    command = ["python", "replay", "--candidate-clusters", str(candidate),
               "--membership-constraints", str(constraints)]
    plan = {"scientific_source": str(source.parent.parent), "scientific_revision": "original",
        "source_root": str(tmp_path / "runner"), "command": command, "aligner": "mafft",
        "tree_builder": "FastTree", "fasta": "/complete-fasta"}
    report = {"execution": {"command": command, "returncode": 0}, "graph_execution": {
        "command": [sys.executable, str(tmp_path / "runner/benchmark_tools/continue_cap400_downstream_worker.py"),
                    "--execution-plan", str(tmp_path / "execution_plan.json")], "returncode": 0}}
    metrics = {"status": "complete", "source": record(source),
        "git": {"commit": "original", "dirty": False},
        "parameters": {"species_tree_mode": "infer", "species_tree": None,
            "aligner": "mafft", "tree_builder": "FastTree", "root_duplication_rule": "species_overlap",
            "pair_orthology_rule": "positive_paralogy", "species_tree_rooting": "min_variance", "cpu": 32,
            "satellite_membership_policy": "high_confidence_pair",
            "explicit_unconstrained_ablation": False, "checkpoint_source": None},
        "input": {"candidate_clusters": record(candidate), "membership_constraints": record(constraints),
                  "fasta_directory": "/complete-fasta"}}
    return plan, report, metrics


def test_completed_original_native_execution_accepted(tmp_path, native):
    audit.require_native_execution(tmp_path, *native)


@pytest.mark.parametrize("stage", ["graph_execution", "execution"])
@pytest.mark.parametrize("change", [{"returncode": 1}, {"command": ["unplanned"]}])
def test_changed_native_or_graph_execution_refused(tmp_path, native, stage, change):
    plan, report, metrics = deepcopy(native)
    report[stage].update(change)
    with pytest.raises(ValueError, match="original command"):
        audit.require_native_execution(tmp_path, plan, report, metrics)


@pytest.mark.parametrize("field,change", [("status", "running"), ("parameters", {"cpu": 16}),
    ("git", {"commit": "other", "dirty": False}), ("source", {})])
def test_changed_or_incomplete_science_refused(tmp_path, native, field, change):
    plan, report, metrics = deepcopy(native)
    metrics[field] = change
    with pytest.raises(ValueError, match="unchanged completed original"):
        audit.require_native_execution(tmp_path, plan, report, metrics)


@pytest.mark.parametrize("field", ["candidate_clusters", "membership_constraints", "fasta_directory"])
def test_changed_native_input_refused(tmp_path, native, field):
    plan, report, metrics = deepcopy(native)
    metrics["input"][field] = "changed"
    with pytest.raises(ValueError, match="different candidate"):
        audit.require_native_execution(tmp_path, plan, report, metrics)


def test_report_writer_creates_new_parent_directory(tmp_path, monkeypatch):
    result = {"status": "verified_fixture", "scores": {}}
    monkeypatch.setattr(audit, "audit", lambda *args: result)
    output = tmp_path / "new/namespace/report.json"
    assert audit.write_report(tmp_path, tmp_path, output) == result
    assert json.loads(output.read_text()) == result


def test_report_writer_preserves_existing_receipt(tmp_path, monkeypatch):
    output = tmp_path / "report.json"
    output.write_text("existing receipt\n")
    monkeypatch.setattr(audit, "audit", lambda *args: pytest.fail("Existing output was reused"))
    with pytest.raises(FileExistsError):
        audit.write_report(tmp_path, tmp_path, output)
    assert output.read_text() == "existing receipt\n"
