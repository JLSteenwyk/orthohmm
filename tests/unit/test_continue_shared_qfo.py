from copy import deepcopy
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest

from benchmark_tools import continue_shared_qfo as continuation
from benchmark_tools.run_full_search_hillclimb import native_command
from benchmark_tools.run_membership_hillclimb import DATASETS


def original(tmp_path):
    source, fasta, folder = tmp_path / "science", tmp_path / "fasta", tmp_path / "original_native"
    command = native_command(source, fasta, folder, ["mafft", "FastTree"], 100, 16)
    command.extend(["--search-result-order", "canonical", "--search-min-total-hits", "1"])
    qfo = {"dataset": "QfO", "label": "shared_hit1", "search_min_total_hits": 1,
           "search_max_candidates_per_query": 100, "command": command, "output": str(folder),
           "fasta": str(fasta), "inputs": [{"path": f"fasta{i}"} for i in range(78)]}
    return {"source_revision": continuation.REVISION, "source_root": str(source), "cpu": 16,
            "cases": [qfo if dataset == "QfO" else {"dataset": dataset} for dataset in DATASETS]}


def skipped_report(plan, identity):
    return {"status": "finished_with_explicit_gaps", "job_id": "24242",
            "source_revision": continuation.REVISION, "execution_plan": identity,
            "arms": [{"dataset": dataset, "label": "shared_hit1", "scores": {},
                      "status": "not_run_unverified_control" if dataset == "QfO" else "other",
                      "error_type": "ValueError",
                      "error": "Canonical control unavailable or not independently verified"} for dataset in DATASETS]}


def test_relocation_changes_only_output_and_metrics_paths(tmp_path):
    plan = original(tmp_path)
    unchanged = deepcopy(plan)
    case = continuation.relocated_case(plan, tmp_path / "new")
    old = continuation.original_case(plan)
    assert case["command"][:3] == old["command"][:3]
    assert case["command"][5:] == old["command"][5:]
    assert case["command"][3:5] == [str(tmp_path / "new/native/output"), str(tmp_path / "new/native/metrics.json")]
    assert {k: v for k, v in case.items() if k not in ("output", "command")} == {
        k: v for k, v in old.items() if k not in ("output", "command")}
    assert plan == unchanged


@pytest.mark.parametrize("field,value", [("search_min_total_hits", 4), ("search_max_candidates_per_query", 400),
                                         ("label", "other"), ("inputs", [])])
def test_changed_original_intervention_rejected(tmp_path, field, value):
    plan = original(tmp_path)
    continuation.original_case(plan)[field] = value
    with pytest.raises(ValueError, match="intervention or inventory"):
        continuation.original_case(plan)


@pytest.mark.parametrize("field,value", [("cpu", 32), ("source_revision", "other")])
def test_wrong_original_scientific_revision_or_resources_rejected(tmp_path, field, value):
    plan = original(tmp_path)
    plan[field] = value
    with pytest.raises(ValueError, match="frozen three-dataset"):
        continuation.original_case(plan)


def test_changed_profile_or_phylogeny_command_rejected(tmp_path):
    plan = original(tmp_path)
    command = continuation.original_case(plan)["command"]
    command[command.index("--species-tree-mode") + 1] = "none"
    with pytest.raises(ValueError, match="prespecified method"):
        continuation.original_case(plan)


def test_only_explicit_unrun_qfo_arm_is_eligible(tmp_path):
    plan = original(tmp_path)
    identity = {"path": "plan", "sha256": "same", "bytes": 1}
    report = skipped_report(plan, identity)
    row = continuation.original_skip(report, plan, identity)
    assert row["dataset"] == "QfO"
    Path(continuation.original_case(plan)["output"]).mkdir()
    with pytest.raises(ValueError, match="do not duplicate"):
        continuation.original_skip(report, plan, identity)


@pytest.mark.parametrize("change", [
    {"status": "failed"}, {"status": "independent_content_and_score_readback_verified"},
    {"error": "Different failure"}, {"error_type": "RuntimeError"}, {"scores": {"GO": .8}},
    {"native": {}}, {"execution": {}}, {"control": {}}, {"cpu_backend_probe": {}}])
def test_inference_attempts_and_other_failures_not_eligible(tmp_path, change):
    plan = original(tmp_path)
    identity = {"path": "plan"}
    report = skipped_report(plan, identity)
    next(row for row in report["arms"] if row["dataset"] == "QfO").update(change)
    with pytest.raises(ValueError, match="do not duplicate"):
        continuation.original_skip(report, plan, identity)


@pytest.mark.parametrize("change", [{"status": "running"}, {"job_id": "other"},
                                    {"execution_plan": {"path": "other"}}])
def test_original_report_must_be_terminal_and_bound_to_plan(tmp_path, change):
    plan = original(tmp_path)
    identity = {"path": "plan"}
    report = skipped_report(plan, identity)
    report.update(change)
    with pytest.raises(ValueError, match="terminal original"):
        continuation.original_skip(report, plan, identity)


@pytest.mark.parametrize("accounting", [
    "24242|FAILED|1:0\n24246|COMPLETED|0:0\n",
    "24246|COMPLETED|0:0\n24242|COMPLETED|0:0\n"])
def test_terminal_dependency_accounting_accepted(accounting):
    continuation.terminal_dependencies(accounting)


@pytest.mark.parametrize("accounting", [
    "24242|RUNNING|0:0\n24246|COMPLETED|0:0\n", "24242|FAILED|1:0\n24246|FAILED|1:0\n",
    "24242|FAILED|1:0\n24246|RUNNING|0:0\n", "24242|CANCELLED|0:15\n24246|COMPLETED|0:0\n",
    "24242|FAILED|1:0\n", "24242|FAILED|1:0\n24242|FAILED|1:0\n", ""])
def test_no_execution_from_live_failed_or_ambiguous_dependencies(accounting):
    with pytest.raises(ValueError):
        continuation.terminal_dependencies(accounting)


def test_scoring_namespace_short_unique_and_collision_guard(tmp_path):
    with TemporaryDirectory(dir="/tmp", prefix="sh-") as temporary:
        repo = Path(temporary)
        name = continuation.short_id(repo, tmp_path)
        assert len(name) == 15
        assert name == continuation.short_id(repo, tmp_path)
        assert name != continuation.short_id(repo, tmp_path / "another")
        (repo / "qfo_benchmark/w" / ("hc1_" + name)).mkdir(parents=True)
        with pytest.raises(FileExistsError):
            continuation.short_id(repo, tmp_path)


def test_worker_cannot_start_native_with_nonterminal_dependency(tmp_path, monkeypatch):
    source = tmp_path / "runner"
    output = tmp_path / "output"
    output.mkdir()
    path = output / "execution_plan.json"
    path.write_text(json.dumps({"repo": str(tmp_path), "source_root": str(source),
        "scientific_source": str(tmp_path / "science"), "output": str(output),
        "parent_root": str(tmp_path / "parent"), "recovery_root": str(tmp_path / "recovery")}))
    monkeypatch.setattr(continuation, "__file__", str(source / "benchmark_tools/continue_shared_qfo.py"))
    monkeypatch.setattr(continuation, "host_activity", lambda: {})
    monkeypatch.setattr(continuation.subprocess, "check_output", lambda *args, **kwargs:
        "24242|RUNNING|0:0\n24246|COMPLETED|0:0\n")
    monkeypatch.setattr(continuation, "timed_stage", lambda *args: pytest.fail("Native inference started"))
    result = continuation.run(path)
    assert result["status"] == "failed_or_unverified_shared_qfo_continuation"
    assert "execution" not in result
    assert result["scores"] == {}
    assert result["goal_completed"] is False
