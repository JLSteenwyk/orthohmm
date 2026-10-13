from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import deferred_qfo_inferred as helper
from benchmark_tools import run_deferred_copy_split as runner
from benchmark_tools.run_membership_hillclimb import argument
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def original_command():
    return ["python", "old/replay.py", "--fasta-directory", "fasta", "--candidate-clusters", "candidates",
        "--membership-constraints", "constraints", "--output-directory", "parent/output", "--json", "parent/metrics.json",
        "--cpu", "16", "--species-tree-mode", "supplied", "--species-tree", "baseline.nwk",
        "--root-rule", "species_overlap", "--pair-rule", "positive_paralogy", "--species-tree-rooting", "min_variance",
        "--satellite-membership-rule", "high_confidence_pair", "--checkpoint-source", "old_cache"]


def test_inferred_command_changes_only_tree_handoff_and_execution_paths():
    original = original_command()
    result = helper.inferred_command(original, Path("frozen"), Path("new"), 16)
    assert original == original_command()
    assert argument(result, "--species-tree-mode") == "infer" and "--species-tree" not in result
    assert argument(result, "--checkpoint-source") == "parent/output"
    assert argument(result, "--output-directory") == "new/output"
    assert argument(result, "--json") == "new/metrics.json"
    for flag in ("--fasta-directory", "--candidate-clusters", "--membership-constraints", "--root-rule",
                 "--pair-rule", "--species-tree-rooting", "--satellite-membership-rule"):
        assert argument(result, flag) == argument(original, flag)
    assert result[1] == "frozen/benchmark_tools/replay_phylogeny.py"


@pytest.mark.parametrize("mutate", [
    lambda c: c.__setitem__(c.index("--species-tree-mode") + 1, "infer"),
    lambda c: c.__setitem__(c.index("--satellite-membership-rule") + 1, "high_confidence_path"),
    lambda c: c.extend(["--species-tree", "second.nwk"]),
    lambda c: c.append("--unconstrained-membership")])
def test_unexpected_parent_handoff_rejected(mutate):
    command = original_command()
    mutate(command)
    with pytest.raises(ValueError):
        helper.inferred_command(command, Path("source"), Path("out"), 16)


def completed_parent():
    plan = {"parent_job": "24198", "parent_plan": {"path": "plan", "sha256": "sha", "bytes": 1}}
    report = {"job_id": "24198", "execution_plan": plan["parent_plan"], "status": "finished_all_arms_verified",
              "arms": [{"label": label, "status": "independent_content_readback_verified",
                        "inputs_reverified_unchanged": True} for label in helper.LABELS]}
    return report, plan


def test_terminal_parent_requires_actual_complete_arm_inventory():
    report, plan = completed_parent()
    helper.require_terminal_parent(report, plan)
    for field, value in (("status", "running"), ("job_id", "other"), ("execution_plan", {}), ("arms", [])):
        broken = deepcopy(report)
        broken[field] = value
        with pytest.raises(ValueError, match="terminal parent"):
            helper.require_terminal_parent(broken, plan)
    for field, value in (("status", "native_outputs_verified_scoring_pending"),
                         ("inputs_reverified_unchanged", False), ("label", "different")):
        broken = deepcopy(report)
        broken["arms"][1][field] = value
        with pytest.raises(ValueError, match="terminal parent"):
            helper.require_terminal_parent(broken, plan)


@pytest.mark.parametrize("markers", [0, -1, True, 1.5, "26", None])
def test_inferred_summary_needs_real_markers_and_complete_raw_tree_reuse(markers):
    summary = {"reconciled_families": 12, "checkpoint_hits": 12, "species_tree_families": markers}
    with pytest.raises(ValueError, match="marker-family"):
        helper.require_inferred_summary(summary)


def test_inferred_summary_accepts_marker_checkpoint_reuse_but_not_missing_gene_tree_reuse():
    summary = {"reconciled_families": 12, "checkpoint_hits": 12, "species_tree_families": 26,
               "species_tree_checkpoint_hit": True}
    helper.require_inferred_summary(summary)
    with pytest.raises(ValueError, match="every reconciled"):
        helper.require_inferred_summary({**summary, "checkpoint_hits": 11})


def test_incomplete_parent_is_retained_without_native_inference(tmp_path, monkeypatch):
    source = Path(runner.__file__).resolve().parent.parent
    runtime = {"files": []}
    plan = {"source_root": str(source), "repo": str(source), "output": str(tmp_path), "source_revision": "test",
            "evaluation_mode": "inferred_parent", "inputs": [], "sources": [], "runtime": runtime,
            "tool_path": "/usr/bin", "scope": "test", "cases": [
                {"label": label, "output": str(tmp_path / label), "command": ["--species-tree-mode", "infer"]}
                for label in helper.LABELS]}
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(plan))
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    monkeypatch.setattr(runner, "runtime_records", lambda: runtime)
    monkeypatch.setattr(runner, "host_activity", lambda: {})
    def failed(*args):
        raise ValueError("Parent incomplete")
    monkeypatch.setattr(helper, "resolve", failed)
    def native(*args):
        pytest.fail("Incomplete parent cannot enter inference")
    monkeypatch.setattr(runner, "timed_stage", native)
    result = runner.run(path)
    assert result["status"] == "failed_preparation_or_input_gate" and result["error"] == "Parent incomplete"
    assert result["arms"] == [] and result["defaults_promoted"] is False


def test_inferred_mode_rejects_a_supplied_tree_before_any_resolution(tmp_path):
    source = Path(runner.__file__).resolve().parent.parent
    plan = {"source_root": str(source), "repo": str(source), "output": str(tmp_path),
            "evaluation_mode": "inferred_parent", "cases": [
                {"label": label, "command": ["--species-tree-mode", "infer", "--species-tree", "baseline"]}
                for label in helper.LABELS]}
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(plan))
    with pytest.raises(ValueError, match="without a supplied tree"):
        runner.run(path)
    assert not (tmp_path / "results.json").exists()


@pytest.mark.parametrize("cpu,parent_job", [(0, "24198"), (17, "24198"), (True, "24198"), (16, "other")])
def test_preparation_rejects_invalid_budget_or_parent_identity(tmp_path, cpu, parent_job):
    with pytest.raises(ValueError, match="explicit numeric parent"):
        helper.prepare(tmp_path, tmp_path, tmp_path / "out", tmp_path / "parent_plan", tmp_path / "parent_results", cpu, parent_job)


def test_completed_parent_resolution_pins_inputs_without_copying_scores(tmp_path, monkeypatch):
    scores = {"VGNC F": .9, "SwissTrees F": .8, "TreeFam-A F": .6, "EC": .96, "GO": .49, "FAS": .76}
    old_cases, new_cases, rows, natives = [], [], [], {}
    for label in helper.LABELS:
        directory = tmp_path / "parent" / label
        directory.mkdir(parents=True)
        candidate = {}
        for key in ("seed_partition", "candidate_partition", "membership_constraints"):
            path = directory / key
            path.write_text("fixture\n")
            candidate[key] = record(path)
        phylo = directory / "output/orthohmm_phylogeny"
        for name in ("alignments", "checkpoints", "gene_trees"):
            (phylo / name).mkdir(parents=True)
            (phylo / name / "input").write_text("immutable input\n")
        old_command = ["--candidate-clusters", candidate["candidate_partition"]["path"],
                       "--membership-constraints", candidate["membership_constraints"]["path"],
                       "--output-directory", str(directory / "output")]
        old_cases.append({"label": label, "output": str(directory), "command": old_command, "scoring_id": label})
        new_cases.append({"label": label, "command": [*old_command, "--checkpoint-source", str(directory / "output")]})
        rows.append({"label": label, "status": "independent_content_readback_verified", "inputs_reverified_unchanged": True,
                     "seed": candidate["seed_partition"], "candidate": candidate, "scores": scores})
        natives[label] = {"genes": 984137, "species": 78, "phylogeny": str(phylo), "summary": {}, "checked": []}
    parent_plan = tmp_path / "parent_plan.json"
    parent_plan.write_text(json.dumps({"cases": old_cases}))
    parent_results = tmp_path / "parent_results.json"
    parent_results.write_text(json.dumps({"status": "finished_all_arms_verified", "job_id": "24198",
        "execution_plan": record(parent_plan), "arms": rows}))
    plan = {"parent_results": str(parent_results), "parent_plan": record(parent_plan), "parent_job": "24198",
            "repo": str(tmp_path), "cases": new_cases, "expected_control_scores": scores}
    monkeypatch.setattr(helper.subprocess, "check_output", lambda *a, **k:
        "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n24198|COMPLETED|0:0|00:01:00|bizon|16\n")
    monkeypatch.setattr(helper, "replay_native", lambda row, case: natives[row["label"]])
    monkeypatch.setattr(helper, "readback_scores", lambda *a, **k: (scores, {}, []))
    output = tmp_path / "followup"
    output.mkdir()
    result = helper.resolve(plan, output)
    assert [r["label"] for r in result["arms"]] == list(helper.LABELS)
    assert all(r["scores"] == {} for r in result["arms"])
    assert all(r["scores"] == scores for r in result["readbacks"])
    assert len({r["path"] for r in result["checked"]}) == len(result["checked"])
    assert record(parent_results) in result["checked"]
    assert len([r for r in result["checked"] if r["path"].endswith("/input")]) == 6
    assert record(output / "resolved_parent_inputs.json") == result["receipt"]
