import json
from pathlib import Path

import numpy as np
import pytest

from benchmark_tools import run_three_kingdoms_membership_supplement as runner


def test_species_labels_need_not_follow_sorted_fasta_order():
    result = runner.validate_ownership(["a", "b", "c"], np.array([4, 4, 2], dtype=np.int32),
                                       {"a": "first", "b": "first", "c": "second"})
    assert result == {"genes": 3, "species": 2, "classes": {4: "first", 2: "second"}}


@pytest.mark.parametrize("names,species,owners", [
    (["a", "a"], [0, 0], {"a": "first"}),
    (["a"], [0], {"a": "first", "b": "second"}),
    (["a", "b"], [0, 0], {"a": "first", "b": "second"}),
    (["a", "b"], [0, 1], {"a": "first", "b": "first"}),
    (["a"], [-1], {"a": "first"}),
    (["a", "b"], [0], {"a": "first", "b": "second"}),
])
def test_wrong_checkpoint_or_species_ownership_rejected(names, species, owners):
    with pytest.raises(ValueError):
        runner.validate_ownership(names, np.array(species, dtype=np.int32), owners)


def test_species_dtype_must_remain_native():
    with pytest.raises(ValueError, match="species vector"):
        runner.validate_ownership(["a"], np.array([0.]), {"a": "first"})


def test_candidate_control_requires_both_exact_products(tmp_path):
    candidate, constraints = tmp_path / "candidate", tmp_path / "constraints"
    candidate.write_text("a b\n")
    constraints.write_text("[]\n")
    expected = runner.record(candidate), runner.record(constraints)
    assert all(runner.candidate_control_matches(candidate, constraints, *expected).values())
    for path in (candidate, constraints):
        original = path.read_text()
        path.write_text(original + "\n")
        with pytest.raises(ValueError, match="exact retained-output"):
            runner.candidate_control_matches(candidate, constraints, *expected)
        path.write_text(original)


def test_commands_preserve_hmm_profile_and_native_phylogeny_settings():
    execution = {"checkpoint": "hits", "plan": {"checkpoint_manifest_sha256": "sha"},
        "fasta_directory": "input", "cpu": 16, "aligner": "/tools/mafft", "tree_builder": "FastTree"}
    graph = runner.replay_command(Path("source"), execution, Path("out"), .08)
    native = runner.phylogeny_command(Path("source"), execution, Path("out"), Path("candidates"),
        Path("constraints"), "high_confidence_path", Path("matching_direct"))
    def argument(command, key):
        return command[command.index(key) + 1]
    assert argument(graph, "--cpm-resolution") == "0.08"
    assert argument(graph, "--profile-iterations") == "1"
    assert argument(graph, "--profile-min-species") == "1"
    assert argument(graph, "--checkpoint-sha256") == "sha"
    assert argument(native, "--checkpoint-source") == "matching_direct"
    assert argument(native, "--species-tree-mode") == "infer"
    assert argument(native, "--pair-rule") == "positive_paralogy"
    assert argument(native, "--satellite-membership-rule") == "high_confidence_path"
    assert "--official-benchmark" not in graph + native


def setup_execution(tmp_path, monkeypatch):
    from orthohmm import accuracy, orthohmm as core
    source = tmp_path / "source"
    (source / "benchmark_tools").mkdir(parents=True)
    (source / "orthohmm").mkdir()
    script = source / "benchmark_tools/run_three_kingdoms_membership_supplement.py"
    script.write_text("frozen executor\n")
    monkeypatch.setattr(runner, "__file__", str(script))
    monkeypatch.setattr(core, "__file__", str(source / "orthohmm/orthohmm.py"))
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    metadata = {"genes": 443217, "species": 12, "classes": {}}
    monkeypatch.setattr(runner, "validate_ownership", lambda *_: metadata)
    names = ["a", "b", "c"]
    monkeypatch.setattr(accuracy, "load_accuracy_checkpoint", lambda *_args, **_kwargs: (
        names, np.array([0, 1, 2], dtype=np.int32), np.array([0], dtype=np.int32),
        np.array([1], dtype=np.int32), np.array([.1])))
    monkeypatch.setattr(runner, "audit", lambda *_: {"status": "numeric_checkpoint_verified"})
    monkeypatch.setattr(runner, "runtime_records", lambda: {"files": []})
    monkeypatch.setattr(runner, "host_activity", lambda: "fixture host snapshot")
    baseline = tmp_path / "baseline"
    phylogeny = baseline / "orthohmm_phylogeny"
    phylogeny.mkdir(parents=True)
    golden = {}
    for name in (*runner.PRODUCTS, "species_tree.rooted.nwk"):
        path = phylogeny / name
        path.write_text("golden\n")
        golden[name] = runner.record(path)
    fasta = tmp_path / "input.fasta"
    fasta.write_text(">a\nACDE\n>b\nACDE\n>c\nACDE\n")
    candidate, constraints = baseline / "candidates", baseline / "constraints"
    candidate.write_text("a b\nc\n")
    constraints.write_text("[]\n")
    output = tmp_path / "output"
    output.mkdir()
    report = {"source_root": str(source), "output": str(output), "repo": str(tmp_path),
        "source_revision": "fixture", "tool_path": "fixture_path", "inputs": [],
        "sources": [runner.record(script)], "runtime": {"files": []}, "cpu": 16,
        "fasta_directory": str(tmp_path), "fastas": [runner.record(fasta)],
        "checkpoint": str(tmp_path / "checkpoint"), "baseline_output": str(baseline),
        "baseline_candidate": runner.record(candidate), "baseline_constraints": runner.record(constraints),
        "baseline_candidate_parameters": {"fixed": True}, "golden": golden,
        "aligner": "mafft", "tree_builder": "FastTree",
        "plan": {"checkpoint_manifest_sha256": "sha", "limitations": [], "arms": [
            {"label": "control_direct", "cpm_resolution": .1, "membership_rule": "high_confidence_pair"},
            {"label": "broader_direct", "cpm_resolution": .08, "membership_rule": "high_confidence_pair"},
            {"label": "broader_path", "cpm_resolution": .08, "membership_rule": "high_confidence_path"}]}}
    path = output / "execution_plan.json"
    runner.write_json(path, report)

    def fake_stage(command, folder, name, _source, _env):
        if name == "graph":
            graph = folder / "graph"
            graph.mkdir()
            seed = graph / "orthogroups_profiles_refined.txt"
            seed.write_text("a b\nc\n")
            runner.write_json(folder / "graph_metrics.json", {"stages": [{"output": runner.record(seed)}]})
        else:
            phylogeny = folder / "native/orthohmm_phylogeny"
            working = folder / "native/orthohmm_working_res"
            phylogeny.mkdir(parents=True)
            working.mkdir()
            for name in (*runner.PRODUCTS, "species_tree.rooted.nwk"):
                (phylogeny / name).write_text("golden\n")
            runner.write_json(phylogeny / "reconciliation_summary.json", {"ortholog_pairs": 1})
            runner.write_json(phylogeny / "provenance_manifest.json", {"membership_reconciliation": {}})
            runner.write_json(folder / "native_metrics.json", {"status": "complete"})
            (working / "orthohmm_edges_clustered.txt").write_text("a b\nc\n")
        return {"returncode": 0, "command": command}

    def fake_candidates(seed, folder, *_):
        working = folder / "orthohmm_working_res"
        working.mkdir(parents=True)
        candidate, constraints = working / "candidates", working / "constraints"
        candidate.write_text(seed.read_text())
        constraints.write_text("[]\n")
        return {"candidate_partition": runner.record(candidate), "membership_constraints": runner.record(constraints),
                "expansion": {"parameters": {"fixed": True}}}
    monkeypatch.setattr(runner, "timed_stage", fake_stage)
    monkeypatch.setattr(runner, "prepare_partition", fake_candidates)
    monkeypatch.setattr(runner, "score_three_kingdoms", lambda *_: {"scores": {"Three Kingdoms F1": .87}})
    return path


def test_runner_validates_and_scores_all_three_arms(tmp_path, monkeypatch):
    path = setup_execution(tmp_path, monkeypatch)
    report = runner.run(path)
    assert report["status"] == "finished_all_arms_scored"
    assert all(row["status"] == "scored_development_result" for row in report["arms"])
    assert report["arms"][0]["candidate_golden_matches"]
    assert report["arms"][0]["golden_matches"]
    assert all(row["coverage"]["all_candidate_genes_preserved"] for row in report["arms"])
    assert report["arms"][2]["candidate"] == report["arms"][1]["candidate"]
    assert not report["defaults_promoted"]
    with pytest.raises(FileExistsError, match="preserved"):
        runner.run(path)


def test_candidate_control_failure_prevents_native_and_variants(tmp_path, monkeypatch):
    path = setup_execution(tmp_path, monkeypatch)
    original = runner.prepare_partition
    calls = []
    def changed(*args):
        calls.append(True)
        result = original(*args)
        candidate = Path(result["candidate_partition"]["path"])
        candidate.write_text("a\nb c\n")
        result["candidate_partition"] = runner.record(candidate)
        return result
    monkeypatch.setattr(runner, "prepare_partition", changed)
    monkeypatch.setattr(runner, "score_three_kingdoms", lambda *_: pytest.fail("Must not score failed control"))
    report = runner.run(path)
    assert len(calls) == 1
    assert [row["status"] for row in report["arms"]] == ["failed", "not_run_failed_control", "not_run_failed_control"]
    assert "native" not in report["arms"][0]


def test_matching_direct_failure_suppresses_path_arm(tmp_path, monkeypatch):
    path = setup_execution(tmp_path, monkeypatch)
    original = runner.timed_stage
    def fail_broader(command, folder, name, *args):
        if folder.name == "broader_direct" and name == "native":
            return {"returncode": 1, "command": command}
        return original(command, folder, name, *args)
    monkeypatch.setattr(runner, "timed_stage", fail_broader)
    report = runner.run(path)
    assert [row["status"] for row in report["arms"]] == ["scored_development_result", "failed", "not_run_failed_broader_direct"]


def test_input_change_invalidates_scored_evidence(tmp_path, monkeypatch):
    path = setup_execution(tmp_path, monkeypatch)
    execution = json.loads(path.read_text())
    script = Path(execution["sources"][0]["path"])
    def mutate(*_):
        script.write_text("changed source\n")
        return {"scores": {"Three Kingdoms F1": .87}}
    monkeypatch.setattr(runner, "score_three_kingdoms", mutate)
    with pytest.raises(ValueError, match="identity changed"):
        runner.run(path)
    report = json.loads((path.parent / "results.json").read_text())
    assert report["status"] == "failed"
    assert report["arms"][0]["status"] == "invalidated_by_execution_failure"
