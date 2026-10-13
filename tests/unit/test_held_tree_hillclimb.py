import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

from benchmark_tools import run_held_tree_hillclimb as runner
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def command():
    return ["python", "replay.py", "--species-tree-mode", "infer", "--output-directory", "old-output",
            "--json", "old-metrics", "--cpu", "32", "--checkpoint-source", "historical-output",
            "--satellite-membership-rule", "high_confidence_pair", "--membership-constraints", "constraints.json",
            "--candidate-clusters", "candidates.txt"]


def test_handoff_changes_only_tree_not_candidates_or_rules(tmp_path):
    original = command()
    own = runner.supplied_command(original, tmp_path, tmp_path / "cell", tmp_path / "own.nwk", 16)
    baseline = runner.supplied_command(original, tmp_path, tmp_path / "cell", tmp_path / "baseline.nwk", 16)
    index = own.index("--species-tree") + 1
    assert own[:index] == baseline[:index]
    assert own[index + 1:] == baseline[index + 1:]
    assert runner.argument(own, "--checkpoint-source") == "old-output"
    for flag in ("--candidate-clusters", "--membership-constraints", "--satellite-membership-rule"):
        assert runner.argument(own, flag) == runner.argument(original, flag)
    assert original == command()


@pytest.mark.parametrize("change", ["mode", "tree", "policy", "unconstrained", "missing", "duplicate"])
def test_unexpected_source_flags_rejected(tmp_path, change):
    argv = command()
    if change == "mode":
        argv[argv.index("--species-tree-mode") + 1] = "supplied"
    elif change == "tree":
        argv.extend(["--species-tree", "tree.nwk"])
    elif change == "policy":
        argv[argv.index("--satellite-membership-rule") + 1] = "high_confidence_path"
    elif change == "unconstrained":
        argv.append("--unconstrained-membership")
    elif change == "missing":
        index = argv.index("--membership-constraints")
        del argv[index:index + 2]
    else:
        argv.extend(["--checkpoint-source", "other"])
    with pytest.raises(ValueError):
        runner.supplied_command(argv, tmp_path, tmp_path / "cell", tmp_path / "tree", 16)


def test_tree_contrast_preserves_rooted_taxon_identity(tmp_path):
    left, right = tmp_path / "left.nwk", tmp_path / "right.nwk"
    left.write_text("((A:1,B:1):1,(C:1,D:1):1);\n")
    right.write_text("((B:2,A:2):2,(D:2,C:2):2);\n")
    assert runner.tree_contrast(left, right)["rooted_clade_distance"] == 0
    right.write_text("((A:1,C:1):1,(B:1,D:1):1);\n")
    assert runner.tree_contrast(left, right)["rooted_clade_distance"] == 4
    right.write_text("((A:1,B:1):1,(C:1,E:1):1);\n")
    with pytest.raises(ValueError, match="taxon universes"):
        runner.tree_contrast(left, right)


@pytest.mark.parametrize("families,hits", [(0, 0), (2, 1), (1, True), (True, 1), (1, 1.0)])
def test_checkpoint_reuse_is_not_assumed(families, hits):
    with pytest.raises(ValueError, match="checkpoint"):
        runner.require_reuse({"reconciled_families": families, "checkpoint_hits": hits})


def test_fas_not_required_to_reproduce_but_other_endpoints_are():
    expected = dict.fromkeys((*runner.DETERMINISTIC, "FAS"), .8)
    runner.require_control_scores({**expected, "FAS": .81}, expected)
    for name in runner.DETERMINISTIC:
        with pytest.raises(ValueError, match="Deterministic"):
            runner.require_control_scores({**expected, name: .81}, expected)


def fake_plan(tmp_path, monkeypatch):
    output = tmp_path / "results"
    output.mkdir()
    plan = {"source_root": str(Path(runner.__file__).resolve().parent.parent), "repo": str(tmp_path),
        "output": str(output), "source_revision": "fixture", "sources": [], "inputs": [],
        "runtime": {"files": []}, "tool_path": os.environ["PATH"], "scope": "fixture",
        "cases": [{"label": label, "output": str(output / label), "command": [label]}
                  for label in runner.LABELS]}
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(plan))
    monkeypatch.setattr(runner, "host_activity", lambda: "fixture")
    return path


def test_failed_control_never_launches_variant_and_cannot_overwrite(tmp_path, monkeypatch):
    path = fake_plan(tmp_path, monkeypatch)
    calls = []
    monkeypatch.setattr(runner, "timed_stage", lambda argv, *args:
                        calls.append(argv) or {"returncode": 1})
    result = runner.run(path)
    assert calls == [[runner.LABELS[0]]]
    assert [row["status"] for row in result["arms"]] == ["failed", "not_run_control_gate_failed"]
    assert result["status"] == "finished_with_explicit_gaps"
    with pytest.raises(ValueError, match="existing results"):
        runner.run(path)


def test_native_control_mismatch_prevents_scoring_and_variant(tmp_path, monkeypatch):
    path = fake_plan(tmp_path, monkeypatch)
    data = json.loads(path.read_text())
    data["golden"] = {}
    tree = tmp_path / "tree.nwk"
    tree.write_text("((A,B),(C,D));")
    for case in data["cases"]:
        case["tree"] = record(tree)
        case["command"] += ["--candidate-clusters", str(tmp_path / "candidate")]
    path.write_text(json.dumps(data))
    monkeypatch.setattr(runner, "timed_stage", lambda *args: {"returncode": 0})
    monkeypatch.setattr(runner, "validate_partition", lambda *args: {})
    monkeypatch.setattr(runner, "record", lambda path: {"path": str(path), "sha256": "fixture"})
    monkeypatch.setattr(runner, "replay_native", lambda *args: {
        "genes": 984137, "species": 78, "summary": {"reconciled_families": 1, "checkpoint_hits": 1}})
    monkeypatch.setattr(runner, "tree_contrast", lambda *args: {"rooted_clade_distance": 0})
    def mismatch(*args):
        raise ValueError("control outputs differ")
    monkeypatch.setattr(runner, "compare_golden", mismatch)
    monkeypatch.setattr(runner, "score_qfo", lambda *args: pytest.fail("Must not score invalid control"))
    # Keep the execution plan identity real while synthetic native products are stubbed.
    monkeypatch.setattr(runner, "check", lambda *args: None)
    result = runner.run(path)
    assert result["arms"][0]["error"] == "control outputs differ"
    assert result["arms"][1]["status"] == "not_run_control_gate_failed"


def test_real_supplied_own_tree_replay_preserves_native_outputs(tmp_path):
    from benchmark_tools.run_full_search_hillclimb import native_command
    from benchmark_tools.verify_frozen_phylogeny_install import fixture
    source = Path(runner.__file__).resolve().parent.parent
    tools = [shutil.which("mafft"), shutil.which("FastTree")]
    if not all(tools):
        pytest.skip("Requires MAFFT and FastTree")
    fasta = tmp_path / "input"
    fasta.mkdir()
    for taxon, genes in fixture().items():
        (fasta / f"{taxon}.faa").write_text("".join(f">{gene}\n{seq}\n" for gene, seq in genes.items()))
    control = tmp_path / "control"
    control.mkdir()
    env = {**os.environ, "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "2",
           "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "PYTHONHASHSEED": "0"}
    child = subprocess.run(native_command(source, fasta, control, tools, 100, 2),
                           cwd=source, env=env, capture_output=True, text=True, timeout=180)
    assert child.returncode == 0, child.stdout + child.stderr
    output = control / "output"
    phylo, working = output / "orthohmm_phylogeny", output / "orthohmm_working_res"
    candidates = working / "phylogeny_candidate_superfamilies.txt"
    original = [str(value) for value in ["python", "replay.py", "--fasta-directory", fasta,
        "--candidate-clusters", candidates, "--membership-constraints", working / "phylogeny_candidate_merges.json",
        "--output-directory", output, "--json", control / "metrics.json", "--cpu", 2,
        "--species-tree-mode", "infer", "--aligner", tools[0], "--tree-builder", tools[1],
        "--root-rule", "species_overlap", "--pair-rule", "positive_paralogy", "--species-tree-rooting", "min_variance",
        "--satellite-membership-rule", "high_confidence_pair", "--checkpoint-source", output]]
    destination = tmp_path / "replayed"
    destination.mkdir()
    golden = {name: record(phylo / name) for name in runner.PRODUCTS}
    raw_inputs = [record(p) for name in ("checkpoints", "gene_trees") for p in (phylo / name).rglob("*") if p.is_file()]
    argv = runner.supplied_command(original, source, destination, phylo / "species_tree.rooted.nwk", 2)
    child = subprocess.run(argv, cwd=source, env=env, capture_output=True, text=True, timeout=180)
    assert child.returncode == 0, child.stdout + child.stderr
    new_phylo = destination / "output/orthohmm_phylogeny"
    assert all(runner.compare_golden(new_phylo, golden).values())
    runner.require_reuse(json.loads((new_phylo / "reconciliation_summary.json").read_text()))
    assert runner.tree_contrast(phylo / "species_tree.rooted.nwk", new_phylo / "species_tree.rooted.nwk")["rooted_clade_distance"] == 0
    row = {"outputs": [record(new_phylo / name) for name in runner.PRODUCTS],
        "coverage": runner.validate_partition(candidates, destination / "output/orthohmm_working_res/orthohmm_edges_clustered.txt")}
    native = runner.replay_native(row, {"command": argv})
    assert (native["genes"], native["species"]) == (16, 4)
    for item in raw_inputs:
        runner.check(item)
