from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest

from benchmark_tools import run_copy_branch_ablation as runner
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def original_command():
    return ["python", "old/replay_phylogeny.py", "--fasta-directory", "fasta", "--candidate-clusters", "candidate",
        "--membership-constraints", "constraints", "--output-directory", "old/output", "--json", "old/metrics.json",
        "--cpu", "16", "--species-tree-mode", "infer", "--root-rule", "species_overlap",
        "--pair-rule", "positive_paralogy", "--species-tree-rooting", "min_variance",
        "--satellite-membership-rule", "high_confidence_pair", "--checkpoint-source", "old/cache",
        "--aligner", "mafft", "--tree-builder", "fasttree"]


def test_native_command_changes_only_execution_paths_with_unchanged_inference():
    original = original_command()
    command = runner.native_command(original, Path("new_source"), Path("new"), "new/candidate", "new/constraints", "new_cache")
    assert original == original_command()
    assert command[1] == "new_source/benchmark_tools/replay_phylogeny.py"
    for flag in ("--fasta-directory", "--cpu", "--species-tree-mode", "--root-rule", "--pair-rule",
                 "--species-tree-rooting", "--satellite-membership-rule", "--aligner", "--tree-builder"):
        assert runner.argument(command, flag) == runner.argument(original, flag)
    assert runner.argument(command, "--checkpoint-source") == "new_cache"


@pytest.mark.parametrize("flag,value", [("--species-tree-mode", "supplied"), ("--root-rule", "mapped_event"),
    ("--pair-rule", "lca"), ("--satellite-membership-rule", "high_confidence_path"), ("--cpu", "32")])
def test_changed_native_strategy_is_rejected(flag, value):
    command = original_command()
    command[command.index(flag) + 1] = value
    with pytest.raises(ValueError, match="unchanged own-tree"):
        runner.native_command(command, Path("source"), Path("out"), "candidate", "constraints", "cache")


def test_terminal_parent_requires_real_raw_job_exit_and_resources():
    header = "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n"
    line = "24203|COMPLETED|0:0|01:00:00|bizon|16\n"
    assert runner.terminal_job(header + line, "24203")["State"] == "COMPLETED"
    assert runner.terminal_job(header + line.replace("COMPLETED|0:0", "FAILED|1:0"), "24203")["State"] == "FAILED"
    for bad in (line.replace("COMPLETED", "RUNNING"), line.replace("0:0", "137:0"),
                line.replace("bizon", "dgx"), line.replace("|16", "|32"), line + line,
                line.replace("24203", "other")):
        with pytest.raises(ValueError):
            runner.terminal_job(header + bad, "24203")


def equivalence_fixture(tmp_path):
    old, new = tmp_path / "old", tmp_path / "new"
    files = {"orthohmm/phylogeny.py": "def reconcile():\n    return 1\n",
        "orthohmm/orthohmm.py": "def _expand_phylogeny_candidates():\n    return 1\n",
        "orthohmm/refinement.py": "DEFAULT_CEILING = 500\ndef helper():\n    return 1\ndef refine_cluster_indices():\n    return 1\ndef _split_high_duplication_index_clusters():\n    return 1\n",
        "benchmark_tools/replay_phylogeny.py": "pass\n"}
    for root in (old, new):
        for name, content in files.items():
            path = root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content)
    return old, new


@pytest.mark.parametrize("mutation", ["downstream", "inventory", "expansion", "helper", "constant"])
def test_reuse_equivalence_rejects_changed_downstream_logic_and_defaults(tmp_path, mutation):
    old, new = equivalence_fixture(tmp_path)
    assert runner.replay_equivalence(old, new)["status"] == "downstream_native_bytes_and_candidate_helpers_unchanged"
    if mutation == "inventory":
        (new / "orthohmm/added.py").write_text("pass\n")
    else:
        filename = {"downstream": "phylogeny.py", "expansion": "orthohmm.py",
                    "helper": "refinement.py", "constant": "refinement.py"}[mutation]
        path = new / "orthohmm" / filename
        content = path.read_text().replace("500", "600") if mutation == "constant" else path.read_text().replace("return 1", "return 2", 1)
        path.write_text(content)
    with pytest.raises(ValueError):
        runner.replay_equivalence(old, new)


def candidate(seed, target, *args):
    target.mkdir(parents=True)
    path, constraints = target / "groups.txt", target / "constraints.json"
    path.write_bytes(seed.read_bytes())
    constraints.write_text(json.dumps({"groups": seed.read_text()}))
    return {"seed_partition": record(seed), "candidate_partition": record(path), "membership_constraints": record(constraints),
            "expansion": {"parameters": {"unchanged": True}}}


def test_construction_proves_narrow_branch_reuse_and_pins_actual_trace(tmp_path, monkeypatch):
    names, species = ["a", "b"], np.array([0, 1], dtype=np.int32)
    empty = np.array([], dtype=np.int32)
    monkeypatch.setattr(runner, "load_graph", lambda *a: (names, species, (empty, empty, np.array([])),
        [[0, 1]], (empty, empty, np.array([])), None, {}, {}))
    def refine(clusters, edges, species, names, profile, seed, **kwargs):
        groups = [[0], [1]] if profile in ("default", "phylogeny_deferred_large") else [[0, 1]]
        seed.write_text("".join(" ".join(names[g] for g in group) + "\n" for group in groups))
        return groups
    monkeypatch.setattr(runner, "regenerate_refinement", refine)
    monkeypatch.setattr(runner, "prepare_partition", candidate)
    monkeypatch.setattr(runner, "refinement_changes", lambda *a: {})
    parents = []
    for index, content in enumerate(("a\nb\n", "a b\n")):
        seed = tmp_path / f"parent_seed{index}"
        seed.write_text(content)
        result = candidate(seed, tmp_path / f"parent{index}")
        parents.append({"candidate": result["candidate_partition"], "constraints": result["membership_constraints"]})
    output = tmp_path / "new"
    output.mkdir()
    rows, _ = runner.construct({"dataset": "OrthoBench", "data": {}, "parent_candidates": parents}, {}, output)
    for arm, reference in (("large_only500", "control"), ("medium_only500", "deferred500")):
        assert rows[arm]["candidate"]["candidate_partition"]["sha256"] == rows[reference]["candidate"]["candidate_partition"]["sha256"]
        assert Path(rows[arm]["copy_split_trace"]["path"]).is_file()


def setup_worker(tmp_path, monkeypatch, dataset="OrthoBench"):
    source = Path(runner.__file__).resolve().parent.parent
    common = tmp_path / "pin.txt"
    common.write_text("immutable\n")
    candidate_info = {k: record(common) for k in ("candidate_partition", "membership_constraints")}
    parent_rows = [{"label": label, "scores": {"OrthoBench F1": score}}
                   for label, score in zip(runner.PARENT_LABELS, (.74, .73))]
    parent_cases = [{"label": label, "output": str(tmp_path / "old" / label)} for label in runner.PARENT_LABELS]
    arms = [{"label": label, "profile": profile, "output": str(tmp_path / label), "command": [], "scoring_id": f"short{index}"}
            for index, (label, profile) in enumerate(runner.ARMS)]
    rows = {label: {"label": label, "profile": profile, "status": "candidates_verified_unscored", "scores": {},
                   "seed": record(common), "copy_split_trace": record(common), "candidate": candidate_info}
            for label, profile in (("control", "default"), ("deferred500", "phylogeny_deferred"), *runner.ARMS)}
    plan = {"source_root": str(source), "repo": str(source), "source_revision": "fixture", "output": str(tmp_path),
            "dataset": dataset, "arms": arms, "parent_cases": parent_cases, "tool_path": "/usr/bin", "scope": "test"}
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(plan))
    monkeypatch.setattr(runner, "verify", lambda *a: None)
    monkeypatch.setattr(runner, "host_activity", lambda: "shared-host fixture")
    monkeypatch.setattr(runner, "resolve_parent", lambda *a: ({"parent_rows": parent_rows}, [record(common)]))
    monkeypatch.setattr(runner, "construct", lambda *a: (rows, {}))
    native = {"genes": runner.EXPECTED_GENES[dataset], "species": 78 if dataset == "QfO" else 12,
        "groups": 1, "pairs": 1, "products": {}, "checked": [], "summary": {"species_tree_families": 2}}
    monkeypatch.setattr(runner, "replay_native", lambda *a: native)
    calls = []
    def readback(repo, dataset, row, native, folder, participant):
        calls.append(("readback", row["label"]))
        return row["scores"], {}, []
    monkeypatch.setattr(runner, "readback_scores", readback)
    return path, calls, plan


def test_narrow_worker_reuses_both_exact_controls_without_native_or_new_scoring(tmp_path, monkeypatch):
    path, calls, _ = setup_worker(tmp_path, monkeypatch)
    monkeypatch.setattr(runner, "timed_stage", lambda *a: pytest.fail("No native inference may run"))
    monkeypatch.setattr(runner, "score_qfo", lambda *a: pytest.fail("No new scoring namespace needed"))
    result = runner.run(path)
    assert result["status"] == "finished_all_branch_arms_verified"
    assert calls == [("readback", "control"), ("readback", "deferred500")]
    assert [r["scores"] for r in result["arms"]] == [{"OrthoBench F1": .74}, {"OrthoBench F1": .73}]
    assert all(r["execution_kind"] == "exact_candidate_structural_reuse_no_native_rerun" for r in result["arms"])
    with pytest.raises(ValueError, match="Existing attempt"):
        runner.run(path)


def test_live_parent_gate_prevents_any_candidate_or_native_work(tmp_path, monkeypatch):
    path, calls, _ = setup_worker(tmp_path, monkeypatch, "QfO")
    def fail(*a):
        raise ValueError("Parent RUNNING")
    monkeypatch.setattr(runner, "resolve_parent", fail)
    monkeypatch.setattr(runner, "construct", lambda *a: pytest.fail("No candidates before parent gate"))
    result = runner.run(path)
    assert result["status"] == "failed_with_explicit_gaps" and result["arms"] == [] and calls == []


def test_qfo_native_failure_is_retained_without_second_arm_or_scoring(tmp_path, monkeypatch):
    path, calls, _ = setup_worker(tmp_path, monkeypatch, "QfO")
    monkeypatch.setattr(runner, "timed_stage", lambda *a: {"returncode": 1})
    monkeypatch.setattr(runner, "score_qfo", lambda *a: pytest.fail("Native failed; scoring forbidden"))
    result = runner.run(path)
    assert result["status"] == "failed_with_explicit_gaps" and len(result["arms"]) == 1
    assert result["arms"][0]["status"] == "failed" and calls == []


@pytest.mark.parametrize("markers", [None, 0, -1, True, 1.5, "2"])
def test_new_own_tree_requires_positive_integer_marker_count(markers):
    with pytest.raises(ValueError, match="markers"):
        runner.require_markers({"species_tree_families": markers})


def test_verify_rejects_changed_inputs_and_cuda_visibility(tmp_path, monkeypatch):
    path = tmp_path / "input.txt"
    path.write_text("old\n")
    plan = {"inputs": [record(path)], "sources": [], "runtime": {"files": []}}
    monkeypatch.setattr(runner, "runtime", lambda: plan["runtime"])
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    runner.verify(plan)
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "0")
    with pytest.raises(ValueError, match="CUDA"):
        runner.verify(plan)
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    path.write_text("new\n")
    with pytest.raises(ValueError):
        runner.verify(plan)


def resolution_fixture(tmp_path, monkeypatch, dataset="QfO"):
    root, output = tmp_path / "parent", tmp_path / "worker"
    root.mkdir()
    output.mkdir()
    parent_path = root / "execution_plan.json"
    parent_path.write_text(json.dumps({"source_revision": "parent_source"}))
    rows = [{"label": label, "status": "independent_content_readback_verified"} for label in runner.PARENT_LABELS]
    top = {"source_revision": "parent_source", "execution_plan": record(parent_path), "job_id": "24203",
           "status": "finished_all_arms_verified", "arms": rows}
    top_path = root / "results.json"
    top_path.write_text(json.dumps(top))
    if dataset == "QfO":
        child_path = top_path
    else:
        child_path = root / dataset / "results.json"
        child_path.parent.mkdir()
        child_path.write_text(json.dumps({"dataset": dataset, "status": "finished_all_arms_verified", "arms": rows}))
    plan = {"parent_job": "24203", "parent_root": str(root), "parent_plan": record(parent_path),
            "parent_results": str(child_path), "dataset": dataset, "repo": str(tmp_path)}
    cells = [{"status": "independent_content_readback_verified", "checked": []} for _ in rows]
    monkeypatch.setattr(runner.subprocess, "check_output", lambda *a, **k:
        "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n24203|COMPLETED|0:0|01:00:00|bizon|16\n")
    monkeypatch.setattr(runner, "collect_deferred_cells", lambda *a: (cells, []))
    return plan, output, top_path, top, cells


def test_parent_resolution_freezes_independently_admitted_controls(tmp_path, monkeypatch):
    plan, output, top_path, _, _ = resolution_fixture(tmp_path, monkeypatch)
    resolution, pins = runner.resolve_parent(plan, output)
    assert resolution["status"] == "terminal_matched_parent_independently_admitted"
    assert pins[0] == record(output / "resolved_parent_inputs.json")
    assert record(top_path) in pins and resolution["scheduler"]["State"] == "COMPLETED"


@pytest.mark.parametrize("mutation", ["live_job", "running_report", "wrong_job", "wrong_revision", "bad_control", "readback"])
def test_parent_resolution_cannot_bypass_terminal_control_gates(tmp_path, monkeypatch, mutation):
    plan, output, top_path, top, cells = resolution_fixture(tmp_path, monkeypatch)
    if mutation == "live_job":
        monkeypatch.setattr(runner.subprocess, "check_output", lambda *a, **k:
            "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n24203|RUNNING|0:0|01:00:00|bizon|16\n")
    elif mutation == "running_report":
        top["status"] = "running"
    elif mutation == "wrong_job":
        top["job_id"] = "other"
    elif mutation == "wrong_revision":
        top["source_revision"] = "other"
    elif mutation == "bad_control":
        top["arms"][0]["status"] = "failed"
    else:
        cells[0]["status"] = "readback_failed"
    top_path.write_text(json.dumps(top))
    with pytest.raises(ValueError):
        runner.resolve_parent(plan, output)


def test_terminal_specific_qfo_variant_failure_is_preserved_but_live_gap_is_not(tmp_path, monkeypatch):
    plan, output, top_path, top, cells = resolution_fixture(tmp_path, monkeypatch)
    monkeypatch.setattr(runner.subprocess, "check_output", lambda *a, **k:
        "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n24203|FAILED|1:0|01:00:00|bizon|16\n")
    top["status"] = "finished_with_explicit_gaps"
    top["arms"][1].update(status="failed", error="Specific scoring failure")
    cells[1]["status"] = "failed"
    top_path.write_text(json.dumps(top))
    resolution, _ = runner.resolve_parent(plan, output)
    assert resolution["parent_rows"][1]["error"] == "Specific scoring failure"
    top["arms"][1]["status"] = "native_running"
    top_path.write_text(json.dumps(top))
    with pytest.raises(ValueError, match="specific terminal"):
        runner.resolve_parent(plan, output)


def test_successful_tk_child_is_reused_despite_explicit_failed_ob_sibling(tmp_path, monkeypatch):
    plan, output, top_path, top, _ = resolution_fixture(tmp_path, monkeypatch, "ThreeKingdoms")
    monkeypatch.setattr(runner.subprocess, "check_output", lambda *a, **k:
        "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n24203|FAILED|1:0|01:00:00|bizon|16\n")
    top["status"] = "finished_with_explicit_gaps"
    top_path.write_text(json.dumps(top))
    result, _ = runner.resolve_parent(plan, output)
    assert result["scheduler"]["State"] == "FAILED"
    assert all(cell["status"] == "independent_content_readback_verified" for cell in result["cells"])


def test_qfo_arms_run_sequentially_after_parent_admission_and_native_validation(tmp_path, monkeypatch):
    path, calls, _ = setup_worker(tmp_path, monkeypatch, "QfO")
    def stage(command, folder, *a):
        calls.append(("native", folder.name))
        phylo = folder / "output/orthohmm_phylogeny"
        phylo.mkdir(parents=True)
        for name in runner.PRODUCTS:
            (phylo / name).write_text("fixture\n")
        return {"returncode": 0}
    def score(repo, phylo, folder, scoring_id):
        calls.append(("score", folder.name))
        return {"scores": {"VGNC F": .9, "SwissTrees F": .8, "TreeFam-A F": .7, "EC": .96, "GO": .49, "FAS": .76}}
    monkeypatch.setattr(runner, "timed_stage", stage)
    monkeypatch.setattr(runner, "validate_partition", lambda *a: {})
    monkeypatch.setattr(runner, "score_qfo", score)
    result = runner.run(path)
    assert result["status"] == "finished_all_branch_arms_verified"
    assert calls == [("native", "large_only500"), ("score", "large_only500"), ("readback", "large_only500"),
                     ("native", "medium_only500"), ("score", "medium_only500"), ("readback", "medium_only500")]
    assert all(row["native_coverage"]["genes"] == 984137 for row in result["arms"])
