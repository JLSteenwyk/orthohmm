from copy import deepcopy
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

from Bio import SeqIO
import numpy as np
import pytest

from benchmark_tools import run_profile_weighting_followup as followup
from benchmark_tools import run_profile_weighting_pilot as pilot
from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_canonical_downstream import commands
from benchmark_tools.run_membership_hillclimb import PRODUCTS, compare_golden, validate_partition, write_json
from benchmark_tools.run_three_kingdoms_membership_supplement import validate_ownership


def pilot_fixture(tmp_path, monkeypatch):
    cells, compact = [], []
    base = {"f1": .7, "precision": .8, "recall": .6}
    for index, mode in enumerate(pilot.ARMS):
        point = {**base, "f1": .7 + (0 if index == 0 else -.01 if index == 3 else .001)}
        gate = pilot.advancement(base, point) if index else None
        cells.append({"label": mode, "status": "independent_content_readback_verified", "source_revision": followup.REVISION,
            "coverage": {"genes": 251378, "species": 12}, "scores": {"OrthoBench: OrthoBench F1": point["f1"]},
            "details": {"precision": point["precision"], "recall": point["recall"]}, "advancement": gate})
        compact.append({"label": mode, **point, "advances": bool(gate and gate["eligible"])})
    plan_path = tmp_path / "pilot_plan.json"
    plan_path.write_text(json.dumps({"source_revision": followup.REVISION}))
    readback = {"status": "profile_pilot_independently_read_back", "accuracy_result_admitted": True,
        "source_revision": followup.REVISION, "cells": cells, "advancing_labels": list(followup.ARMS[1:]),
        "execution_plan": record(plan_path), "producer": record(plan_path), "checked": []}
    path = tmp_path / "readback.json"
    path.write_text(json.dumps(readback))
    admission = {"status": "profile_weighting_orthobench_independently_admitted", "scientific_revision": followup.REVISION,
        "readback": record(path), "score_table": record(plan_path), "arms": compact, "advancing_labels": list(followup.ARMS[1:])}
    api = SimpleNamespace(check=check, record=record, advancement=pilot.advancement)
    monkeypatch.setattr(followup, "committed", lambda *args: (admission, record(plan_path)))
    monkeypatch.setattr(followup, "successful_job", lambda *args: {"fixture": list(args)})
    return api, admission, readback, path


def test_both_passing_arms_are_required_not_only_the_largest_gain(tmp_path, monkeypatch):
    api, admission, _, _ = pilot_fixture(tmp_path, monkeypatch)
    result, readback, plan, pins, accounting = followup.pilot_inputs(api, tmp_path)
    assert result["advancing_labels"] == ["henikoff", "henikoff_cluster"]
    assert len(readback["cells"]) == 4 and plan["source_revision"] == followup.REVISION
    assert pins and len(accounting) == 2


@pytest.mark.parametrize("kind", ["partial", "source", "missing_arm", "score", "gate", "chosen_best_only", "compact_eligibility", "universe"])
def test_pilot_prerequisites_reject_cherry_picking_or_changed_evidence(tmp_path, monkeypatch, kind):
    api, admission, report, path = pilot_fixture(tmp_path, monkeypatch)
    if kind == "partial":
        report["status"] = "running"
    elif kind == "source":
        report["source_revision"] = "other"
    elif kind == "missing_arm":
        report["cells"].pop()
    elif kind == "score":
        admission["arms"][1]["f1"] = .9
    elif kind == "gate":
        report["cells"][1]["advancement"]["eligible"] = False
    elif kind == "chosen_best_only":
        report["advancing_labels"] = ["henikoff"]
    elif kind == "compact_eligibility":
        admission["arms"][3]["advances"] = True
    else:
        report["cells"][1]["coverage"]["genes"] -= 1
    path.write_text(json.dumps(report))
    admission["readback"] = record(path)
    with pytest.raises(ValueError):
        followup.pilot_inputs(api, tmp_path)


def test_committed_prerequisite_cannot_be_replaced_by_current_worktree_content(tmp_path, monkeypatch):
    path = tmp_path / "evidence.json"
    original = b'{"status":"admitted"}'
    path.write_bytes(original)
    monkeypatch.setattr(followup.subprocess, "check_output", lambda *args, **kwargs: original)
    api = SimpleNamespace(record=record)
    assert followup.committed(api, tmp_path, "evidence.json")[0]["status"] == "admitted"
    path.write_text('{"status":"changed"}')
    with pytest.raises(ValueError, match="committed prerequisite"):
        followup.committed(api, tmp_path, "evidence.json")


def execution_fixture(tmp_path, monkeypatch, *, failed=None, candidate_mismatch=False, dataset="ThreeKingdoms"):
    source = tmp_path / "science"
    (source / "benchmark_tools").mkdir(parents=True)
    slurm = source / "benchmark_tools/run_profile_weighting_pilot.slurm"
    slurm.write_text("original science\n")
    output = tmp_path / "output"
    output.mkdir()
    fasta = tmp_path / "input"
    fasta.mkdir()
    (fasta / "sp.faa").write_text(">a\nAAAA\n>b\nAAAA\n>c\nAAAA\n")
    golden = tmp_path / "control/orthohmm_phylogeny"
    golden.mkdir(parents=True)
    for name in (*PRODUCTS, "species_tree.rooted.nwk"):
        (golden / name).write_text("native fixture\n")
    control_candidate = tmp_path / "control_candidates.txt"
    control_candidate.write_text("a b c\n")
    control_constraints = tmp_path / "control_constraints.json"
    control_constraints.write_text("[]\n")
    checkpoint = tmp_path / "checkpoint"
    checkpoint.mkdir()
    (checkpoint / "manifest.json").write_text("fixture\n")
    extra = tmp_path / "admitted.json"
    extra.write_text("fixture\n")
    case = {"dataset": dataset, "checkpoint": str(checkpoint), "checkpoint_manifest": record(checkpoint / "manifest.json"),
        "checkpoint_sha": record(checkpoint / "manifest.json")["sha256"],
        "case": {"fasta_directory": str(fasta), "fastas": [record(fasta / "sp.faa")],
                 "baseline_output": str(golden.parent), "aligner": "mafft", "tree_builder": "FastTree"},
        "control": {"coverage": {"genes": 3, "species": 1}, "scores": {dataset + ": score": .7}},
        "golden": {name: record(golden / name) for name in (*PRODUCTS, "species_tree.rooted.nwk")},
        "candidate": {"candidate_partition": record(control_candidate), "membership_constraints": record(control_constraints)}}
    sources = [record(slurm)]
    runtime = {"files": []}
    inputs = [record(extra), case["checkpoint_manifest"], *case["case"]["fastas"], *case["golden"].values(),
              *case["candidate"].values(), *sources]
    driver = Path(followup.__file__).resolve().parents[1]
    plan = {"status": "prepared_weighting_followup_unrun", "repo": str(tmp_path), "output": str(output),
        "scientific_revision": followup.REVISION, "scientific_source": str(source), "driver_root": str(driver),
        "driver_revision": "fixture", "driver_sources": [], "cpu": 16, "dataset": dataset, "arms": list(followup.ARMS),
        "scoring_ids": {mode: followup.scoring_id(output, mode) for mode in followup.ARMS[1:]}, "case": case,
        "pilot_admission": record(extra), "pilot_readback": record(extra), "pilot_plan": record(extra),
        "prerequisite_accounting": [], "scientific_sources": sources, "runtime": runtime, "inputs": inputs,
        "tool_path": "/tools", "scope": "mocked driver control-flow test"}
    monkeypatch.setattr(followup, "UNIVERSES", {"QfO": (3, 1), "ThreeKingdoms": (3, 1)})
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    monkeypatch.setattr(followup.subprocess, "check_output", lambda *args, **kwargs: followup.REVISION + "\n")
    monkeypatch.setattr(followup, "pilot_inputs", lambda *args: (
        {"readback": record(extra)}, {}, {"sources": sources, "runtime": runtime, "tool_path": "/tools"}, [], []))
    monkeypatch.setattr(followup, "selected_case", lambda *args: (case, []))
    api = SimpleNamespace(record=record, check=check, runtime=lambda: runtime, sources=lambda *args: [],
        load_accuracy_checkpoint=lambda *args, **kwargs: (["a", "b", "c"], np.zeros(3, dtype=np.int32), [], [], []),
        SeqIO=SeqIO, validate_ownership=validate_ownership, host_activity=lambda: {}, write_json=write_json,
        commands=commands, _expand_phylogeny_candidates=None, load_membership_constraints=None,
        validate_partition=validate_partition, PRODUCTS=PRODUCTS, compare_golden=compare_golden,
        endpoint_vector=lambda dataset, values: {dataset + ": " + key: value for key, value in values.items()},
        validate_activation=pilot.validate_activation, read_partition=read_partition)
    executed, scored = [], []
    def timed(command, folder, stage, science, env):
        mode = folder.name
        executed.append((mode, stage))
        assert env["NUMBA_NUM_THREADS"] == "16" and env["CUDA_VISIBLE_DEVICES"] == ""
        if stage == "graph":
            if mode == "none":
                assert "--profile-sequence-weighting" not in command
            else:
                assert command[-2:] == ["--profile-sequence-weighting", mode]
            graph = folder / "graph"
            graph.mkdir()
            seed = graph / "orthogroups_profiles_refined.txt"
            seed.write_text("a b c\n")
            metrics = {"parameters": {"profile_sequence_weighting": mode, "profile_iterations": 1,
                "profile_expansion": True, "matrix": "BLOSUM62", "cpm_resolution": .1, "leiden_seed": 4,
                "jackknife_profile_thresholds": False, "jackknife_single_copy_profiles": False,
                "profile_min_species": 1, "singleton_tie_policy": "first_hit"},
                "counts": {"profiles_built": 1, "strict_profile_edges": 1}, "profile_iterations": [{"profiles_built": 1,
                    "profile_effective_counts": [] if mode == "none" else [{"cluster_id": 0, "effective_sequence_count": 3.}]}],
                "stages": [{"output": record(seed)}]}
            write_json(folder / "graph_metrics.json", metrics)
        else:
            phylo = folder / "native/orthohmm_phylogeny"
            phylo.mkdir(parents=True)
            for name in (*PRODUCTS, "species_tree.rooted.nwk"):
                (phylo / name).write_text("native fixture\n")
            target = folder / "native/orthohmm_working_res/orthohmm_edges_clustered.txt"
            target.parent.mkdir()
            target.write_text("a b c\n")
        return {"returncode": int((mode, stage) == failed), "command": command}
    api.timed_stage = timed
    def partition(seed, folder, *args):
        working = folder / "orthohmm_working_res"
        working.mkdir(parents=True)
        path = working / "orthohmm_edges_clustered.txt"
        path.write_text("a b\nc\n" if candidate_mismatch and folder.parent.name == "none" else "a b c\n")
        constraints = working / "phylogeny_candidate_merges.json"
        constraints.write_text("[]\n")
        return {"candidate_partition": record(path), "membership_constraints": record(constraints)}
    api.prepare_partition = partition
    api.replay_native = lambda row, command: {"genes": 3, "species": 1, "summary": {}, "checked": [],
        "prediction": str(Path(command["command"][command["command"].index("--output-directory") + 1]) / "orthohmm_working_res/orthohmm_edges_clustered.txt")}
    def controls(api, plan, row, native, folder):
        if failed == ("none", "score"):
            raise ValueError("controlled assessment mismatch")
        row["scores"] = {"score": .7}
        return row["scores"], {}, []
    api.control_scores = controls
    def score(*args):
        identifier = args[-1] if dataset == "QfO" else None
        scored.append(identifier)
        return {"scores": {"score": .71}}
    api.score_qfo = score
    api.score_three_kingdoms = score
    api.readback_scores = lambda repo, dataset, row, *args: (row["scores"], {}, [])
    monkeypatch.setattr(followup, "science_api", lambda source: api)
    path = output / "execution_plan.json"
    write_json(path, plan)
    return plan, path, executed, scored


@pytest.mark.parametrize("dataset", ["QfO", "ThreeKingdoms"])
def test_all_declared_arms_run_in_order_with_distinct_qfo_assessments(tmp_path, monkeypatch, dataset):
    plan, path, executed, scored = execution_fixture(tmp_path, monkeypatch, dataset=dataset)
    result = followup.run(path)
    assert result["status"] == "all_weighting_followup_arms_scored"
    assert executed == [(mode, stage) for mode in followup.ARMS for stage in ("graph", "native")]
    assert set(result["deltas"]) == {"henikoff", "henikoff_cluster"}
    if dataset == "QfO":
        assert scored == list(plan["scoring_ids"].values()) and len(set(scored)) == 2
    else:
        assert scored == [None, None]
    with pytest.raises(ValueError, match="preserve"):
        followup.run(path)


@pytest.mark.parametrize("failure", [("none", "graph"), ("none", "native"), ("none", "score")])
def test_failed_control_blocks_all_active_arms(tmp_path, monkeypatch, failure):
    _, path, executed, _ = execution_fixture(tmp_path, monkeypatch, failed=failure)
    result = followup.run(path)
    assert result["status"] == "finished_with_explicit_gaps"
    assert [row["status"] for row in result["arms"]] == ["failed", "not_run_failed_control_gate", "not_run_failed_control_gate"]
    assert all(mode == "none" for mode, stage in executed) and result["deltas"] == {}


def test_changed_canonical_candidate_control_blocks_variants(tmp_path, monkeypatch):
    _, path, executed, _ = execution_fixture(tmp_path, monkeypatch, candidate_mismatch=True)
    result = followup.run(path)
    assert result["arms"][0]["status"] == "failed" and "candidate control differs" in result["arms"][0]["error"]
    assert executed == [("none", "graph")]


def test_one_weighted_failure_does_not_suppress_the_other_passing_arm(tmp_path, monkeypatch):
    _, path, executed, _ = execution_fixture(tmp_path, monkeypatch, failed=("henikoff", "graph"))
    result = followup.run(path)
    assert result["arms"][1]["status"] == "failed"
    assert result["arms"][2]["status"] == "native_and_score_independently_verified"
    assert ("henikoff_cluster", "native") in executed and set(result["deltas"]) == {"henikoff_cluster"}


@pytest.mark.parametrize("kind", ["source", "arms", "budget", "identifier", "unpinned", "fasta", "checkpoint"])
def test_changed_execution_plan_rejected_before_inference(tmp_path, monkeypatch, kind):
    plan, _, _, _ = execution_fixture(tmp_path, monkeypatch)
    if kind == "source":
        plan["scientific_revision"] = "other"
    elif kind == "arms":
        plan["arms"].pop()
    elif kind == "budget":
        plan["cpu"] = 32
    elif kind == "identifier":
        plan["scoring_ids"]["henikoff"] = "other"
    elif kind == "unpinned":
        plan["inputs"].pop(0)
    elif kind == "fasta":
        plan["case"]["case"]["fasta_directory"] = "/other"
    else:
        plan["case"]["checkpoint_sha"] = "other"
    with pytest.raises(ValueError):
        followup.validate_plan(plan)


def test_scoring_identifiers_are_short_deterministic_and_arm_specific():
    output = Path("/a/very/long/execution/namespace")
    ids = [followup.scoring_id(output, mode) for mode in followup.ARMS[1:]]
    assert len(set(ids)) == 2 and all(len(value) == 15 and value.startswith("pw_") for value in ids)
    assert followup.scoring_id(output, "henikoff") == ids[0]


def test_existing_preparation_and_mixed_scientific_modules_are_not_overwritten(tmp_path):
    with pytest.raises(FileExistsError):
        followup.prepare(tmp_path, tmp_path, tmp_path, "QfO")
    with pytest.raises(ValueError, match="another source"):
        followup.science_api(tmp_path / "different_source")


def test_real_standalone_bridge_resolves_every_scientific_helper_to_original_pilot():
    source = Path(__file__).resolve().parents[2] / "benchmarks/work/profile_weighting_pilot_v1_20261018"
    if not source.is_dir():
        pytest.skip("Requires retained immutable pilot checkout")
    program = """import json, pathlib, runpy, sys
scope = runpy.run_path(sys.argv[1], run_name='standalone_bridge_probe')
source = pathlib.Path(sys.argv[2]).resolve()
api = scope['science_api'](source)
paths = {name: str(pathlib.Path(sys.modules[getattr(api, name).__module__].__file__).resolve())
         for name in ('validate_activation','advancement','commands','readback_scores','prepare_partition','control_scores')}
assert all(pathlib.Path(path).is_relative_to(source) for path in paths.values())
print(json.dumps(paths))
"""
    env = {**os.environ, "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
           "MKL_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "1"}
    completed = subprocess.run([sys.executable, "-c", program, str(Path(followup.__file__).resolve()), str(source)],
        capture_output=True, text=True, env=env, cwd=source, timeout=60)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert len(json.loads(completed.stdout)) == 6
