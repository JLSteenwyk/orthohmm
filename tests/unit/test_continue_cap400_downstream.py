from copy import deepcopy
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from benchmark_tools.continue_cap400_downstream import (
    REVISION, phylogeny_command, require_graph_success, require_original_settings, run, validate_graph_stage,
)
from benchmark_tools.continue_cap400_downstream_worker import UnusedSearchHits, verify_metadata


def test_complete_metadata_order_and_species_classes():
    rows = [{"name": "a", "spp": "s1"}, {"name": "b", "spp": "s2"}, {"name": "c", "spp": "s1"}]
    result = verify_metadata(rows, ["a", "b", "c"], np.array([0, 1, 0], dtype="int32"), 2)
    assert result == {"genes": 3, "species": 2, "classes": {"s1": 0, "s2": 1}}
    with pytest.raises(ValueError):
        verify_metadata(rows, ["b", "a", "c"], np.array([0, 1, 0]), 2)
    with pytest.raises(ValueError):
        verify_metadata(rows, ["a", "b", "c"], np.array([0, 1, 1]), 2)
    with pytest.raises(ValueError):
        verify_metadata(rows, ["a", "b", "c"], np.array([0, 1, 0]), 3)
    with pytest.raises(RuntimeError):
        UnusedSearchHits().values()


def command_plan():
    return {"scientific_source": "/frozen", "output": "/recovery", "fasta": "/fasta",
            "aligner": "/bin/mafft", "tree_builder": "/bin/FastTree"}


def test_original_phylogeny_command_is_explicit_and_has_no_fake_reuse():
    command = phylogeny_command(command_plan())
    assert command[:2] == [sys.executable, "/frozen/benchmark_tools/replay_phylogeny.py"]
    flags = dict(zip(command[2::2], command[3::2]))
    assert flags == {
        "--fasta-directory": "/fasta",
        "--candidate-clusters": "/recovery/graph/orthohmm_working_res/phylogeny_candidate_superfamilies.txt",
        "--membership-constraints": "/recovery/graph/orthohmm_working_res/phylogeny_candidate_merges.json",
        "--output-directory": "/recovery/native", "--json": "/recovery/native_metrics.json",
        "--cpu": "32", "--species-tree-mode": "infer", "--aligner": "/bin/mafft",
        "--tree-builder": "/bin/FastTree", "--root-rule": "species_overlap",
        "--pair-rule": "positive_paralogy", "--species-tree-rooting": "min_variance",
        "--satellite-membership-rule": "high_confidence_pair"}
    assert "--checkpoint-source" not in command and "--accuracy-checkpoint" not in command


def graph_fixture():
    partition = {"path": "/partition", "sha256": "hash", "bytes": 3}
    admission = {"status": "original_frozen_final_partition_independently_read_back_downstream_pending",
                 "scores": {}, "scientific_revision": REVISION, "partition": partition,
                 "genes": 984137, "groups": 276617}
    plan = {"scientific_revision": REVISION, "output": "/graph",
            "metadata": {"cpm_resolution": .1, "seed": 4, "include_isolates": True,
                         "output_directory": "/graph/native"}}
    report = {"status": "final_graph_partition_verified_downstream_pending", "optimizer_attempted": True,
              "postflight_pins_verified": True, "execution": {"returncode": 0}, "scores": {},
              "partition": {"product": deepcopy(partition)}, "independent_rows": {
                  "status": "all_saved_rows_and_float64_bits_verified", "edges": 33585180,
                  "genes": 984137, "isolates_preserved_by_complete_gene_table": True,
                  "edge_rows_deduplicated_or_reordered": False}}
    scheduler = {"JobIDRaw": "24347", "State": "COMPLETED", "ExitCode": "0:0",
                 "AllocCPUS": "1", "NodeList": "bizon"}
    return admission, plan, report, scheduler


def test_accept_exact_graph_success():
    require_graph_success(*graph_fixture())


@pytest.mark.parametrize("index,key,value", [
    (0, "status", "pending"), (0, "scores", {"F1": 1}), (0, "groups", 1),
    (1, "scientific_revision", "new"), (1, "metadata", {}),
    (2, "postflight_pins_verified", False), (2, "optimizer_attempted", False),
    (2, "execution", {"returncode": -11}), (2, "independent_rows", {}),
    (3, "State", "RUNNING"), (3, "ExitCode", "0:11"), (3, "JobIDRaw", "24245"),
])
def test_graph_success_gate_rejects_changes(index, key, value):
    args = list(graph_fixture())
    args[index][key] = value
    with pytest.raises(ValueError):
        require_graph_success(*args)


def test_attempted_namespace_cannot_relaunch(tmp_path, monkeypatch):
    import benchmark_tools.continue_cap400_downstream as runner
    plan = {"source_root": str(Path(runner.__file__).resolve().parent.parent),
            "repo": str(tmp_path), "output": str(tmp_path), "scientific_source": "/frozen"}
    path = tmp_path / "execution_plan.json"
    path.write_text(json.dumps(plan))
    (tmp_path / "results.json").write_text("{}")
    def forbidden(*args, **kwargs):
        pytest.fail("Attempted namespace must not launch subprocesses")
    monkeypatch.setattr(runner.subprocess, "run", forbidden)
    with pytest.raises(ValueError, match="no silent retry"):
        run(path)


def settings_fixture():
    plan = {**command_plan(), "scientific_revision": REVISION, "checkpoint": "/checkpoint",
            "partition": {"path": "/partition"}, "payload": "/graph/payload",
            "runtime": {"original": True}, "tool_path": "/bin", "control": {"scores": {"F1": .5}}}
    plan["command"] = phylogeny_command(plan)
    parent = {"scientific_source": "/frozen", "runtime": plan["runtime"], "tool_path": "/bin",
              "case": {"fresh_case": {"fasta": "/fasta", "command": ["--aligner", "/bin/mafft",
                                       "--tree-builder", "/bin/FastTree"]}}}
    admission = {"complete_gene_order": {"path": "/checkpoint/gene_names.txt"},
                 "partition": {"path": "/partition"}}
    graph_plan = {"output": "/graph"}
    return plan, parent, admission, graph_plan, {"control": deepcopy(plan["control"])}


def test_accept_authenticated_original_settings():
    require_original_settings(*settings_fixture())


@pytest.mark.parametrize("field", ["scientific_source", "scientific_revision", "fasta", "aligner",
    "tree_builder", "checkpoint", "partition", "payload", "runtime", "tool_path", "control", "command"])
def test_reject_changed_original_settings(field):
    args = list(settings_fixture())
    args[0][field] = "changed"
    with pytest.raises(ValueError):
        require_original_settings(*args)


def stage_fixture():
    graph = "/recovery/graph/"
    working = graph + "orthohmm_working_res/"
    paths = [graph + "orthogroups_refined.txt", *[working + n for n in (
        "phylogeny_candidate_superfamilies.txt", "phylogeny_candidate_merges.json",
        "phylogeny_candidate_seeds.tsv", "orthohmm_edges_clustered.txt")]]
    return {"status": "original_refinement_and_satellite_candidates_verified", "execution_plan": {},
        "scientific_revision": REVISION, "ownership": {"genes": 984137, "species": 78},
        "initial_hmm_search_rerun": False, "profile_expansion_rerun": False,
        "graph_clustering_rerun": False, "scores": {}, "candidate": {
            "profile": "satellite_v2", "candidate_checkpoint": paths[1],
            "membership_policy": "high_confidence_pair", "seed_families": 2, "candidate_families": 1},
        "refined_groups": 2, "candidate_groups": 1,
        "products": [{"path": p} for p in paths],
        "modules": [{"path": "/frozen/orthohmm/" + n} for n in ("orthohmm.py", "accuracy.py", "helpers.py")]}


def test_accept_original_graph_stage(monkeypatch):
    import benchmark_tools.continue_cap400_downstream as runner
    checked = []
    monkeypatch.setattr(runner, "check", checked.append)
    validate_graph_stage(stage_fixture(), command_plan(), {})
    assert len(checked) == 8


@pytest.mark.parametrize("change", [{"initial_hmm_search_rerun": True},
    {"profile_expansion_rerun": True}, {"graph_clustering_rerun": True},
    {"scores": {"F1": 1}}, {"scientific_revision": "new"},
    {"candidate_groups": 9}, {"refined_groups": 9}, {"modules": []}, {"products": []}])
def test_graph_stage_rejects_unplanned_behavior(change, monkeypatch):
    import benchmark_tools.continue_cap400_downstream as runner
    monkeypatch.setattr(runner, "check", lambda pin: None)
    stage = stage_fixture()
    stage.update(change)
    with pytest.raises(ValueError):
        validate_graph_stage(stage, command_plan(), {})


def test_original_frozen_refinement_uses_broad_branch(tmp_path):
    source = Path(os.environ.get("ORTHOHMM_RECOVERY_SCIENCE",
        str(Path(__file__).resolve().parents[2] / "benchmarks/work/full_search_hillclimb_r1_v3_20261012")))
    if not (source / "orthohmm/orthohmm.py").is_file():
        pytest.skip("Original scientific worktree required for native integration fixture")
    code = r'''
import pathlib, sys
import numpy as np
sys.path.insert(0, sys.argv[1])
import orthohmm.orthohmm as engine
from orthohmm.helpers import IndexedEdges
assert pathlib.Path(engine.__file__).resolve() == pathlib.Path(sys.argv[1]) / "orthohmm/orthohmm.py"
class NoHits:
    def values(self):
        raise RuntimeError("unexpected hit access")
names = [f"g{i:02d}" for i in range(78)]
rows = np.array([(f"s{i:02d}", n, 100) for i,n in enumerate(names)],
                dtype=[("spp", object), ("name", object), ("len", int)])
output = pathlib.Path(sys.argv[2])
working = output / "orthohmm_working_res"
working.mkdir()
partition = working / "orthohmm_edges_clustered.txt"
partition.write_text(" ".join(names) + "\n")
edges = IndexedEdges(names, np.arange(77, dtype="int32"), np.arange(1,78,dtype="int32"), np.ones(77))
result = engine._refine_cluster_file(str(output), rows, NoHits(), edges, 1e-4, refinement_profile="default")
assert set(partition.read_text().split()) == set(names)
assert set(result) == {"setup_and_cluster_read", "hit_and_edge_collection", "numeric_refinement", "cluster_write"}
try:
    engine._refine_cluster_file(str(output), rows[:49], NoHits(), edges, 1e-4, refinement_profile="default")
except RuntimeError as error:
    assert str(error) == "unexpected hit access"
else:
    raise AssertionError("Narrow branch must not silently accept the sentinel")
print("original broad-panel production refinement verified")
'''
    env = {**os.environ, "PYTHONPATH": str(source), "CUDA_VISIBLE_DEVICES": "",
           "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    done = subprocess.run([sys.executable, "-c", code, str(source), str(tmp_path)],
                          cwd=source, env=env, text=True, capture_output=True, timeout=90)
    assert done.returncode == 0, done.stdout + done.stderr
    assert "original broad-panel production refinement verified" in done.stdout
