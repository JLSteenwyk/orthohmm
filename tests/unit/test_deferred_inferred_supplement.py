from pathlib import Path
import json

import numpy as np
import pytest

from benchmark_tools import run_deferred_inferred_supplement as runner
from benchmark_tools.run_deferred_copy_split import regenerate_refinement
from orthohmm.helpers import IndexedEdges, write_indexed_edges


def test_native_writer_reader_roundtrip_is_full_precision_and_indexed(tmp_path):
    names = ["z", "alpha", "gene_with_a_long_name"]
    weights = np.array([np.nextafter(1., 2.), .013465874648943589, 2.7])
    original = IndexedEdges(names, np.array([0, 0, 1], dtype=np.int32),
                           np.array([1, 2, 2], dtype=np.int32), weights)
    (tmp_path / "orthohmm_working_res").mkdir()
    write_indexed_edges(original, str(tmp_path))
    result = runner.load_native_edges(tmp_path / "orthohmm_working_res/orthohmm_edges.txt", names, 3)
    assert result.gene_names == names
    for left, right in ((original.sources, result.sources), (original.targets, result.targets), (weights, result.weights)):
        assert np.array_equal(left, right) and left.dtype == right.dtype


@pytest.mark.parametrize("content,count", [
    ("a\tb\t1\n", 2), ("a\tb\t1\na\tc\t2\n", 1),
    ("a\tb\t1\na\tb\t1\n", 2), ("a\tc\t1\na\tb\t1\n", 2),
    ("b\ta\t1\n", 1), ("a\ta\t1\n", 1), ("a\tx\t1\n", 1),
    ("a\tb\tnan\n", 1), ("a\tb\tinf\n", 1), ("a\tb\t0\n", 1),
    ("a\tb\t-1\n", 1), ("a\tb\tinvalid\n", 1), ("a b 1\n", 1),
    ("a\tb\t1\textra\n", 1)])
def test_invalid_native_graphs_are_not_silently_repaired(tmp_path, content, count):
    path = tmp_path / "edges"
    path.write_text(content)
    with pytest.raises(ValueError):
        runner.load_native_edges(path, ["a", "b", "c"], count)


@pytest.mark.parametrize("count,names", [(0, ["a"]), (-1, ["a"]), (True, ["a"]), (1.5, ["a"]), (1, ["a", "a"])])
def test_invalid_native_dimensions_rejected_before_reading(tmp_path, count, names):
    with pytest.raises(ValueError):
        runner.load_native_edges(tmp_path / "absent", names, count)


def test_refinement_helper_passes_native_narrow_panel_hit_arrays(tmp_path, monkeypatch):
    from orthohmm import refinement
    queries, targets, scores = np.array([0], dtype=np.int32), np.array([1], dtype=np.int32), np.array([2.])
    species = np.array([0, 1], dtype=np.int32)
    empty = np.array([], dtype=np.int32)
    captured = []
    def refine(clusters, q, t, s, *_args, **kwargs):
        captured.append((q, t, s, kwargs))
        return clusters
    monkeypatch.setattr(refinement, "refine_cluster_indices", refine)
    hits = runner.production_refinement_hits(queries, targets, scores, species)
    result = regenerate_refinement([[0], [1]], (empty, empty, np.array([])), species, ["a", "b"],
                                   "phylogeny_deferred", tmp_path / "seed", hits=hits)
    assert result == [[0], [1]]
    assert all(actual is expected for actual, expected in zip(captured[0][:3], hits))
    assert captured[0][3]["defer_copy_split_max_size"] == 500


def test_lossless_graph_can_reproduce_native_leiden_partition(tmp_path):
    from orthohmm.externals import execute_leiden
    from orthohmm.accuracy import read_index_clusters
    original = IndexedEdges(["a", "b", "c", "d", "isolated"], np.array([0, 2], dtype=np.int32),
                           np.array([1, 3], dtype=np.int32), np.array([2., 2.]))
    roots = [tmp_path / "original", tmp_path / "readback"]
    for root in roots:
        (root / "orthohmm_working_res").mkdir(parents=True)
    write_indexed_edges(original, str(roots[0]))
    readback = runner.load_native_edges(roots[0] / "orthohmm_working_res/orthohmm_edges.txt", original.gene_names, 2)
    for root, edges in zip(roots, (original, readback)):
        execute_leiden(.1, str(root), edges=edges, include_isolates=True, seed=4)
    partitions = [read_index_clusters(str(root / "orthohmm_working_res/orthohmm_edges_clustered.txt"), original.gene_names)
                  for root in roots]
    assert partitions[0] == partitions[1]
    assert [4] in partitions[0]


@pytest.mark.parametrize("cpu", [0, 17, True])
def test_preparation_checks_budget_before_inputs(tmp_path, cpu):
    with pytest.raises(ValueError, match="1-16 CPUs"):
        runner.prepare(tmp_path, tmp_path, tmp_path / "output", tmp_path / "comparison", cpu)


def setup_runner(tmp_path, monkeypatch):
    from benchmark_tools.prepare_ob_candidate_neighborhood import record
    source = Path(runner.__file__).resolve().parent.parent
    runtime = {"files": []}
    plan = {"source_root": str(source), "repo": str(source), "output": str(tmp_path),
            "source_revision": "test", "scope": "test", "inputs": [], "sources": [],
            "runtime": runtime, "tool_path": "/usr/bin", "cases": []}
    for dataset in runner.DATASETS:
        plan["cases"].append({"dataset": dataset, "arms": [
            {"label": label, "output": str(tmp_path / dataset / label),
             "command": ["--species-tree-mode", "infer", dataset, label]} for label in runner.LABELS]})
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(plan))
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    monkeypatch.setattr(runner, "runtime", lambda: runtime)
    monkeypatch.setattr(runner, "host_activity", lambda: {})
    def construct(case, report, folder):
        for arm in case["arms"]:
            directory = Path(arm["output"])
            directory.mkdir()
            candidate = {}
            for name in ("candidate_partition", "membership_constraints", "seed_partition"):
                file = directory / name
                file.write_text("fixture\n")
                candidate[name] = record(file)
            report["arms"].append({"label": arm["label"], "status": "candidates_verified_unscored", "candidate": candidate})
    monkeypatch.setattr(runner, "construct", construct)
    return path


def test_failed_control_gates_variant_without_hiding_other_dataset(tmp_path, monkeypatch):
    path = setup_runner(tmp_path, monkeypatch)
    calls = []
    def failed(command, *_args):
        calls.append(command[-2:])
        return {"returncode": 1}
    monkeypatch.setattr(runner, "timed_stage", failed)
    result = runner.run(path)
    assert calls == [[dataset, "control"] for dataset in runner.DATASETS]
    assert result["status"] == "finished_with_explicit_gaps"
    assert len(result["datasets"]) == 2
    for report in result["datasets"]:
        assert report["arms"][0]["status"] == "failed"
        assert report["arms"][1]["status"] == "not_run_control_gate_failed"
    with pytest.raises(ValueError, match="existing attempt"):
        runner.run(path)


def test_supplied_species_tree_cannot_enter_inferred_followup(tmp_path, monkeypatch):
    path = setup_runner(tmp_path, monkeypatch)
    plan = json.loads(path.read_text())
    plan["cases"][0]["arms"][1]["command"].extend(["--species-tree", "baseline.nwk"])
    path.write_text(json.dumps(plan))
    with pytest.raises(ValueError, match="candidate-specific inferred trees"):
        runner.run(path)
    assert not (tmp_path / "results.json").exists()


def test_preparation_failure_retained_and_next_dataset_attempted(tmp_path, monkeypatch):
    path = setup_runner(tmp_path, monkeypatch)
    attempted = []
    def failed(case, *_args):
        attempted.append(case["dataset"])
        raise ValueError("Default candidate control mismatch")
    monkeypatch.setattr(runner, "construct", failed)
    result = runner.run(path)
    assert attempted == list(runner.DATASETS)
    assert all(r["status"] == "failed_preparation_or_gate" and r["arms"] == [] for r in result["datasets"])
    assert result["defaults_promoted"] is False and result["goal_completed"] is False


def test_only_failed_orthobench_dataset_is_attempted_in_retry(tmp_path, monkeypatch):
    path = setup_runner(tmp_path, monkeypatch)
    plan = json.loads(path.read_text())
    plan["datasets"] = ["OrthoBench"]
    plan["cases"] = plan["cases"][:1]
    path.write_text(json.dumps(plan))
    calls = []
    def failed(command, *_args):
        calls.append(command[-2:])
        return {"returncode": 1}
    monkeypatch.setattr(runner, "timed_stage", failed)
    result = runner.run(path)
    assert calls == [["OrthoBench", "control"]]
    assert len(result["datasets"]) == 1
    assert not (tmp_path / "ThreeKingdoms").exists()
