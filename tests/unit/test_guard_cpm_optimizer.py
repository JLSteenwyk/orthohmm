import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import igraph
import leidenalg
import numpy as np
import pytest

from benchmark_tools.guard_cpm_optimizer import require_exact_constructor, run_frozen_worker
from benchmark_tools.probe_leiden_boundary import observe_partition


def payload(path):
    path.mkdir()
    (path / "gene_names.txt").write_text("a\nb\nc\nd\n")
    for name, data in (("sources", np.array([0, 2], dtype=np.int32)),
                       ("targets", np.array([1, 3], dtype=np.int32)),
                       ("weights", np.array([1., 2.]))):
        np.save(path / f"{name}.npy", data)
    return path


@pytest.mark.parametrize("problem", [None, "constructor", "native", "weight", "missing", "repeated", "existing", "raises"])
def test_constructor_guard(tmp_path, problem):
    path = payload(tmp_path / "payload")
    original_calls = []
    def original(graph, *args, **kwargs):
        original_calls.append((args, kwargs))
        if problem == "raises":
            raise RuntimeError("native failure")
        return "result"
    module = SimpleNamespace(find_partition=original)
    graph = igraph.Graph(n=4, edges=[(0, 1), (2, 3)])
    graph.es["weight"] = [1., 2.]
    if problem == "native":
        graph.delete_edges(0)
        graph.add_edge(0, 2, weight=1.)
    elif problem == "weight":
        graph.es[0]["weight"] = 9.
    elif problem == "existing":
        (path / "constructor_parity.json").write_text("preserve")
    def invoke():
        graph_edges = np.array([[0, 1], [2, 3]], dtype=np.int32)
        if problem == "constructor":
            graph_edges[0, 1] = 2
        elif problem == "missing":
            del graph_edges
        return module.find_partition(graph, "partition", seed=4, resolution_parameter=.12)
    if problem:
        with pytest.raises((ValueError, FileExistsError, RuntimeError)):
            with require_exact_constructor(module, path):
                invoke()
                if problem == "repeated":
                    invoke()
        assert len(original_calls) == (1 if problem in ("repeated", "raises") else 0)
        if problem == "existing":
            assert (path / "constructor_parity.json").read_text() == "preserve"
    else:
        with require_exact_constructor(module, path):
            assert invoke() == "result"
        assert original_calls == [(("partition",), {"seed": 4, "resolution_parameter": .12})]
        report = json.loads((path / "constructor_parity.json").read_text())
        assert report["fingerprint"] == report["saved"]
        assert all(report["differences"][name]["different_edges"] == 0 for name in (
            "native_vs_saved", "constructor_vs_saved", "native_vs_constructor"))
    assert module.find_partition is original


def test_real_optimizer_observer_composition(tmp_path):
    path = payload(tmp_path / "payload")
    graph_edges = np.array([[0, 1], [2, 3]], dtype=np.int32)
    graph = igraph.Graph(n=4, edges=[tuple(map(int, edge)) for edge in graph_edges])
    graph.es["weight"] = [1., 2.]
    original = leidenalg.find_partition
    with observe_partition(leidenalg, path):
        with require_exact_constructor(leidenalg, path):
            result = leidenalg.find_partition(graph, leidenalg.CPMVertexPartition,
                weights="weight", seed=4, resolution_parameter=.12)
    assert len(result.membership) == 4
    assert leidenalg.find_partition is original
    call = json.loads((path / "native_boundary.json").read_text())["calls"][0]
    assert call["status"] == "optimizer_returned"
    assert call["before"] == call["saved"] == call["after"]
    assert call["arguments"]["kwargs"] == {"resolution_parameter": .12}


def test_real_frozen_worker_subprocess(tmp_path):
    root = Path(__file__).resolve().parents[2]
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    if not launcher.exists():
        pytest.skip("Requires retained frozen publication worker")
    path = payload(tmp_path / "payload")
    (path / "metadata.json").write_text(json.dumps(dict(cpm_resolution=.12, seed=4,
        include_isolates=True, output_directory=str(tmp_path))))
    code = '''
import os, sys
from pathlib import Path
root, launcher, payload = map(Path, sys.argv[1:])
sys.path.insert(0, str(root))
sys.path.insert(0, str(root/'benchmark_tools'))
from benchmark_tools.guard_cpm_optimizer import run_frozen_worker
run_frozen_worker(launcher, payload)
'''
    done = subprocess.run([sys.executable, "-I", "-B", "-c", code, str(root), str(launcher), str(path)],
                          cwd=launcher, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert json.loads((path / "native_boundary.json").read_text())["calls"][0]["status"] == "optimizer_returned"
    assert json.loads((path / "constructor_parity.json").read_text())["status"] == "constructor_parity_verified_before_optimizer"
    groups = [set(line.split()) for line in (tmp_path / "orthohmm_working_res/orthohmm_edges_clustered.txt").read_text().splitlines()]
    assert set.union(*groups) == {"a", "b", "c", "d"}


@pytest.mark.parametrize("problem", ["cwd", "resolution", "output", "constructor_parity.json",
    "constructor_adapter.json", "native_boundary.json", "worker_before.json", "working"])
def test_worker_rejects_stale_or_changed_setup(tmp_path, monkeypatch, problem):
    path = payload(tmp_path / "payload")
    metadata = dict(cpm_resolution=.12, seed=4, include_isolates=True, output_directory=str(tmp_path))
    if problem == "resolution":
        metadata["cpm_resolution"] = .1
    elif problem == "output":
        metadata["output_directory"] = str(tmp_path / "elsewhere")
    elif problem == "working":
        (tmp_path / "orthohmm_working_res").mkdir()
    elif problem.endswith(".json"):
        (path / problem).write_text("preserve")
    (path / "metadata.json").write_text(json.dumps(metadata))
    monkeypatch.chdir(path if problem == "cwd" else tmp_path)
    with pytest.raises((ValueError, FileExistsError)):
        run_frozen_worker(tmp_path, path)
    if problem.endswith(".json"):
        assert (path / problem).read_text() == "preserve"
