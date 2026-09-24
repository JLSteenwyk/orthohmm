import json
from types import SimpleNamespace

import numpy as np
import pytest

from benchmark_tools.stop_before_leiden import DiagnosticStop, stop_before_partition


def setup(tmp_path):
    igraph = pytest.importorskip("igraph")
    (tmp_path / "gene_names.txt").write_text("a\nb\nc\n")
    for name, values in (("sources", np.array([0, 1], dtype=np.int32)),
                         ("targets", np.array([1, 2], dtype=np.int32)),
                         ("weights", np.array([.5, .8], dtype=np.float64))):
        np.save(tmp_path / (name + ".npy"), values)
    graph = igraph.Graph(n=3, edges=[(0, 1), (1, 2)], directed=False)
    graph.es["weight"] = [.5, .8]
    def forbidden(*args, **kwargs):
        raise AssertionError("Optimizer must never execute")
    return graph, SimpleNamespace(find_partition=forbidden)


def invoke(module, graph, edges=None):
    graph_edges = np.array([[0, 1], [1, 2]], dtype=np.int32) if edges is None else edges
    module.find_partition(graph)


def test_native_small_graph_stops_before_optimizer(tmp_path):
    graph, module = setup(tmp_path)
    original = module.find_partition
    with stop_before_partition(module, tmp_path) as calls:
        with pytest.raises(DiagnosticStop):
            invoke(module, graph)
    assert module.find_partition is original and calls == [True]
    result = json.loads((tmp_path / "preoptimizer_stop.json").read_text())
    assert result["optimizer_called"] is False
    assert result["fingerprint"] == result["saved"]
    assert result["fingerprint"]["vertices"] == 3
    with pytest.raises(FileExistsError):
        with stop_before_partition(module, tmp_path):
            pass


@pytest.mark.parametrize("change", ["weights", "endpoints", "constructor", "missing_constructor"])
def test_boundary_mismatch_does_not_write_success(tmp_path, change):
    graph, module = setup(tmp_path)
    original = module.find_partition
    if change == "weights":
        graph.es["weight"] = [.4, .8]
    elif change == "endpoints":
        graph.delete_edges(0)
        graph.add_edge(0, 2, weight=.5)
    with pytest.raises(ValueError):
        with stop_before_partition(module, tmp_path):
            if change == "missing_constructor":
                module.find_partition(graph)
            else:
                edges = np.array([[0, 2], [1, 2]], dtype=np.int32) if change == "constructor" else None
                invoke(module, graph, edges)
    assert module.find_partition is original
    assert not (tmp_path / "preoptimizer_stop.json").exists()
