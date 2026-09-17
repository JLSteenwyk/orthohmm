import json
from types import SimpleNamespace

import igraph
import leidenalg
import numpy as np
import pytest

from benchmark_tools.probe_leiden_boundary import graph_fingerprint, saved_fingerprint, observe_partition, endpoint_differences


@pytest.fixture
def native(tmp_path):
    pairs = [(2, 0), (1, 1), (0, 2)]
    weights = [1., 2., 3.]
    graph = igraph.Graph(n=4, edges=pairs, directed=False)
    graph.es["weight"] = weights
    np.save(tmp_path / "sources.npy", [pair[0] for pair in pairs])
    np.save(tmp_path / "targets.npy", [pair[1] for pair in pairs])
    np.save(tmp_path / "weights.npy", weights)
    (tmp_path / "gene_names.txt").write_text("a\nb\nc\nd\n")
    return graph, tmp_path


def test_chunked_native_equals_saved(native):
    graph, payload = native
    assert graph_fingerprint(graph, 1) == saved_fingerprint(payload, 2)
    assert graph_fingerprint(graph, 20) == graph_fingerprint(graph, 1)


def test_weight_and_order_sensitive(native):
    graph, _ = native
    before = graph_fingerprint(graph)
    graph.es["weight"] = [3., 2., 1.]
    assert graph_fingerprint(graph) != before
    other = igraph.Graph(n=4, edges=[(1, 1), (0, 2), (0, 2)])
    other.es["weight"] = [2., 1., 3.]
    assert graph_fingerprint(other)["ordered_endpoints_sha256"] != before["ordered_endpoints_sha256"]


def test_actual_optimizer_defaults_and_restoration(native):
    graph, payload = native
    original = leidenalg.find_partition
    with observe_partition(leidenalg, payload):
        partition = leidenalg.find_partition(graph, leidenalg.CPMVertexPartition,
                                            weights="weight", resolution_parameter=.1, seed=4)
    assert leidenalg.find_partition is original
    assert sum(map(len, partition)) == 4
    calls = json.loads((payload / "native_boundary.json").read_text())["calls"]
    assert len(calls) == 1
    row = calls[0]
    assert row["before"] == row["after"] == row["saved"]
    assert row["arguments"] == {"initial_membership": None, "weights": "weight",
        "n_iterations": 2, "max_comm_size": 0, "seed": 4,
        "kwargs": {"resolution_parameter": .1},
        "partition_type": "leidenalg.VertexPartition.CPMVertexPartition"}


def test_failure_restores_and_records(native):
    graph, payload = native
    def failing(graph, partition_type):
        raise RuntimeError("native failure")
    module = SimpleNamespace(find_partition=failing)
    with pytest.raises(RuntimeError, match="native failure"):
        with observe_partition(module, payload):
            module.find_partition(graph, leidenalg.CPMVertexPartition)
    assert module.find_partition is failing
    assert json.loads((payload / "native_boundary.json").read_text())["calls"][0]["status"] == "optimizer_failed"


def test_mismatch_rejected_before_optimizer(native):
    graph, payload = native
    graph.es[0]["weight"] = 10.
    with pytest.raises(ValueError, match="differs"):
        with observe_partition(leidenalg, payload):
            leidenalg.find_partition(graph, leidenalg.CPMVertexPartition, weights="weight")


def test_nonfinite_and_invalid_chunks(native):
    graph, payload = native
    with pytest.raises(ValueError, match="positive"):
        saved_fingerprint(payload, 0)
    with pytest.raises(ValueError, match="positive"):
        graph_fingerprint(graph, 0)
    graph.es[0]["weight"] = float("nan")
    with pytest.raises(ValueError, match="Nonfinite"):
        graph_fingerprint(graph)


def test_constructor_and_native_differences_are_separated(native):
    graph, payload = native
    constructor = np.array([(2, 0), (1, 1), (0, 2)], dtype=np.int32)
    original = constructor.copy()
    matching = endpoint_differences(graph, payload, constructor, 1)
    assert all(matching[key]["different_edges"] == 0 for key in
               ("native_vs_saved", "constructor_vs_saved", "native_vs_constructor"))
    constructor[0] = [3, 0]
    result = endpoint_differences(graph, payload, constructor, 2)
    assert result["constructor_vs_saved"]["different_edges"] == 1
    assert result["native_vs_constructor"]["different_edges"] == 1
    assert result["native_vs_saved"]["different_edges"] == 0
    assert result["constructor_vs_saved"]["examples"] == [{"edge_index": 0, "left": [0, 3], "right": [0, 2]}]
    assert np.array_equal(constructor[1:], original[1:])


def test_mismatch_captures_callers_constructor_input(native):
    _, payload = native
    graph_edges = np.array([(3, 0), (1, 1), (0, 2)], dtype=np.int32)
    graph = igraph.Graph(n=4, edges=graph_edges)
    graph.es["weight"] = [1., 2., 3.]
    with pytest.raises(ValueError, match="differs"):
        with observe_partition(leidenalg, payload):
            leidenalg.find_partition(graph, leidenalg.CPMVertexPartition, weights="weight")
    row = json.loads((payload / "native_boundary.json").read_text())["calls"][0]
    assert row["status"] == "native_graph_mismatch_before_optimizer"
    assert row["endpoint_differences"]["native_vs_saved"]["different_edges"] == 1
    assert row["endpoint_differences"]["constructor_vs_saved"]["different_edges"] == 1
    assert row["endpoint_differences"]["native_vs_constructor"]["different_edges"] == 0
