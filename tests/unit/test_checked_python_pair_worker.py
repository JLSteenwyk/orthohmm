import json

import igraph
import leidenalg
import numpy as np
import pytest

from benchmark_tools.checked_python_pair_worker import python_pair_constructor, require_admission
from benchmark_tools.probe_leiden_boundary import observe_partition


@pytest.mark.parametrize("problem", [None, "pending", "missing", "order", "format", "mode", "before", "after", "optimized"])
def test_admission_gate(problem):
    report = {"status": "constructor_format_observations_verified", "accuracy_evaluated": False,
              "optimizer_called": False, "workers": [{"index": i, "edge_format": fmt,
                "mode": "minimal_imports", "before_mismatched_edges": 0, "after_mismatched_edges": 0}
                for i in range(3) for fmt in ("numpy", "python_pairs")]}
    if problem == "pending":
        report["status"] = "running"
    elif problem == "missing":
        report["workers"].pop()
    elif problem == "order":
        report["workers"].reverse()
    elif problem == "format":
        report["workers"][1]["edge_format"] = "numpy"
    elif problem == "mode":
        report["workers"][1]["mode"] = "frozen_imports"
    elif problem in {"before", "after"}:
        report["workers"][1][problem + "_mismatched_edges"] = 6
    elif problem == "optimized":
        report["optimizer_called"] = True
    if problem:
        with pytest.raises(ValueError):
            require_admission(report)
    else:
        require_admission(report)


@pytest.mark.parametrize("corrupt", [False, True])
def test_checked_constructor_and_optimizer(tmp_path, corrupt):
    pairs = np.array([[2, 0], [1, 1], [0, 2]], dtype=np.int32)
    np.save(tmp_path / "sources.npy", pairs[:, 0])
    np.save(tmp_path / "targets.npy", pairs[:, 1])
    np.save(tmp_path / "weights.npy", np.array([1., 2., 3.]))
    (tmp_path / "gene_names.txt").write_text("a\nb\nc\nd\n")
    original = igraph.Graph.__init__
    if corrupt:
        pairs[0] = [3, 0]
    with python_pair_constructor(igraph, tmp_path / "constructor_adapter.json") as calls:
        graph = igraph.Graph(n=4, edges=pairs, directed=False)
        graph.es["weight"] = [1., 2., 3.]
        with observe_partition(leidenalg, tmp_path):
            if corrupt:
                with pytest.raises(ValueError, match="differs"):
                    leidenalg.find_partition(graph, leidenalg.CPMVertexPartition, weights="weight",
                                            resolution_parameter=.1, seed=4)
            else:
                partition = leidenalg.find_partition(graph, leidenalg.CPMVertexPartition, weights="weight",
                                                    resolution_parameter=.1, seed=4)
                assert sum(map(len, partition)) == 4
        assert len(calls) == 1 and calls[0]["status"] == "constructor_returned"
    assert igraph.Graph.__init__ is original
    boundary = json.loads((tmp_path / "native_boundary.json").read_text())["calls"]
    assert len(boundary) == 1
    if corrupt:
        assert boundary[0]["status"] == "native_graph_mismatch_before_optimizer"
        assert "after" not in boundary[0]
    else:
        assert boundary[0]["before"] == boundary[0]["after"] == boundary[0]["saved"]
        assert boundary[0]["arguments"]["seed"] == 4


def test_constructor_restores_on_failure_and_rejects_second_array(tmp_path):
    original = igraph.Graph.__init__
    pairs = np.array([[0, 1]], dtype=np.int32)
    with pytest.raises(ValueError, match="one keyword"):
        with python_pair_constructor(igraph, tmp_path / "adapter.json"):
            igraph.Graph(n=2, edges=pairs, directed=False)
            igraph.Graph(n=2, edges=pairs, directed=False)
    assert igraph.Graph.__init__ is original
