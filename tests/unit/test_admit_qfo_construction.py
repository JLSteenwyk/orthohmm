from copy import deepcopy

import igraph
import numpy as np
import pytest

from benchmark_tools.admit_qfo_construction import check_observation
from benchmark_tools.probe_leiden_boundary import graph_fingerprint


def observation(mismatch, mode):
    sources, targets = np.array([0, 1, 2], dtype=np.int32), np.array([1, 2, 3], dtype=np.int32)
    graph = igraph.Graph(n=4, edges=list(zip(sources.tolist(), targets.tolist())), directed=False)
    graph.es["weight"] = [1., 2., 3.]
    saved = graph_fingerprint(graph)
    native = igraph.Graph(n=4, edges=[(0, 1), (0, 2) if mismatch else (1, 2), (2, 3)], directed=False)
    native.es["weight"] = [1., 2., 3.]
    examples = [{"edge_index": 1, "left": [0, 2], "right": [1, 2]}] if mismatch else []
    diff = {"different_edges": len(examples), "examples": examples}
    original = {"constructor_dtype": "int32", "constructor_c_contiguous": True,
                "constructor_vs_saved": {"different_edges": 0, "examples": []},
                "native_vs_constructor": diff, "native_vs_saved": diff}
    dtype = "int64" if mode == "explicit_int64" else "int32"
    result = {"status": "construction_observed_without_optimization", "mode": mode,
              "accuracy_evaluated": False, "optimizer_called": False, "saved": saved,
              "native": graph_fingerprint(native), "differences": original,
              "converted_differences": {**deepcopy(original), "constructor_dtype": dtype},
              "converted_dtype": dtype, "witnesses": [{"edge_index": 1, "tuple": [0, 2],
                  "source_target": [0, 2], "expected_pair_edge_id": -1}] if mismatch else []}
    return {"mode": mode, "index": 0, "result": result}, saved, sources, targets


@pytest.mark.parametrize("mode", ["original_int32", "explicit_int64"])
@pytest.mark.parametrize("mismatch", [False, True])
def test_full_hash_reconstruction(mode, mismatch):
    args = observation(mismatch, mode)
    before = args[2].copy(), args[3].copy()
    result = check_observation(*args)
    assert result["mismatched_edges"] == int(mismatch)
    assert result["complete_endpoint_hash_reconstructed"] is True
    assert np.array_equal(before[0], args[2]) and np.array_equal(before[1], args[3])


@pytest.mark.parametrize("problem", ["optimized", "dtype", "input", "hash", "weights", "missing",
                                     "index", "expected", "native_access", "edge_lookup"])
def test_corrupt_evidence_rejected(problem):
    args = observation(True, "explicit_int64")
    result = args[0]["result"]
    if problem == "optimized":
        result["optimizer_called"] = True
    elif problem == "dtype":
        result["converted_dtype"] = "int32"
    elif problem == "input":
        result["converted_differences"]["constructor_vs_saved"]["different_edges"] = 1
    elif problem in {"hash", "weights"}:
        result["native"]["ordered_endpoints_sha256" if problem == "hash" else "ordered_weights_sha256"] = "bad"
    elif problem == "missing":
        result["witnesses"] = []
    elif problem in {"index", "expected"}:
        for key in ("differences", "converted_differences"):
            example = result[key]["native_vs_saved"]["examples"][0]
            example["edge_index" if problem == "index" else "right"] = 99 if problem == "index" else [1, 3]
    elif problem == "native_access":
        result["witnesses"][0]["source_target"] = [1, 2]
    elif problem == "edge_lookup":
        result["witnesses"][0]["expected_pair_edge_id"] = 1
    with pytest.raises(ValueError):
        check_observation(*args)
