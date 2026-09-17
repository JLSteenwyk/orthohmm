from copy import deepcopy
import hashlib

import numpy as np
import pytest

from benchmark_tools.admit_qfo_direct_graph import check_panel, check_stages


def fixture(mismatch):
    sources, targets = np.array([0, 1, 2], dtype=np.int32), np.array([1, 2, 3], dtype=np.int32)
    expected = np.column_stack((sources, targets)).astype("<i8")
    saved = {"vertices": 4, "edges": 3, "directed": False,
             "ordered_endpoints_sha256": hashlib.sha256(expected.tobytes()).hexdigest(),
             "ordered_weights_sha256": "weights"}
    observed = expected.copy()
    if mismatch:
        observed[1] = [0, 2]
    examples = [{"edge_index": 1, "left": [0, 2], "right": [1, 2]}] if mismatch else []
    delta = {"different_edges": len(examples), "examples": examples}
    differences = {"constructor_dtype": "int32", "constructor_c_contiguous": True,
                   "constructor_vs_saved": {"different_edges": 0, "examples": []},
                   "native_vs_saved": delta, "native_vs_constructor": deepcopy(delta)}
    result = {"status": "direct_construction_observed", "mode": "minimal_imports", "accuracy_evaluated": False,
              "optimizer_called": False, "saved": saved,
              "before_weights": {"vertices": 4, "edges": 3, "directed": False, "differences": differences},
              "after_weights": {"fingerprint": {**saved, "ordered_endpoints_sha256": hashlib.sha256(observed.tobytes()).hexdigest()},
                                "differences": deepcopy(differences)}}
    return result, "minimal_imports", saved, sources, targets


@pytest.mark.parametrize("mismatch", [False, True])
def test_both_stages_reconstruct(mismatch):
    args = fixture(mismatch)
    sources, targets = args[3].copy(), args[4].copy()
    report = check_stages(*args)
    assert report["before_mismatched_edges"] == report["after_mismatched_edges"] == int(mismatch)
    assert report["stage_witnesses_equal"] is True
    assert report["before_witness_implied_hash"] == report["after_native_hash_reconstructed"]
    assert np.array_equal(args[3], sources) and np.array_equal(args[4], targets)


@pytest.mark.parametrize("problem", ["mode", "optimized", "shape", "input", "count", "index", "endpoint", "hash", "weights"])
def test_inconsistent_stages_rejected(problem):
    args = fixture(True)
    result = args[0]
    if problem == "mode":
        result["mode"] = "frozen_imports"
    elif problem == "optimized":
        result["optimizer_called"] = True
    elif problem == "shape":
        result["before_weights"]["edges"] = 4
    elif problem == "input":
        result["before_weights"]["differences"]["constructor_vs_saved"]["different_edges"] = 1
    elif problem == "count":
        result["before_weights"]["differences"]["native_vs_saved"]["different_edges"] = 2
    elif problem in {"index", "endpoint"}:
        for key in ("native_vs_saved", "native_vs_constructor"):
            row = result["before_weights"]["differences"][key]["examples"][0]
            row["edge_index" if problem == "index" else "right"] = 99 if problem == "index" else [1, 3]
    else:
        result["after_weights"]["fingerprint"]["ordered_endpoints_sha256" if problem == "hash" else "ordered_weights_sha256"] = "bad"
    with pytest.raises(ValueError):
        check_stages(*args)


@pytest.mark.parametrize("problem", [None, "missing", "duplicate", "running", "exit", "job", "scored"])
def test_exact_terminal_panel(problem):
    report = {"status": "six_direct_workers_complete_unscored", "job_id": "21326",
              "accuracy_evaluated": False, "optimizer_called": False,
              "workers": [{"index": i, "mode": mode, "exit_code": 0} for i in range(3)
                          for mode in ("minimal_imports", "frozen_imports")]}
    if problem == "missing":
        report["workers"].pop()
    elif problem == "duplicate":
        report["workers"][-1] = report["workers"][0]
    elif problem == "running":
        report["status"] = "running"
    elif problem == "exit":
        report["workers"][0]["exit_code"] = 1
    elif problem == "job":
        report["job_id"] = "21327"
    elif problem == "scored":
        report["accuracy_evaluated"] = True
    if problem:
        with pytest.raises(ValueError):
            check_panel(report)
    else:
        check_panel(report)
