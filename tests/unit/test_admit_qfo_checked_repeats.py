import json

import igraph
import leidenalg
import numpy as np
import pytest

from benchmark_tools.admit_qfo_checked_repeats import check_panel, check_native, reconstruct, admit
from benchmark_tools.checked_python_pair_worker import python_pair_constructor
from benchmark_tools.probe_leiden_boundary import observe_partition, saved_fingerprint
from benchmark_tools.run_qfo_checked_repeats import constructor_digest


@pytest.mark.parametrize("problem", [None, "running", "job", "accuracy", "planned", "missing", "order", "exit"])
def test_panel(problem):
    report = {"status": "three_checked_repeats_complete_unscored", "accuracy_evaluated": False,
              "job_id": "21328", "planned_repeats": 3,
              "repeats": [{"index": i, "execution": {"exit_code": 0}} for i in range(3)]}
    if problem == "running":
        report["status"] = "running"
    elif problem == "job":
        report["job_id"] = "21327"
    elif problem == "accuracy":
        report["accuracy_evaluated"] = True
    elif problem == "planned":
        report["planned_repeats"] = 4
    elif problem == "missing":
        report["repeats"].pop()
    elif problem == "order":
        report["repeats"].reverse()
    elif problem == "exit":
        report["repeats"][1]["execution"]["exit_code"] = 1
    if problem:
        with pytest.raises(ValueError):
            check_panel(report)
    else:
        check_panel(report)


def save_graph(tmp_path, pairs):
    np.save(tmp_path / "sources.npy", pairs[:, 0])
    np.save(tmp_path / "targets.npy", pairs[:, 1])
    np.save(tmp_path / "weights.npy", np.ones(len(pairs), dtype=np.float64))
    (tmp_path / "gene_names.txt").write_text("a\nb\nc\nd\n")


def test_independent_reconstruction_chunked(tmp_path):
    pairs = np.tile(np.array([[2, 0], [1, 1], [0, 2]], dtype=np.int32), (25000, 1))
    save_graph(tmp_path, pairs)
    saved, oriented, universe = reconstruct(tmp_path)
    assert saved == saved_fingerprint(tmp_path)
    assert oriented == constructor_digest(tmp_path)
    assert universe == {"a", "b", "c", "d"}


@pytest.mark.parametrize("problem", [None, "before", "after", "seed", "input", "duplicate", "failed"])
def test_native_observation(tmp_path, problem):
    pairs = np.array([[2, 0], [1, 1], [0, 2]], dtype=np.int32)
    save_graph(tmp_path, pairs)
    with python_pair_constructor(igraph, tmp_path / "constructor_adapter.json"):
        graph = igraph.Graph(n=4, edges=pairs, directed=False)
        graph.es["weight"] = [1., 1., 1.]
        with observe_partition(leidenalg, tmp_path):
            leidenalg.find_partition(graph, leidenalg.CPMVertexPartition, weights="weight", resolution_parameter=.1, seed=4)
    boundary = json.loads((tmp_path / "native_boundary.json").read_text())
    adapter = json.loads((tmp_path / "constructor_adapter.json").read_text())
    saved, oriented, _ = reconstruct(tmp_path)
    if problem in ("before", "after"):
        boundary["calls"][0][problem]["ordered_endpoints_sha256"] = "changed"
    elif problem == "seed":
        boundary["calls"][0]["arguments"]["seed"] = 7
    elif problem == "input":
        adapter["calls"][0]["ordered_input_bytes_sha256"] = "changed"
    elif problem == "duplicate":
        adapter["calls"].append(adapter["calls"][0])
    elif problem == "failed":
        boundary["calls"][0]["status"] = "failed"
    if problem:
        with pytest.raises(ValueError):
            check_native(boundary, adapter, saved, oriented)
    else:
        check_native(boundary, adapter, saved, oriented)


@pytest.mark.parametrize("problem", ["negative", "outside", "nonfinite", "names", "dtype"])
def test_invalid_saved_graph(tmp_path, problem):
    pairs = np.array([[2, 0], [1, 1]], dtype=np.int32)
    if problem == "negative":
        pairs[0, 0] = -1
    elif problem == "outside":
        pairs[0, 0] = 4
    elif problem == "dtype":
        pairs = pairs.astype(np.int64)
    save_graph(tmp_path, pairs)
    if problem == "nonfinite":
        np.save(tmp_path / "weights.npy", np.array([np.nan, 1.]))
    elif problem == "names":
        (tmp_path / "gene_names.txt").write_text("a\na\nc\nd\n")
    with pytest.raises(ValueError):
        reconstruct(tmp_path)


def test_refuses_existing_admission(tmp_path):
    with pytest.raises(FileExistsError):
        admit(tmp_path, tmp_path, "unused")
