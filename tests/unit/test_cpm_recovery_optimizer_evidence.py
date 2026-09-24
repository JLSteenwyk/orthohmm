import copy
import json
from pathlib import Path
import sys

import pytest

from benchmark_tools import run_cpm_checkpoint_recovery as module
from benchmark_tools import probe_leiden_boundary, repeat_qfo_saved_graph
from benchmark_tools.prepare_ob_candidate_neighborhood import record


@pytest.fixture(scope="module")
def universe(tmp_path_factory):
    root = tmp_path_factory.mktemp("optimizer_evidence")
    output = root / "recovery"
    payload = output / "payload"
    payload.mkdir(parents=True)
    names = [f"g{i}" for i in range(984137)]
    (payload / "gene_names.txt").write_text("\n".join(names) + "\n")
    partition = output / "orthohmm_working_res/orthohmm_edges_clustered.txt"
    partition.parent.mkdir()
    partition.write_text(" ".join(names) + "\n")
    for name in ("sources.npy", "targets.npy", "weights.npy"):
        (payload / name).write_text("fixture")
    return root, output, payload, partition, names


@pytest.mark.parametrize("problem", [None, "snapshot_gate", "affinity", "inputs", "observer", "python",
    "fingerprint", "unfinished", "settings", "after", "extra_call", "adapter", "parity",
    "parity_weights", "duplicate", "missing", "unknown", "missing_artifact", "mutated_record"])
def test_optimizer_evidence(universe, monkeypatch, problem):
    # This is an evidence-validator test, not native graph execution. Native
    # construction/optimization has a separate real frozen-worker subprocess test.
    root, output, payload, partition, names = universe
    partition.write_text(" ".join(names) + "\n")
    saved = dict(vertices=984137, edges=25501180, directed=False,
                 ordered_endpoints_sha256="endpoints", ordered_weights_sha256="weights")
    environment = dict(PYTHONPATH=str(root), PYTHONHASHSEED="0", OMP_NUM_THREADS="1",
                       OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    snapshot = dict(cpu_affinity=[0], inputs=[record(payload / name) for name in (
        "gene_names.txt", "sources.npy", "targets.npy", "weights.npy")],
        observer=record(Path(module.__file__).with_name("repeat_qfo_saved_graph.py")),
        python=record(sys.executable), modules={}, native_libraries=[])
    arguments = dict(initial_membership=None, weights="weight", n_iterations=2, max_comm_size=0,
        seed=4, kwargs={"resolution_parameter": .12}, partition_type="leidenalg.VertexPartition.CPMVertexPartition")
    boundary = dict(accuracy_evaluated=False, calls=[dict(arguments=arguments, before=saved, saved=saved,
                                                       after=saved, status="optimizer_returned")])
    adapter = dict(format="python_pairs", accuracy_evaluated=False, calls=[dict(status="constructor_returned",
        shape=[25501180, 2], dtype="int32", n=984137, directed=False,
        ordered_input_bytes_sha256="93b36aad4916b85195c4dcc7dccb37d944cb911f7b051cf989ad1479729fded4")])
    parity = dict(status="constructor_parity_verified_before_optimizer", accuracy_evaluated=False,
        fingerprint=saved, saved=saved, differences={key: {"different_edges": 0} for key in (
            "native_vs_saved", "constructor_vs_saved", "native_vs_constructor")})
    snapshot, boundary, adapter, parity = copy.deepcopy((snapshot, boundary, adapter, parity))
    def snapshot_gate(observed, launcher, directory, overrides, **kwargs):
        assert launcher == root and directory == payload and overrides == environment
        assert kwargs == {"expected_cpm_resolution": .12}
        if problem == "snapshot_gate":
            raise ValueError("worker snapshot rejected")
    monkeypatch.setattr(repeat_qfo_saved_graph, "check_worker", snapshot_gate)
    def fingerprint(path):
        if problem == "mutated_record":
            (payload / "worker_before.json").write_text("changed after reading")
        return {**saved, "edges": 1} if problem == "fingerprint" else saved
    monkeypatch.setattr(probe_leiden_boundary, "saved_fingerprint", fingerprint)
    if problem in ("affinity", "inputs", "observer", "python"):
        snapshot[{"affinity": "cpu_affinity"}.get(problem, problem)] = [] if problem in ("affinity", "inputs") else {}
    elif problem == "unfinished":
        boundary["calls"][0]["status"] = "before_optimizer"
    elif problem == "settings":
        boundary["calls"][0]["arguments"]["kwargs"] = {"resolution_parameter": .1}
    elif problem == "after":
        boundary["calls"][0]["after"] = {}
    elif problem == "extra_call":
        boundary["calls"].append(copy.deepcopy(boundary["calls"][0]))
    elif problem == "adapter":
        adapter["calls"][0]["ordered_input_bytes_sha256"] = "changed"
    elif problem == "parity":
        parity["differences"]["constructor_vs_saved"]["different_edges"] = 1
    elif problem == "parity_weights":
        parity["fingerprint"] = {**saved, "ordered_weights_sha256": "changed"}
    elif problem == "duplicate":
        partition.write_text(" ".join(names) + " g0\n")
    elif problem == "missing":
        partition.write_text(" ".join(names[:-1]) + "\n")
    elif problem == "unknown":
        partition.write_text(" ".join(names) + " foreign\n")
    for name, value in zip(("worker_before.json", "native_boundary.json", "constructor_adapter.json", "constructor_parity.json"),
                           (snapshot, boundary, adapter, parity)):
        (payload / name).write_text(json.dumps(value))
    if problem == "missing_artifact":
        (payload / "native_boundary.json").unlink()
    preflight = {"context": {"cwd": str(root)}, "saved_graph": saved}
    if problem:
        with pytest.raises((ValueError, FileNotFoundError)):
            module.optimizer_evidence(root, output, preflight, environment)
    else:
        result = module.optimizer_evidence(root, output, preflight, environment)
        assert result["genes"] == 984137 and result["groups"] == 1
        assert result["partition"] == record(partition)
        assert result["accuracy_evaluated"] is False
