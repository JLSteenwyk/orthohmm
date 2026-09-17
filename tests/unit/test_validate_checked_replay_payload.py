import json
from pathlib import Path

import igraph
import leidenalg
import numpy as np
import pytest

from benchmark_tools.checked_python_pair_worker import python_pair_constructor
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.probe_leiden_boundary import observe_partition
from benchmark_tools.validate_checked_replay_payload import validate
from benchmark_tools.run_qfo_checked_full_replay import run


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


@pytest.mark.parametrize("problem", [None, "boundary", "input_hash", "metadata", "module", "coverage", "duplicate", "provenance"])
def test_native_result_gate(tmp_path, problem):
    root = tmp_path
    executor = root / "executor"
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    payload = root / "stage/payload"
    payload.mkdir(parents=True)
    pairs = np.array([[2, 0], [0, 2]], dtype=np.int32)
    for name, data in (("sources", pairs[:, 0]), ("targets", pairs[:, 1]), ("weights", np.ones(2))):
        np.save(payload / (name + ".npy"), data)
    (payload / "gene_names.txt").write_text("a\nb\nc\nd\n")
    metadata = {"cpm_resolution": .1, "seed": 4, "include_isolates": True, "output_directory": str(root / "output")}
    write_json(payload / "metadata.json", metadata)
    manifest = {"stage": "initial", "inputs": [record(payload / name) for name in
                ("gene_names.txt", "sources.npy", "targets.npy", "weights.npy", "metadata.json")]}
    write_json(payload.parent / "payload_manifest.json", manifest)
    with python_pair_constructor(igraph, payload / "constructor_adapter.json"):
        graph = igraph.Graph(n=4, edges=pairs, directed=False)
        graph.es["weight"] = [1., 1.]
        with observe_partition(leidenalg, payload):
            leidenalg.find_partition(graph, leidenalg.CPMVertexPartition, weights="weight", seed=4, resolution_parameter=.1)
    modules = {}
    for name in ("orthohmm.leiden_worker", "orthohmm.externals", "orthohmm.helpers"):
        path = launcher / (name.replace(".", "/") + ".py")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("# fixture\n")
        modules[name] = record(path)
    helpers = []
    for name in ("repeat_qfo_saved_graph.py", "checked_python_pair_worker.py", "probe_leiden_boundary.py", "checked_replay_payload_worker.py"):
        path = executor / "benchmark_tools" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("# fixture\n")
        helpers.append(record(path))
    admission = root / "benchmark_tools/results/qfo_checked_repeats_verified_20260917.json"
    write_json(admission, {})
    provenance = {"source": helpers[-1], "helpers": helpers[:3], "manifest": record(payload.parent / "payload_manifest.json"),
                  "inputs": manifest["inputs"], "stage": "initial", "accuracy_evaluated": False, "admission": record(admission)}
    worker = {"status": "before_native_clustering", "accuracy_evaluated": False, "metadata": dict(metadata),
              "cwd": str(launcher), "inputs": manifest["inputs"][:4], "cpu_affinity": [0], "modules": modules,
              "native_libraries": [record(igraph._igraph.__file__)], "python": record(__import__("sys").executable), "observer": helpers[0],
              "environment": {"PYTHONPATH": str(launcher), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                              "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}}
    partition = root / "output/orthohmm_working_res/orthohmm_edges_clustered.txt"
    partition.parent.mkdir(parents=True)
    partition.write_text("a c\nb\nd\n")
    if problem == "boundary":
        boundary = json.loads((payload / "native_boundary.json").read_text())
        boundary["calls"][0]["after"]["ordered_endpoints_sha256"] = "changed"
        write_json(payload / "native_boundary.json", boundary)
    elif problem == "input_hash":
        adapter = json.loads((payload / "constructor_adapter.json").read_text())
        adapter["calls"][0]["ordered_input_bytes_sha256"] = "changed"
        write_json(payload / "constructor_adapter.json", adapter)
    elif problem == "metadata":
        worker["metadata"]["seed"] = 5
    elif problem == "module":
        worker["modules"]["orthohmm.externals"]["sha256"] = "changed"
    elif problem == "coverage":
        partition.write_text("a c\nb\n")
    elif problem == "duplicate":
        partition.write_text("a c\nb\nd a\n")
    elif problem == "provenance":
        provenance["stage"] = "multipass"
    write_json(payload / "checked_payload_provenance.json", provenance)
    write_json(payload / "worker_before.json", worker)
    if problem:
        with pytest.raises(ValueError):
            validate(payload, manifest, root, executor)
    else:
        result = validate(payload, manifest, root, executor)
        assert result["genes"] == 4 and result["groups"] == 3


def test_parent_refuses_existing_output(tmp_path):
    with pytest.raises(FileExistsError):
        run(tmp_path, tmp_path)
