import hashlib
import json

import numpy as np
import pytest

from benchmark_tools.audit_matched_graph import checkpoint, execution_bindings, partition, reference_edges, verify_edges
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def numeric():
    return dict(gene_names=["a", "b", "c"], gene_to_species=[0, 1, 2],
                hit_queries=[0, 1, 2, 2], hit_targets=[1, 0, 0, 1], hit_scores=[3., 5., 2., 1.])


def test_independent_reciprocal_threshold_and_singleton_attachment():
    rbnh, multipass = reference_edges(numeric(), [["a", "b"], ["c"]])
    assert rbnh == {(0, 1): 5.}
    assert multipass == {(0, 1): 5., (0, 2): 2., (1, 2): 1.}


def test_empty_edges_and_self_hits():
    data = dict(gene_names=["a"], gene_to_species=[0], hit_queries=[0], hit_targets=[0], hit_scores=[4])
    assert reference_edges(data, [["a"]]) == ({}, {})


def test_partition_coverage(tmp_path):
    path = tmp_path / "partition"
    path.write_text("a\tb\nc\n")
    assert partition(path, ["a", "b", "c"]) == [["a", "b"], ["c"]]
    for text in ("a\tb\nb\tc\n", "a\tb\n", "a\tb\nunknown\n", "a\tb\nc\n\n"):
        path.write_text(text)
        with pytest.raises(ValueError):
            partition(path, ["a", "b", "c"])


def test_edge_weights_and_duplicates(tmp_path):
    path = tmp_path / "edges.npz"
    np.savez(path, sources=np.array([0], dtype=np.int32), targets=np.array([1], dtype=np.int32), weights=np.array([5.]))
    verify_edges(path, {(0, 1): 5.})
    with pytest.raises(ValueError):
        verify_edges(path, {(0, 1): 4.})
    np.savez(path, sources=np.array([0, 0], dtype=np.int32), targets=np.array([1, 1], dtype=np.int32), weights=np.array([5., 5.]))
    with pytest.raises(ValueError):
        verify_edges(path, {(0, 1): 5.})


def write_checkpoint(path, data):
    (path / "gene_names.txt").write_text("\n".join(data["gene_names"]) + "\n")
    for key in data.keys() - {"gene_names"}:
        np.save(path / (key + ".npy"), np.asarray(data[key], dtype=np.float64 if key == "hit_scores" else np.int32))
    files = {p.name: dict(bytes=p.stat().st_size, sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in path.iterdir()}
    (path / "manifest.json").write_text(json.dumps(dict(complete=True, genes=3, hits=4, files=files)))


def test_checkpoint_exact_mapping(tmp_path):
    data = numeric()
    write_checkpoint(tmp_path, data)
    checkpoint(tmp_path, data)
    altered = dict(data, hit_scores=[3., 5., 2., 2.])
    with pytest.raises(ValueError, match="numeric array"):
        checkpoint(tmp_path, altered)
    (tmp_path / "extra").write_text("")
    with pytest.raises(ValueError, match="inventory"):
        checkpoint(tmp_path, data)


def test_checkpoint_tamper(tmp_path):
    data = numeric()
    write_checkpoint(tmp_path, data)
    np.save(tmp_path / "hit_scores.npy", np.array([1., 2., 3., 4.]))
    with pytest.raises(ValueError, match="checksum"):
        checkpoint(tmp_path, data)


@pytest.fixture
def bindings(tmp_path):
    python = tmp_path / "venv/bin/python"
    for name in ("numeric.json", "graph/receipt.json", "graph.log", "graph.time.txt",
                 "graph/rbnh_edges.npz", "graph/multipass_edges.npz", "graph/initial.tsv",
                 "graph/multipass.tsv", "graph/final.tsv",
                 "graph/orthohmm_working_res/high_sensitivity_checkpoint/manifest.json",
                 "venv/lib/orthohmm/accuracy.py", "venv/lib/orthohmm/externals.py",
                 "venv/lib/orthohmm/refinement.py"):
        file = tmp_path / name
        file.parent.mkdir(parents=True, exist_ok=True)
        file.write_text("fixture")
    modules = [record(tmp_path / "venv/lib/orthohmm" / name)
               for name in ("accuracy.py", "externals.py", "refinement.py")]
    execution = dict(environment=dict(PATH="/usr/bin:/bin", HOME=str(tmp_path), LC_ALL="C",
                                     OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                                     MKL_NUM_THREADS="1", PYTHONHASHSEED="0"),
                     private_numeric=record(tmp_path / "numeric.json"),
                     native_receipt=record(tmp_path / "graph/receipt.json"),
                     logs=[record(tmp_path / name) for name in ("graph.log", "graph.time.txt")])
    native = dict(numeric=record(tmp_path / "numeric.json"),
                  status="native_graph_completed_pending_independent_readback",
                  executable=str(python), prefix=str(python.parent.parent), modules=modules,
                  outputs=[record(tmp_path / "graph" / name) for name in
                           ("rbnh_edges.npz", "multipass_edges.npz", "initial.tsv", "multipass.tsv", "final.tsv")],
                  checkpoint_manifest=record(tmp_path / "graph/orthohmm_working_res/high_sensitivity_checkpoint/manifest.json"))
    return tmp_path, execution, native, python, dict(checked_records=list(modules))


def test_execution_bindings_valid(bindings):
    execution_bindings(*bindings)


@pytest.mark.parametrize("key", ["PATH", "HOME", "LC_ALL", "OMP_NUM_THREADS",
                                "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "PYTHONHASHSEED", "PYTHONPATH"])
def test_changed_or_extra_environment_rejected(bindings, key):
    bindings[1]["environment"][key] = "changed"
    with pytest.raises(ValueError, match="environment"):
        execution_bindings(*bindings)


@pytest.mark.parametrize("field", ["private_numeric", "native_receipt", "logs"])
def test_execution_wrong_artifact_rejected(bindings, field):
    item = bindings[1][field]
    if isinstance(item, list):
        item.pop()
    else:
        item["path"] += ".elsewhere"
    with pytest.raises(ValueError, match="binding|inventory"):
        execution_bindings(*bindings)


@pytest.mark.parametrize("field", ["numeric", "status", "executable", "prefix", "checkpoint_manifest"])
def test_native_wrong_binding_rejected(bindings, field):
    item = bindings[2][field]
    if isinstance(item, dict):
        item["path"] += ".elsewhere"
    else:
        bindings[2][field] = "different"
    with pytest.raises(ValueError, match="binding|inventory"):
        execution_bindings(*bindings)


@pytest.mark.parametrize("field", ["outputs", "modules"])
@pytest.mark.parametrize("mutation", ["missing", "duplicate", "empty", "outside"])
def test_native_inventory_rejected(bindings, field, mutation):
    items = bindings[2][field]
    if mutation == "missing":
        items.pop()
    elif mutation == "duplicate":
        items[0] = items[1]
    elif mutation == "empty":
        items.clear()
    else:
        items[0] = dict(items[0], path="/outside/accuracy.py")
        if field == "modules":
            bindings[4]["checked_records"].append(items[0])
    with pytest.raises(ValueError, match="inventory"):
        execution_bindings(*bindings)
