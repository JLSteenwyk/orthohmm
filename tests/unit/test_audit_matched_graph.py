import hashlib
import json

import numpy as np
import pytest

from benchmark_tools.audit_matched_graph import checkpoint, partition, reference_edges, verify_edges


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
