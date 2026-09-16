from pathlib import Path
from types import SimpleNamespace

import numpy as np

from benchmark_tools.diagnose_qfo_initial_graph import (
    CORE, OLD, array_evidence, compare_edges, load_revision,
)


def edges(weight=1., names=None):
    return SimpleNamespace(gene_names=names or ["a", "b"], sources=np.array([0], dtype=np.int32),
                           targets=np.array([1], dtype=np.int32), weights=np.array([weight]))


def test_equal_counts_do_not_hide_weight_or_name_differences():
    assert compare_edges(edges(), edges())["byte_equal"]
    assert not compare_edges(edges(), edges(np.nextafter(1., 2.)))["byte_equal"]
    assert not compare_edges(edges(), edges(names=["b", "a"]))["byte_equal"]


def test_array_dtype_and_layout_are_recorded():
    a = np.arange(8, dtype=np.int32)[::2]
    assert array_evidence(a) == array_evidence(a.copy())
    assert array_evidence(a) != array_evidence(a.astype(np.int64))


def test_historical_and_frozen_rbnh_on_ties_and_self_hits():
    root = Path(__file__).resolve().parents[2]
    old, old_source = load_revision(root, OLD, "test_old")
    new, new_source = load_revision(root, CORE, "test_new")
    names, species = ["a", "b", "c"], np.array([0, 1, 1], dtype=np.int32)
    q, t = np.array([0, 0, 0, 1, 2]), np.array([0, 1, 2, 0, 0])
    scores = np.array([9., 2., 2., 2., 1.])
    assert old_source["sha256"] != new_source["sha256"]
    assert compare_edges(old.build_rbnh_edges(names, species, q, t, scores),
                         new.build_rbnh_edges(names, species, q, t, scores))["byte_equal"]
