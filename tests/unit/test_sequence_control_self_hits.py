import sqlite3

import numpy as np
import pytest

from benchmark_tools.convert_sequence_search_control import initialize, ingest, selected_sql
from orthohmm.accuracy import build_rbnh_edges, build_singleton_assignment_edges
from orthohmm.refinement import _cluster_pair_arrays


@pytest.mark.parametrize("self_score,expected_self", [(200., True), (100., False), (50., False)])
def test_self_hits_compete_for_cap_without_special_retention(self_score, expected_self):
    # Query ID 100 sorts after 100 tied same-species alternatives.
    rows = [(100, target, 0, 100., 1.) for target in range(100)] + [(100, 100, 0, self_score, self_score / 100)]
    with sqlite3.connect(":memory:") as database:
        initialize(database)
        ingest(database, reversed(rows))
        full = list(database.execute(selected_sql(None)))
        capped = list(database.execute(selected_sql(100)))
    assert len(full) == 101 and len(capped) == 100
    assert any(q == t for q, t, _ in full)
    assert any(q == t for q, t, _ in capped) is expected_self
    assert sum(q != t for q, t, _ in capped) == (99 if expected_self else 100)


def test_lexically_early_tied_self_hit_consumes_one_slot():
    with sqlite3.connect(":memory:") as database:
        initialize(database)
        ingest(database, [(0, target, 0, 100., 1.) for target in reversed(range(101))])
        capped = list(database.execute(selected_sql(100)))
    assert [target for _, target, _ in capped] == list(range(100))


def arrays(include_self):
    q, t, s = [0, 1, 2, 0], [1, 0, 0, 2], [2., 2., 1., 1.]
    if include_self:
        q += [0, 1, 2]
        t += [0, 1, 2]
        s += [100., 200., 300.]
    return np.asarray(q), np.asarray(t), np.asarray(s)


def assert_edges_equal(a, b):
    for name in ("sources", "targets", "weights"):
        np.testing.assert_array_equal(getattr(a, name), getattr(b, name))


def test_rbnh_ignores_direct_self_rows():
    names, species = ["a", "b", "c"], [0, 1, 0]
    assert_edges_equal(build_rbnh_edges(names, species, *arrays(False)),
                       build_rbnh_edges(names, species, *arrays(True)))


def test_singleton_assignment_ignores_direct_self_rows():
    names, clusters = ["a", "b", "c"], [[0, 1], [2]]
    assert_edges_equal(build_singleton_assignment_edges(names, clusters, *arrays(False)),
                       build_singleton_assignment_edges(names, clusters, *arrays(True)))


def test_cross_cluster_refinement_ignores_direct_self_rows():
    first = _cluster_pair_arrays([[0, 1], [2]], *arrays(False), total_genes=3)
    second = _cluster_pair_arrays([[0, 1], [2]], *arrays(True), total_genes=3)
    for a, b in zip(first, second):
        np.testing.assert_array_equal(a, b)
