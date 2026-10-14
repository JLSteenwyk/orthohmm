import numpy as np
import pytest

from benchmark_tools.diagnose_component_pair_structure import direct_support, feature_key, partition_membership


def test_partition_preserves_source_order_and_full_gene_universe(tmp_path):
    path = tmp_path / "groups.txt"
    path.write_text("g3 g1\ng2\n")
    members, sizes = partition_membership(path, ["g1", "g2", "g3"])
    np.testing.assert_array_equal(members, [0, 1, 0])
    np.testing.assert_array_equal(sizes, [2, 1])
    path.write_text("g1\ng2\n")
    with pytest.raises(ValueError, match="universe"):
        partition_membership(path, ["g1", "g2", "g3"])


def test_chunked_undirected_support_preserves_large_integer_keys_and_duplicate_maximum():
    n = 984137
    a, b = 900000, 900001
    wanted = np.array([1 * n + 2, a * n + b], dtype=np.int64)
    sources = np.array([2, b, a, 3], dtype=np.int32)
    targets = np.array([1, a, b, 4], dtype=np.int32)
    weights = np.array([.4, 1.6, 2.3, 5.])
    before = [arr.tobytes() for arr in (wanted, sources, targets, weights)]
    np.testing.assert_array_equal(direct_support(wanted, sources, targets, weights, n, chunk_size=1), [.4, 2.3])
    assert before == [arr.tobytes() for arr in (wanted, sources, targets, weights)]


def test_absent_and_empty_watched_pairs_are_not_imputed_as_positive_support():
    sources, targets, weights = np.array([0]), np.array([1]), np.array([2.])
    np.testing.assert_array_equal(direct_support(np.array([2], dtype=np.int64), sources, targets, weights, 3), [0.])
    assert direct_support(np.array([], dtype=np.int64), sources, targets, weights, 3).shape == (0,)


@pytest.mark.parametrize("mutation", ["shape", "mismatch", "float_endpoint", "negative", "overflow", "nan", "zero", "order", "duplicate", "bool_size", "key_limit", "chunk"])
def test_invalid_graph_data_and_key_domains_are_rejected(mutation):
    wanted, sources, targets, weights = np.array([1, 2]), np.array([0]), np.array([1]), np.array([2.])
    n, chunk = 3, 1
    if mutation == "shape":
        sources = sources.reshape(1, 1)
    elif mutation == "mismatch":
        targets = np.array([1, 2])
    elif mutation == "float_endpoint":
        sources = np.array([0.])
    elif mutation == "negative":
        sources[0] = -1
    elif mutation == "overflow":
        targets[0] = 3
    elif mutation == "nan":
        weights[0] = np.nan
    elif mutation == "zero":
        weights[0] = 0
    elif mutation == "order":
        wanted = wanted[::-1]
    elif mutation == "duplicate":
        wanted = np.array([1, 1])
    elif mutation == "bool_size":
        n = True
    elif mutation == "key_limit":
        n = 3_037_000_500
    else:
        chunk = 0
    with pytest.raises(ValueError):
        direct_support(wanted, sources, targets, weights, n, chunk)


def test_feature_classification_is_label_independent_and_threshold_is_inclusive():
    parent = np.array([0, 0])
    seed = np.array([0, 1])
    candidate = np.array([3, 3])
    retained = np.array([True, False])
    assert feature_key(0, 1, parent, seed, candidate, retained, 1.5) == (True, False, True, 1, "at_least_1_5")
    assert feature_key(0, 1, parent, seed, candidate, retained, .5)[-1] == "positive_below_1_5"
    assert feature_key(0, 1, parent, seed, candidate, retained, 0)[-1] == "absent"
