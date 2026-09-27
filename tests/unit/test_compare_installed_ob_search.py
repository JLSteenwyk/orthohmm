import numpy as np
import pytest

from benchmark_tools.compare_installed_ob_search import compare_hits, partition


def run(old, q, t, s, names=("a", "b", "c")):
    return compare_hits(old, names, np.array(q, dtype=int), np.array(t, dtype=int), np.array(s))


def test_reordering_duplicates_and_self_hits():
    result = run({("a", "b"): 3.0, ("c", "a"): 2.0},
                 [2, 0, 0, 0], [0, 1, 0, 1], [2., 1., 9., 3.])
    assert result["normalized_scores_exact"]
    assert result["counts"]["fresh_duplicate_nonself_rows"] == 1
    assert result["counts"]["fresh_self_rows_excluded"] == 1


def test_presence_and_score_changes():
    result = run({("a", "b"): 1., ("b", "a"): 2.}, [0, 2], [1, 1], [1.5, 3.])
    assert result["counts"]["shared_changed"] == 1
    assert result["counts"]["historical_only"] == 1
    assert result["counts"]["fresh_only"] == 1
    assert result["max_shared_absolute_difference"] == .5
    assert not result["directed_nonself_presence_equal"]


def test_tolerance_is_not_exactness():
    result = run({("a", "b"): 1.}, [0], [1], [1. + 1e-13])
    assert result["directed_nonself_presence_equal"]
    assert not result["normalized_scores_exact"]
    assert result["counts"].get("shared_outside_1e12_tolerance", 0) == 0


@pytest.mark.parametrize("q,t,s", [([3], [1], [1.]), ([-1], [1], [1.]),
    ([0], [1], [float("nan")]), ([0], [1], [0.]), ([0], [], [1.])])
def test_invalid_arrays(q, t, s):
    with pytest.raises(ValueError):
        run({}, q, t, s)


@pytest.mark.parametrize("old", [{("z", "b"): 1.}, {("a", "a"): 1.},
                                    {("a", "b"): float("inf")}])
def test_invalid_historical(old):
    with pytest.raises(ValueError):
        run(old, [], [], [])


@pytest.mark.parametrize("text", ["a a\nb\n", "a\na b\n", "a\n", "a b z\n", "\na b\n"])
def test_partition_rejects_invalid(tmp_path, text):
    path = tmp_path / "groups"
    path.write_text(text)
    with pytest.raises(ValueError):
        partition(path, {"a", "b"})


def test_partition_order_independent(tmp_path):
    path = tmp_path / "groups"
    path.write_text("b a\nc\n")
    assert set(partition(path, {"a", "b", "c"})) == {frozenset(("a", "b")), frozenset(("c",))}
