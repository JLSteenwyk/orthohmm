import pytest

from benchmark_tools.search_decision_trace import watched_candidate_arrays


def test_sorted_indices_and_empty_query_offsets():
    ids, offsets = watched_candidate_arrays(["b", "a", "c"], ["z", "y"],
                                            [("c", "y"), ("b", "y"), ("b", "z")])
    assert ids.tolist() == [0, 1, 1]
    assert offsets.tolist() == [0, 2, 2, 3]


@pytest.mark.parametrize("queries,targets,pairs", [
    (["a", "a"], ["b"], []), (["a"], ["b", "b"], []),
    (["a"], ["b"], [("a", "b"), ("a", "b")]),
    (["a"], ["b"], [("missing", "b")]),
    (["a"], ["b"], [("a", "missing")]),
])
def test_invalid_candidates(queries, targets, pairs):
    with pytest.raises(ValueError):
        watched_candidate_arrays(queries, targets, pairs)
