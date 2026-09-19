import numpy as np
import pytest

from benchmark_tools.audit_ob_search_decisions import reconstruct


def raw():
    return dict(query_ids=np.array(["q"]), target_ids=np.array(["a", "b", "c"]),
                query_indices=np.array([0, 0]), target_indices=np.array([0, 1]),
                scores=np.array([8., -1.]), evalues=np.array([0., 1e-4]),
                candidate_count=np.array(2))


def watched():
    return {("q", t): dict(historical_present=t != "b", families={"f2", "f1"}) for t in "abc"}


def test_independent_reconstruction():
    rows = reconstruct(raw(), ["q"], list("abc"), watched())
    assert [r["decision"] for r in rows] == [
        "accepted", "scored_not_significant", "not_selected_by_prefilter"]
    assert rows[2]["score"] is None and rows[2]["historical_present"] is True
    assert all(r["families"] == "f1,f2" for r in rows)


@pytest.mark.parametrize("key,value", [
    ("candidate_count", np.array(1)),
    ("candidate_count", np.array(2.)),
    ("query_ids", np.array(["different"])),
    ("target_ids", np.array(["b", "a", "c"])),
    ("query_indices", np.array([0., 0.])),
    ("query_indices", np.array([-1, 0])),
    ("target_indices", np.array([0, 3])),
    ("target_indices", np.array([0, 0])),
    ("scores", np.array([np.nan, 0.])),
    ("scores", np.array([0.])),
    ("evalues", np.array([np.inf, 0.])),
    ("evalues", np.array([-1., 0.])),
])
def test_reject_corrupt_raw_arrays(key, value):
    data = raw()
    data[key] = value
    with pytest.raises(ValueError):
        reconstruct(data, ["q"], list("abc"), watched())
