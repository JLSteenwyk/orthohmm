import numpy as np
import pytest

from benchmark_tools.audit_accuracy_checkpoint import check_arrays


def arrays():
    return ["b", "a"], np.array([4, 9], dtype="int32"), np.array([0, 1], dtype="int32"), np.array([1, 1], dtype="int32"), np.array([2., 3.])


def test_chunked_checks_preserve_unsorted_native_order():
    data = arrays()
    result = check_arrays(*data, chunk_size=1)
    assert result["species"] == 2 and result["self_hits"] == 1
    assert result["gene_names_lexically_sorted"] is False
    assert data[0] == ["b", "a"]


@pytest.mark.parametrize("change", ["duplicate", "bounds", "negative", "shape", "dtype", "nonfinite"])
def test_rejects_invalid_native_arrays(change):
    names, species, q, t, s = arrays()
    if change == "duplicate":
        names = ["a", "a"]
    elif change == "bounds":
        q[0] = 2
    elif change == "negative":
        species[0] = -1
    elif change == "shape":
        t = t[:1]
    elif change == "dtype":
        q = q.astype("int64")
    else:
        s[1] = np.nan
    with pytest.raises(ValueError):
        check_arrays(names, species, q, t, s, chunk_size=1)
