import numpy as np
import pytest

from benchmark_tools.probe_ob_candidate_order_scores import factorial_arrays, LABELS
from benchmark_tools.probe_ob_canonical_candidates import canonical_factorial, array_identity


def inputs():
    old = np.array([1,0]),np.array([0,1]),np.array([2.,1.])
    fresh = np.array([0,0,1]),np.array([1,0,0]),np.array([1.1,9.,2.1])
    return factorial_arrays(2,old,fresh)[0]


def test_canonical_factorial_preserves_score_contrast():
    actual,evidence = canonical_factorial(["a","b"],inputs())
    assert all(evidence["checks"].values()) and evidence["self_hits"] == 1
    assert evidence["canonical_arrays"][LABELS[0]] == evidence["canonical_arrays"][LABELS[2]]
    assert evidence["canonical_arrays"][LABELS[1]] == evidence["canonical_arrays"][LABELS[3]]
    assert evidence["canonical_arrays"][LABELS[0]] != evidence["canonical_arrays"][LABELS[1]]
    assert actual[LABELS[0]][2].tolist() == [1.,2.]
    assert actual[LABELS[1]][2].tolist() == [1.1,2.1]


@pytest.mark.parametrize("problem", ["inventory","historical_scores","fresh_scores","self_control"])
def test_invalid_factorial_identity(problem):
    arrays = {k:tuple(a.copy() for a in v) for k,v in inputs().items()}
    if problem == "inventory":
        arrays.pop(LABELS[0])
    else:
        label = {"historical_scores":LABELS[2],"fresh_scores":LABELS[3],"self_control":LABELS[4]}[problem]
        arrays[label][2][0] += .5
    with pytest.raises(ValueError):
        canonical_factorial(["a","b"],arrays)


def test_identity_tracks_values_dtype_shape_and_handles_noncontiguous():
    a = np.arange(6,dtype=np.int32)[::2]
    assert array_identity(a) == array_identity(a.copy())
    assert array_identity(a) != array_identity(a.astype(np.int64))
    assert array_identity(a) != array_identity(a.reshape(1,3))
    assert array_identity(a) != array_identity(a+1)
