import numpy as np
import pytest

from benchmark_tools.probe_ob_candidate_order_scores import factorial_arrays, LABELS


def fixture():
    old = np.array([0,1]), np.array([1,2]), np.array([1.,2.])
    fresh = np.array([1,0,0]), np.array([2,0,1]), np.array([2.1,9.,1.1])
    return old, fresh


def test_factorial_isolates_order_scores_and_self_hits():
    old, fresh = fixture()
    arms, removed = factorial_arrays(3, old, fresh)
    assert tuple(arms) == LABELS and removed == 1
    expected = ([1.,2.], [1.1,2.1], [2.,1.], [2.1,1.1], [2.1,9.,1.1])
    for label, scores in zip(LABELS, expected):
        assert arms[label][2].tolist() == scores
    assert arms[LABELS[0]][0].tolist() == arms[LABELS[1]][0].tolist() == [0,1]
    assert arms[LABELS[2]][0].tolist() == arms[LABELS[3]][0].tolist() == [1,0]
    assert all(x is y for x,y in zip(arms[LABELS[4]],fresh))


@pytest.mark.parametrize("problem", ["missing", "duplicate", "historical_self", "nan"])
def test_invalid_alignment_rejected(problem):
    old, fresh = fixture()
    if problem == "missing":
        fresh[0][0] = 2
        fresh[1][0] = 1
    elif problem == "duplicate":
        fresh[0][0] = 0
        fresh[1][0] = 1
    elif problem == "historical_self":
        old[1][0] = 0
    else:
        old[2][0] = np.nan
    with pytest.raises(ValueError):
        factorial_arrays(3, old, fresh)
