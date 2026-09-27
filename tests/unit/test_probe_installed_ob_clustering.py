import numpy as np
import pytest

from benchmark_tools.probe_installed_ob_clustering import validate_graph, ARMS


def test_planned_repeat_and_single_score_contrast():
    assert len(ARMS) == 3
    assert ARMS[0][1] == ARMS[2][1]
    assert ARMS[0][1] != ARMS[1][1]


def test_valid_graph():
    validate_graph(["a", "b", "c"], np.array([0, 1]), np.array([1, 2]), np.array([1., 2.]))


@pytest.mark.parametrize("q,t,w", [([1], [0], [1.]), ([0], [3], [1.]),
    ([-1], [1], [1.]), ([0], [0], [1.]), ([0, 0], [1, 1], [1., 1.]),
    ([1, 0], [2, 1], [1., 1.]), ([0], [1], [float("nan")]),
    ([0], [1], [-1.]), ([0], [], [1.])])
def test_invalid_graph(q, t, w):
    with pytest.raises(ValueError):
        validate_graph(["a", "b", "c"], np.array(q), np.array(t), np.array(w))


@pytest.mark.parametrize("names", [["a", "a"], ["a", "b c"], []])
def test_invalid_names(names):
    with pytest.raises(ValueError):
        validate_graph(names, np.array([0]), np.array([1]), np.array([1.]))
