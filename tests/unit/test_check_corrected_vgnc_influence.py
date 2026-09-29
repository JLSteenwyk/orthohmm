from fractions import Fraction

import pytest

from benchmark_tools.check_corrected_vgnc_influence import close, ratios, reconstruct


def test_exact_ratios():
    assert ratios([3, 2, 1]) == (Fraction(3, 5), Fraction(3, 4), Fraction(2, 3))


@pytest.mark.parametrize("counts", [[-1, 0, 0], [0, 0, 0], [True, 1, 1], [1., 1, 1]])
def test_invalid_counts(counts):
    with pytest.raises(ValueError):
        ratios(counts)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), True, "0.5", .6])
def test_comparison_fails_closed(value):
    with pytest.raises(ValueError):
        close(value, Fraction(1, 2))


def test_independent_incident_accumulation(tmp_path):
    path = tmp_path / "cells.tsv"
    path.write_text("block_left\tblock_right\tTP\tFP\tFN\nA\tA\t2\t1\t1\nA\tB\t0\t4\t0\nB\tB\t1\t0\t1\n")
    total, full, deleted = reconstruct(path, {"A", "B"})
    assert total == [3, 5, 2]
    assert full[2] == Fraction(6, 13)
    assert deleted["A"] == ([2, 5, 1], (Fraction(1), Fraction(1, 2), Fraction(2, 3)))
    assert deleted["B"][0] == [1, 4, 1]


@pytest.mark.parametrize("row", ["A\tB\t1\t1\t0", "A\tC\t0\t1\t0", "B\tA\t0\t1\t0"])
def test_invalid_sparse_cell(tmp_path, row):
    path = tmp_path / "cells.tsv"
    path.write_text("block_left\tblock_right\tTP\tFP\tFN\n" + row + "\n")
    with pytest.raises(ValueError):
        reconstruct(path, {"A", "B"})
