import pytest

from benchmark_tools.diagnose_corrected_vgnc_influence import contrast, removals
from benchmark_tools.diagnose_vgnc_block_influence import deleted_scores


def test_incident_cells_count_once_per_endpoint(tmp_path):
    path = tmp_path / "cells.tsv"
    path.write_text("block_left\tblock_right\tTP\tFP\tFN\nA\tA\t2\t1\t1\nA\tB\t0\t4\t0\nB\tB\t1\t0\t1\n")
    total, removed = removals(path, {"A": {"asserted_pairs": 3}, "B": {"asserted_pairs": 2}})
    assert total == dict(TP=3, FP=5, FN=2)
    assert removed == {"A": dict(TP=2, FP=5, FN=1), "B": dict(TP=1, FP=4, FN=1)}
    assert deleted_scores(total, removed["A"])[1]["f1"] == pytest.approx(2/3)
    assert deleted_scores(total, removed["B"])[1]["f1"] == pytest.approx(2/3)


@pytest.mark.parametrize("row", ["A\tB\t1\t4\t0", "A\tA\t-2\t1\t1", "A\tC\t0\t1\t0"])
def test_invalid_cells_rejected(tmp_path, row):
    path = tmp_path / "cells.tsv"
    path.write_text("block_left\tblock_right\tTP\tFP\tFN\n"+row+"\n")
    with pytest.raises(ValueError):
        removals(path, {"A": {"asserted_pairs": 3}, "B": {"asserted_pairs": 2}})


def test_all_signs_and_deterministic_extrema():
    result = contrast({"A": .2, "B": .5, "C": .7}, {"A": .5, "B": .5, "C": .5}, -.1)
    assert [result[k] for k in ("negative", "zero", "positive")] == [1, 1, 1]
    assert result["minimum_block"] == "A"
    assert result["maximum_block"] == "C"
    assert result["full_difference"] == -.1


def test_mismatched_blocks_rejected():
    with pytest.raises(ValueError):
        contrast({"A": .5}, {"B": .5}, 0.)
