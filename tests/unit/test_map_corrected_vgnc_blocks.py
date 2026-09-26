import pytest

from benchmark_tools.map_corrected_vgnc_blocks import aggregate, write_cells


def test_cross_block_fp_preserved_once_and_truth_merged(tmp_path):
    mapping = {"a": "a", "b": "a", "c": "c"}
    rows = [["p", "q", "TP", "a", "b", "S1", "S2"],
            ["p", "r", "FP", "a", "c", "S1", "S2"],
            ["s", "t", "FN", "c", "c", "S1", "S2"]]
    result = aggregate(rows, mapping)
    assert result["a", "a"]["TP"] == 1
    assert result["a", "c"]["FP"] == 1
    assert result["c", "c"]["FN"] == 1
    output = tmp_path / "cells.tsv"
    write_cells(output, result)
    assert "a\tc\t0\t1\t0\n" in output.read_text()
    with pytest.raises(FileExistsError):
        write_cells(output, result)


def test_category_overlap_is_not_silently_deduplicated():
    rows = [["p", "q", category, "a", "b", "S1", "S2"] for category in ("TP", "FP")]
    result = aggregate(rows, {"a": "a", "b": "a"})
    assert result["a", "a"] == {"TP": 1, "FP": 1}


@pytest.mark.parametrize("rows", [
    [["p", "q", "TP", "a", "b", "S1", "S2"]],
    [["p", "p", "FP", "a", "b", "S1", "S2"]],
    [["p", "q", "bad", "a", "b", "S1", "S2"]],
    [["p", "q", "FP", "missing", "b", "S1", "S2"]],
    [["p", "q"]],
    [["p", "q", "FP", "a", "b", "S1", "S2"],
     ["q", "p", "FP", "b", "a", "S2", "S1"]],
])
def test_invalid_rows_rejected(rows):
    with pytest.raises(ValueError):
        aggregate(rows, {"a": "a", "b": "b"})
