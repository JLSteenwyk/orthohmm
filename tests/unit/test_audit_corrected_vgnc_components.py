import pytest

from benchmark_tools.audit_corrected_vgnc_components import components, read_cells, read_reference


REF = {"A": dict(proteins=2, asserted_pairs=1), "B": dict(proteins=3, asserted_pairs=3),
       "C": dict(proteins=2, asserted_pairs=1), "D": dict(proteins=2, asserted_pairs=1)}


def test_union_links_create_transitive_dependency():
    first = {("A", "B")}
    second = {("B", "C")}
    assert components(REF, first)["largest_component_blocks"] == 2
    result = components(REF, first | second)
    assert result["largest_component_blocks"] == 3
    assert result["largest_component_reference_pairs"] == 5
    assert result["isolated_blocks"] == 1
    assert result["component_size_histogram"] == {1: 1, 3: 1}


def test_no_links_preserves_all_blocks():
    result = components(REF, set())
    assert result["connected_components"] == result["isolated_blocks"] == 4
    assert result["linked_reference_pairs"] == 0


def table(tmp_path, extra=""):
    path = tmp_path / "cells.tsv"
    path.write_text("block_left\tblock_right\tTP\tFP\tFN\n"
                    "A\tA\t1\t2\t0\nB\tB\t2\t0\t1\nC\tC\t0\t0\t1\nD\tD\t1\t0\t0\n" + extra)
    return path


def test_diagonal_fp_not_graph_link(tmp_path):
    links, counts = read_cells(table(tmp_path, "A\tB\t0\t8\t0\n"), REF)
    assert links == {("A", "B")}
    assert counts == dict(TP=4, FP=10, FN=2)


@pytest.mark.parametrize("extra", ["A\tB\t1\t0\t0\n", "B\tA\t0\t1\t0\n",
    "A\tZ\t0\t1\t0\n", "A\tB\t0\t-1\t0\n", "A\tB\t0\t0\t0\n",
    "A\tB\t0\t1\t0\nA\tB\t0\t2\t0\n"])
def test_invalid_cell(tmp_path, extra):
    with pytest.raises(ValueError):
        read_cells(table(tmp_path, extra), REF)


def test_missing_reference_truth_rejected(tmp_path):
    changed = dict(REF, A=dict(proteins=2, asserted_pairs=2))
    with pytest.raises(ValueError, match="truth counts changed"):
        read_cells(table(tmp_path), changed)


def test_duplicate_reference_rejected(tmp_path):
    path = tmp_path / "ref"
    path.write_text("block\tproteins\tasserted_pairs\tspecies\nA\t2\t1\tx,y\nA\t2\t1\tx,y\n")
    with pytest.raises(ValueError, match="duplicate"):
        read_reference(path)
