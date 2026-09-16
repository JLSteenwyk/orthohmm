import csv

import pytest

from benchmark_tools.verify_ygob_reference import reconstruct


def inputs(tmp_path, rows, fasta):
    table, sequences = tmp_path / "pillars", tmp_path / "aa"
    with table.open("w") as handle:
        csv.writer(handle, delimiter="\t").writerows(rows)
    sequences.write_text(fasta)
    return table, sequences


def row(cells):
    values = ["---"] * 33
    for column, gene in cells.items():
        values[column] = gene
    return values


def test_column_ownership_filters_and_terminal_stop(tmp_path):
    paths = inputs(tmp_path, [row({0: "a", 32: "b", 11: "excluded", 13: "c", 14: "off"})],
                   ">a {ON}\nacde*\n>b {ON}\nACDE\n>excluded {ON}\nACDE\n>c {ON}\nACDE\n>off {OFF}\nACDE\n")
    seqs, ref, excluded = reconstruct(*paths)
    assert seqs == {"a": ("Vpolyspora", "ACDE"), "b": ("Vpolyspora", "ACDE"), "c": ("Zrouxii", "ACDE")}
    assert ref == {"Pillar00001": ["a", "b", "c"]}
    assert excluded == []


def test_entire_duplicate_rows_excluded_but_inputs_retained(tmp_path):
    paths = inputs(tmp_path, [row({0: "a", 13: "b"}), row({0: "a", 14: "c"})],
                   ">a {ON}\nACD\n>b {ON}\nACD\n>c {ON}\nACD\n")
    seqs, ref, excluded = reconstruct(*paths)
    assert set(seqs) == {"a", "b", "c"}
    assert ref == {}
    assert excluded == [1, 2]


@pytest.mark.parametrize("fasta", [">a {ON}\nACD\n>a {ON}\nACD\n", ">a\nACD\n", ">a {ON}\nAC?\n"])
def test_invalid_fasta_rejected(tmp_path, fasta):
    with pytest.raises(ValueError):
        reconstruct(*inputs(tmp_path, [row({0: "a"})], fasta))


def test_wrong_column_count_rejected(tmp_path):
    with pytest.raises(ValueError, match="33"):
        reconstruct(*inputs(tmp_path, [["a"]], ">a {ON}\nACD\n"))
