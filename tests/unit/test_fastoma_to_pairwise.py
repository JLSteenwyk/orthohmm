import gzip
import io

import pytest

from benchmark_tools.fastoma_to_pairwise import input_owners, iter_pairs, write_pairs


@pytest.mark.parametrize("compressed", [False, True])
def test_orientation_order_and_multiplicity(compressed, tmp_path):
    path = tmp_path / ("pairs.tsv.gz" if compressed else "pairs.tsv")
    opener = gzip.open if compressed else open
    with opener(path, "wt") as stream:
        stream.write("b\ta\r\na\tb\nc\tb\n")
    output = io.StringIO()
    assert write_pairs(path, {"a": "s1", "b": "s2", "c": "s1"}, output) == 3
    assert output.getvalue() == "a\tb\na\tb\nb\tc\n"


@pytest.mark.parametrize("text", ["\n", "a\n", "a\tb\tc\n", "\tb\n", "a \tb\n",
                                  "foreign\tb\n", "a\ta\n", "a\tc\n"])
def test_rejects_every_invalid_row(text, tmp_path):
    path = tmp_path / "pairs.tsv"
    path.write_text(text)
    with pytest.raises(ValueError):
        list(iter_pairs(path, {"a": "s1", "b": "s2", "c": "s1"}))


def test_empty_output_rejected(tmp_path):
    path = tmp_path / "pairs.tsv"
    path.touch()
    with pytest.raises(ValueError, match="Empty"):
        write_pairs(path, {}, io.StringIO())


def test_truncated_gzip_rejected(tmp_path):
    path = tmp_path / "pairs.gz"
    path.write_bytes(gzip.compress(b"a\tb\n")[:-5])
    with pytest.raises(EOFError):
        list(iter_pairs(path, {"a": "s1", "b": "s2"}))


def test_accession_ownership(tmp_path):
    first, second = tmp_path / "s1.fasta", tmp_path / "s2.fasta"
    first.write_text(">sp|a|name details\nACDE\n")
    second.write_text(">b\nACDE\n")
    assert input_owners([first, second]) == {"a": "s1", "b": "s2"}
    second.write_text(">tr|a|other\nACDE\n")
    with pytest.raises(ValueError, match="duplicate"):
        input_owners([first, second])


def test_empty_species_rejected(tmp_path):
    path = tmp_path / "empty.fasta"
    path.touch()
    with pytest.raises(ValueError, match="Empty species"):
        input_owners([path])
