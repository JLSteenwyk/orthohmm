import gzip
from io import StringIO

import pytest

from benchmark_tools.fastoma_distinct_pairs import write_pairs


OWNERS = {"a": "species1", "b": "species2", "c": "species1", "d": "species3"}


def convert(tmp_path, text, owners=OWNERS):
    source = tmp_path / "orthologs.tsv.gz"
    source.write_bytes(gzip.compress(text.encode()))
    output = StringIO()
    stats = write_pairs(source, owners, output, tmp_path)
    assert not list(tmp_path.glob("fastoma-pairs-*"))
    return output.getvalue(), stats


def test_native_pairs_sorted_and_duplicates_counted(tmp_path):
    output, stats = convert(tmp_path, "d\tb\nb\ta\na\tb\nb\td\nc\td\n")
    assert output == "a\tb\nb\td\nc\td\n"
    assert stats == {"native_rows": 5, "distinct_pairs": 3, "duplicate_relations": 2}
    # No clique expansion: a-d and b-c were never predicted.
    assert "a\td\n" not in output and "b\tc\n" not in output


@pytest.mark.parametrize("bad", ["", "\n", "a\tb\textra\n", "a\n", "a\t\n",
                                  " a\tb\n", "a\tb x\n", "x\tb\n", "a\ta\n",
                                  "a\tc\n", "gene1\tgene2\n", "sp|a|id\tb\n"])
def test_bad_rows_fail_without_partial_output(tmp_path, bad):
    source = tmp_path / "orthologs.tsv.gz"
    source.write_bytes(gzip.compress(("a\tb\n" + bad if bad else "").encode()))
    output = StringIO()
    with pytest.raises(ValueError):
        write_pairs(source, OWNERS, output, tmp_path)
    assert output.getvalue() == ""
    assert not list(tmp_path.glob("fastoma-pairs-*"))


@pytest.mark.parametrize("text", ["b\ta", "b\ta\r\n", "b\ta\n"])
def test_line_endings(tmp_path, text):
    assert convert(tmp_path, text)[0] == "a\tb\n"


@pytest.mark.parametrize("owners", [{}, {"a": ""}, {"a b": "species"}])
def test_invalid_owner_map(tmp_path, owners):
    with pytest.raises(ValueError):
        convert(tmp_path, "a\tb\n", owners)


def test_truncated_gzip_does_not_emit_output(tmp_path):
    source = tmp_path / "orthologs.tsv.gz"
    source.write_bytes(gzip.compress(b"a\tb\n")[:-5])
    output = StringIO()
    with pytest.raises(EOFError):
        write_pairs(source, OWNERS, output, tmp_path)
    assert output.getvalue() == ""
    assert not list(tmp_path.glob("fastoma-pairs-*"))


def test_transaction_boundary_preserves_counts(tmp_path):
    output, stats = convert(tmp_path, "a\tb\n" * 100001 + "d\tc\n")
    assert output == "a\tb\nc\td\n"
    assert stats == {"native_rows": 100002, "distinct_pairs": 2,
                     "duplicate_relations": 100000}


def test_output_failure_cleans_scratch(tmp_path):
    class BrokenOutput:
        def write(self, value):
            raise OSError("full output filesystem")

    source = tmp_path / "orthologs.tsv.gz"
    source.write_bytes(gzip.compress(b"a\tb\n"))
    with pytest.raises(OSError, match="full output"):
        write_pairs(source, OWNERS, BrokenOutput(), tmp_path)
    assert not list(tmp_path.glob("fastoma-pairs-*"))
