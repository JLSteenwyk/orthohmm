import pytest

from benchmark_tools.validate_sonicparanoid_groups import record, validate


@pytest.fixture
def files(tmp_path):
    first, second = tmp_path / "a.fasta", tmp_path / "b.fasta"
    first.write_text(">a1\nACD\n>a2\nACE\n")
    second.write_text(">tr|b1|name\nACD\n")
    snapshot = tmp_path / "snapshot.tsv"
    snapshot.write_text(f"1\ta.fasta\t{record(first)['sha256']}\t2\t6\n"
                        f"2\tb.fasta\t{record(second)['sha256']}\t1\t3\n")
    table = tmp_path / "groups.tsv"
    table.write_text("group_id\tgroup_size\tsp_in_grp\tseed_ortholog_cnt\ta.fasta\tb.fasta\n"
                     "1\t2\t2\t2\ta1\ttr|b1|name\n")
    normalized = tmp_path / "normalized.txt"
    normalized.write_text("a1 tr|b1|name\n")
    return table, [first, second], snapshot, normalized


def test_full_identifiers_and_unassigned_count(files):
    result = validate(*files)
    assert result["grouped_genes"] == 2
    assert result["unassigned_input_genes"] == 1
    assert not result["accuracy_admitted"]


@pytest.mark.parametrize("old,new", [
    ("a.fasta\tb.fasta", "a.fasta\ta.fasta"),
    ("a1\ttr|b1|name", "tr|b1|name\ta1"),
    ("a1\ttr|b1|name", "a1,a1\ttr|b1|name"),
    ("a1\ttr|b1|name", "a1,\ttr|b1|name"),
    ("\t2\t2\t2\t", "\t3\t2\t2\t"),
    ("\t2\t2\t2\t", "\t2\t1\t2\t"),
    ("\t2\t2\t2\t", "\t2\t2\t3\t"),
    ("tr|b1|name\n", "b1\n"),
    ("tr|b1|name\n", "tr|b1|name\textra\n"),
])
def test_native_corruption_rejected(files, old, new):
    table = files[0]
    table.write_text(table.read_text().replace(old, new))
    with pytest.raises(ValueError):
        validate(*files)


@pytest.mark.parametrize("text", ["a1\ntr|b1|name\n", "a1 tr|b1|name a2\n", "a1 tr|b1|name\na2\n"])
def test_normalized_group_corruption_rejected(files, text):
    files[3].write_text(text)
    with pytest.raises(ValueError):
        validate(*files)


def test_snapshot_mismatch(files):
    files[2].write_text(files[2].read_text().replace("\t2\t6", "\t3\t6"))
    with pytest.raises(ValueError, match="Snapshot input"):
        validate(*files)
