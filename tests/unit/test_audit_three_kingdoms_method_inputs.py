import pytest

from benchmark_tools.audit_three_kingdoms_method_inputs import classify, parse_sha256, snapshot_rows


def test_sha_inventory_preserves_full_paths():
    assert parse_sha256("a" * 64 + "  /tmp/a b.fasta\n") == {"/tmp/a b.fasta": "a" * 64}


@pytest.mark.parametrize("text", ["", "z" * 64 + "  /tmp/a", "a" * 64 + " /tmp/a", ("a" * 64 + "  /tmp/a\n") * 2])
def test_bad_sha_inventory(text):
    with pytest.raises(ValueError):
        parse_sha256(text)


def test_raw_and_staged_discrimination():
    assert classify("a", "a", "a") == "matches_staged"
    assert classify("b", "a", "b") == "matches_raw_not_staged"
    assert classify("c", "a", "b") == "matches_neither"


def test_native_snapshot_parser():
    text = "1\tA.fasta\t" + "a" * 64 + "\t12\t123\n"
    assert snapshot_rows(text) == {"A.fasta": {"sha256": "a" * 64, "proteins": 12, "reported_residues": 123}}


@pytest.mark.parametrize("text", ["2\tA.fasta\t" + "a" * 64 + "\t1\t1", "1\tA.fasta\tbad\t1\t1",
                                  "1\tA.fasta\t" + "a" * 64 + "\t-1\t1", "1\tA.fasta\n"])
def test_bad_native_snapshot(text):
    with pytest.raises(ValueError):
        snapshot_rows(text)
