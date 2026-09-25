import hashlib

import pytest

from benchmark_tools.trace_native_residue_hits import scan


def row(q, s):
    return f"{q}\t{s}\t100\t10\t0\t0\t1\t10\t1\t10\t1e-8\t40\n".encode()


def test_directionality_hsps_reciprocity_and_exact_identifiers(tmp_path):
    data = row("A", "A") + row("A", "B") * 2 + row("B", "A") + row("C", "A") + row("AA", "D")
    table = tmp_path / "all.m8"
    table.write_bytes(data)
    output = tmp_path / "selected.m8"
    result = scan(table, {"A", "B", "absent"}, output, hashlib.sha256(data).hexdigest(), 6)
    assert result["retained_hsp_rows"] == 5
    assert output.read_bytes() == data[:-len(row("AA", "D"))]
    a, b, absent = result["proteins"]
    assert a == dict(gene="A", outgoing_hsps=3, incoming_hsps=3, self_hsps=1,
                     outgoing_partners=["A", "B"], incoming_partners=["A", "B", "C"],
                     reciprocal_partners=["A", "B"])
    assert b["incoming_hsps"] == 2
    assert absent["outgoing_partners"] == absent["incoming_partners"] == []


@pytest.mark.parametrize("wrong_hash,wrong_rows", [(True, False), (False, True)])
def test_reject_changed_identity(tmp_path, wrong_hash, wrong_rows):
    data = row("A", "B")
    table = tmp_path / "all.m8"
    table.write_bytes(data)
    with pytest.raises(ValueError, match="digest or row count"):
        scan(table, {"A"}, tmp_path / "selected.m8",
             "0" * 64 if wrong_hash else hashlib.sha256(data).hexdigest(), 2 if wrong_rows else 1)


@pytest.mark.parametrize("data", [b"bad\n", b"A\tB\tbad\n"])
def test_reject_malformed_rows(tmp_path, data):
    table = tmp_path / "all.m8"
    table.write_bytes(data)
    with pytest.raises(ValueError, match="Invalid"):
        scan(table, {"A"}, tmp_path / "selected.m8", hashlib.sha256(data).hexdigest(), 1)


def test_never_overwrite_subset(tmp_path):
    table = tmp_path / "all.m8"
    table.write_bytes(row("A", "B"))
    output = tmp_path / "selected.m8"
    output.write_bytes(b"preserve")
    with pytest.raises(FileExistsError):
        scan(table, {"A"}, output, "0" * 64, 1)
    assert output.read_bytes() == b"preserve"
