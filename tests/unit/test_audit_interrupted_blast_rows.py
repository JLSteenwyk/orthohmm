import hashlib
import io
import json

import pytest

from benchmark_tools.audit_interrupted_blast_rows import ObservedPrefix
from benchmark_tools.audit_orthomcl_search_table import audit_table


def row(q):
    return q + b"\ts\t100.00\t4\t0\t0\t1\t4\t1\t4\t1e-8\t20\n"


def make(tmp_path, data, end=None, sha=None):
    path = tmp_path / "partial"
    path.write_bytes(data)
    blocks = io.StringIO()
    prefix = ObservedPrefix(path, len(data) if end is None else end,
                            sha or hashlib.sha256(data).hexdigest(), {b"q": 0, b"r": 1}, blocks)
    return prefix, blocks


def test_real_validator_and_excluded_tail(tmp_path):
    complete = row(b"q") + row(b"r")
    prefix, blocks = make(tmp_path, complete + b"r\0\0", len(complete))
    fasta, log = tmp_path / "fa", tmp_path / "log"
    fasta.write_bytes(b">q\nAAAA\n>r\nAAAA\n>s\nAAAA\n")
    log.write_text("")
    assert audit_table(prefix, fasta, log)["hsp_rows"] == 2
    inventory = [json.loads(line) for line in blocks.getvalue().splitlines()]
    assert [x["final_observed_query"] for x in inventory] == [False, True]
    assert all(not x["reuse_authorized"] for x in inventory)
    assert prefix.summary["last_query_start"] == len(row(b"q"))
    assert prefix.summary["excluded_tail_bytes"] == 3


@pytest.mark.parametrize("data", [row(b"r") + row(b"q"), row(b"q") + row(b"r") + row(b"q"),
    b"unknown\tvalue\n", b"q\0\n", b"q\ttruncated", b""])
def test_bad_order_or_bytes(tmp_path, data):
    prefix, _ = make(tmp_path, data)
    with pytest.raises(ValueError):
        with prefix.open("rb") as stream:
            list(stream)


def test_checksum_mismatch(tmp_path):
    prefix, _ = make(tmp_path, row(b"q"), sha="0" * 64)
    with pytest.raises(ValueError, match="checksum"):
        with prefix.open("rb") as stream:
            list(stream)


def test_numeric_error_propagates(tmp_path):
    prefix, _ = make(tmp_path, row(b"q").replace(b"100.00", b"nan"))
    fasta, log = tmp_path / "fa", tmp_path / "log"
    fasta.write_bytes(b">q\nAAAA\n>s\nAAAA\n")
    log.write_text("")
    with pytest.raises(ValueError, match="range"):
        audit_table(prefix, fasta, log)
