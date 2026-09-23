import hashlib
import io

import pytest

from benchmark_tools.audit_interrupted_blast_bytes import audit, scan
from benchmark_tools import audit_interrupted_blast_bytes as module


@pytest.mark.parametrize("data,zeros,first,last,trailing", [
    (b"", 0, None, None, False),
    (b"q\ts\n", 0, None, None, False),
    (b"q\ts\nq\0\0\0", 3, 5, 7, True),
    (b"\0\0\0\0", 4, 0, 3, True),
    (b"q\0x\0", 2, 1, 3, False),
    (b"q\0x\n", 1, 1, 1, False),
])
@pytest.mark.parametrize("chunk_size", [1, 2, 3, 8, 1024])
def test_scan(data, zeros, first, last, trailing, chunk_size):
    result = scan(io.BytesIO(data), chunk_size)
    assert result["bytes"] == len(data)
    assert result["sha256"] == hashlib.sha256(data).hexdigest()
    assert result["nul_bytes"] == zeros
    assert result["first_nul_offset"] == first
    assert result["last_nul_offset"] == last
    assert result["nuls_form_single_trailing_run"] is trailing
    assert result["newline_count"] == data.count(b"\n")
    assert result["last_newline_offset"] == (data.rfind(b"\n") if b"\n" in data else None)
    assert result["ends_with_newline"] == data.endswith(b"\n")


def test_invalid_chunk_size():
    with pytest.raises(ValueError):
        scan(io.BytesIO(b"x"), 0)


def test_audit_is_read_only_and_never_admits(tmp_path):
    path, output = tmp_path / "partial", tmp_path / "audit.json"
    path.write_bytes(b"query\tsubject\npartial\0")
    original = path.read_bytes()
    report = audit(path, output)
    assert report["search_admitted"] is False
    assert report["reuse_authorized"] is False
    assert path.read_bytes() == original
    with pytest.raises(FileExistsError):
        audit(path, output)
    with pytest.raises(FileExistsError):
        audit(path, path)


def test_changed_input_is_rejected(tmp_path, monkeypatch):
    path, output = tmp_path / "partial", tmp_path / "audit.json"
    path.write_bytes(b"original\n")

    def mutate(stream):
        result = scan(stream)
        path.write_bytes(b"changed and longer\n")
        return result

    monkeypatch.setattr(module, "scan", mutate)
    with pytest.raises(ValueError, match="Input changed"):
        audit(path, output)
    assert not output.exists()
