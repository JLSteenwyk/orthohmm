import io
import json

import pytest

from benchmark_tools import audit_three_kingdoms_provider_checksums as module


def test_parse():
    assert module.parse_checksums("22850 6967 input.fa.gz\n") == {
        "input.fa.gz": dict(bsd_checksum=22850, blocks_1024=6967)}


@pytest.mark.parametrize("text", ["", "x 3 file", "65536 3 file", "1 -1 file", "1 2 ../file",
                                  "1 2 file\n1 2 file", "<html>unavailable</html>"])
def test_invalid(text):
    with pytest.raises(ValueError):
        module.parse_checksums(text)


@pytest.mark.parametrize("url", ["http://ftp.ebi.ac.uk/file", "https://example.org/file",
                                 "https://ftp.ebi.ac.uk/current_release/file"])
def test_untrusted_or_moving_source(url):
    with pytest.raises(ValueError):
        module.checksum_url(url)


def test_offline_audit_matches_and_skips_moving(tmp_path, monkeypatch):
    compressed = tmp_path / "local.gz"
    compressed.write_bytes(b"abc")
    import subprocess
    values = subprocess.check_output(["/usr/bin/sum", "-r", str(compressed)], text=True).split()
    source = tmp_path / "source.json"
    source.write_text(json.dumps(dict(inputs=[dict(code=code, files=dict(compressed=module.record(compressed)),
        source_url="https://ftp.ebi.ac.uk/release-47/input.fa.gz", moving_release_url=moving)
        for code, moving in (("fixed", False), ("moving", True))])))
    class Response(io.BytesIO):
        url = "https://ftp.ebi.ac.uk/release-47/CHECKSUMS"
        headers = {}
    calls = []
    def fetch(url, timeout):
        calls.append(url)
        return Response(f"{values[0]} {values[1]} input.fa.gz\n".encode())
    monkeypatch.setattr(module, "urlopen", fetch)
    output = tmp_path / "out"
    result = module.audit(source, module.record(source)["sha256"], output)
    assert [r["status"] for r in result["rows"]] == ["matched", "unresolved_moving_release"]
    assert len(calls) == 1
    assert result["redistribution_cleared"] is False
    with pytest.raises(FileExistsError):
        module.audit(source, module.record(source)["sha256"], output)


@pytest.mark.parametrize("data,status", [(b"abc", True), (b"xyz", False), (b"ab", False)])
def test_fresh_download_compares_exact_bytes(tmp_path, monkeypatch, data, status):
    expected = tmp_path / "retained"
    expected.write_bytes(b"abc")
    class Response(io.BytesIO):
        url = "https://ftp.ebi.ac.uk/input.fa.gz"
        headers = {}
    monkeypatch.setattr(module, "urlopen", lambda *a, **k: Response(data))
    destination = tmp_path / "fresh"
    result = module.reacquire(Response.url, module.record(expected), destination)
    assert result["exact_retained_bytes"] is status
    assert expected.read_bytes() == b"abc"
    with pytest.raises(FileExistsError):
        module.reacquire(Response.url, module.record(expected), destination)


def test_download_size_bound_retains_partial(tmp_path, monkeypatch):
    class Response(io.BytesIO):
        url = "https://ftp.ebi.ac.uk/input.fa.gz"
        headers = {}
    monkeypatch.setattr(module, "urlopen", lambda *a, **k: Response(b"oversized"))
    destination = tmp_path / "fresh"
    with pytest.raises(ValueError, match="size or time bound"):
        module.reacquire(Response.url, dict(bytes=3), destination)
    assert destination.read_bytes() == b"over"
