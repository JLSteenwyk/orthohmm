import io
import json
import subprocess
import sys

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


@pytest.fixture
def offline_source(tmp_path, monkeypatch):
    compressed = tmp_path / "local.gz"
    compressed.write_bytes(b"abc")
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
        return Response(b"16556 1 input.fa.gz\n")
    monkeypatch.setattr(module, "urlopen", fetch)
    return source, compressed, calls


def test_offline_audit_matches_and_skips_moving(tmp_path, offline_source):
    source, compressed, calls = offline_source
    output = tmp_path / "out"
    result = module.audit(source, module.record(source)["sha256"], output)
    assert [r["status"] for r in result["rows"]] == ["matched", "unresolved_moving_release"]
    assert len(calls) == 1
    assert result["redistribution_cleared"] is False
    assert result["rows"][0]["command"] == ["/usr/bin/sum", str(compressed)]
    assert result["rows"][0]["observed"] == dict(bsd_checksum=16556, blocks_1024=1)
    with pytest.raises(FileExistsError):
        module.audit(source, module.record(source)["sha256"], output)


@pytest.mark.parametrize("data,blocks", [(b"", 0), (b"abc", 1), (b"x" * 1023, 1),
                                        (b"x" * 1024, 1), (b"x" * 1025, 2)])
def test_native_default_bsd_sum_and_block_rounding(tmp_path, data, blocks):
    path = tmp_path / "payload with spaces"
    path.write_bytes(data)
    values = subprocess.check_output(["/usr/bin/sum", str(path)], text=True).split(maxsplit=2)
    assert 0 <= int(values[0]) <= 65535 and int(values[1]) == blocks
    assert values[2].strip() == str(path)
    if sys.platform.startswith("linux"):
        assert subprocess.check_output(["/usr/bin/sum", "-r", str(path)], text=True).split(maxsplit=2) == values


@pytest.mark.parametrize("stdout,status", [
    ("16556 1 {path}\n", "matched"), ("16555 1 {path}\n", "mismatch"),
    ("16556", "unresolved"), ("65536 1 {path}\n", "unresolved"),
    ("16556 2 {path}\n", "unresolved"), ("16556 -1 {path}\n", "unresolved"),
    ("16556 1 wrong-file\n", "unresolved"), ("not-a-sum 1 {path}\n", "unresolved")])
def test_native_command_output_validation(tmp_path, monkeypatch, offline_source, stdout, status):
    source, compressed, calls = offline_source
    commands = []
    def run(command, **kwargs):
        commands.append(command)
        assert kwargs == dict(capture_output=True, text=True, check=True)
        return subprocess.CompletedProcess(command, 0, stdout.format(path=compressed), "")
    monkeypatch.setattr(module.subprocess, "run", run)
    result = module.audit(source, module.record(source)["sha256"], tmp_path / "out")
    assert [r["status"] for r in result["rows"]] == [status, "unresolved_moving_release"]
    assert result["rows"][0]["command_stdout"] == stdout.format(path=compressed)
    assert result["rows"][0]["expected"] == dict(bsd_checksum=16556, blocks_1024=1)
    assert commands == [["/usr/bin/sum", str(compressed)]] and len(calls) == 1
    assert compressed.read_bytes() == b"abc"


@pytest.mark.parametrize("error", [subprocess.CalledProcessError(1, ["sum"], output="partial", stderr="failure"),
                                  PermissionError("not executable")])
def test_checksum_command_failure_is_retained_without_retry(tmp_path, monkeypatch, offline_source, error):
    source, compressed, calls = offline_source
    commands = []
    def run(command, **kwargs):
        commands.append(command)
        raise error
    monkeypatch.setattr(module.subprocess, "run", run)
    result = module.audit(source, module.record(source)["sha256"], tmp_path / "out")
    row = result["rows"][0]
    assert row["status"] == "unresolved" and row["error_type"] == type(error).__name__
    assert commands == [["/usr/bin/sum", str(compressed)]] and len(calls) == 1
    if isinstance(error, subprocess.CalledProcessError):
        assert (row["command_returncode"], row["command_stdout"], row["command_stderr"]) == (1, "partial", "failure")
    assert json.loads((tmp_path / "out/report.json").read_text()) == result


def test_postflight_retained_mutation_still_fails(tmp_path, monkeypatch, offline_source):
    source, compressed, _ = offline_source
    def run(command, **kwargs):
        compressed.write_bytes(b"xyz")
        return subprocess.CompletedProcess(command, 0, f"16556 1 {compressed}\n", "")
    monkeypatch.setattr(module.subprocess, "run", run)
    with pytest.raises(ValueError):
        module.audit(source, module.record(source)["sha256"], tmp_path / "out")
    assert not (tmp_path / "out/report.json").exists()


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
