import hashlib
import io
import json

import pytest

from benchmark_tools import acquire_publication_fasttree as module


class Response(io.BytesIO):
    def geturl(self):
        return "https://example.test/file"


def test_download_pins_and_preserves_mismatch(tmp_path, monkeypatch):
    monkeypatch.setattr(module, "urlopen", lambda *a, **kw: Response(b"changed"))
    target = tmp_path / "download"
    with pytest.raises(ValueError, match="SHA-256"):
        module.download("https://example.test/file", target, "0" * 64)
    assert target.read_bytes() == b"changed"


def test_download_rejects_oversize(tmp_path, monkeypatch):
    monkeypatch.setattr(module, "MAX_BYTES", 3)
    monkeypatch.setattr(module, "urlopen", lambda *a, **kw: Response(b"large"))
    with pytest.raises(ValueError, match="size limit"):
        module.download("https://example.test/file", tmp_path / "file", "0" * 64)


def test_download_rejects_insecure_redirect(tmp_path, monkeypatch):
    monkeypatch.setattr(Response, "geturl", lambda self: "http://example.test/file")
    monkeypatch.setattr(module, "urlopen", lambda *a, **kw: Response(b"file"))
    with pytest.raises(ValueError, match="Non-HTTPS"):
        module.download("https://example.test/file", tmp_path / "file", "0" * 64)


def test_symlink_rejected(tmp_path):
    original = tmp_path / "file"
    original.write_bytes(b"test")
    link = tmp_path / "link"
    link.symlink_to(original)
    with pytest.raises(ValueError, match="nonsymlink"):
        module.identity(link)


@pytest.fixture
def acquisition(tmp_path, monkeypatch):
    binary = tmp_path / "installed"
    binary.write_bytes(b"binary")
    pins = {"FastTree": hashlib.sha256(b"binary").hexdigest()}
    monkeypatch.setattr(module, "FILES", pins)
    monkeypatch.setattr(module, "urlopen", lambda *a, **kw: Response(b"binary"))
    return tmp_path / "output", binary


def test_success_and_existing_directory_refusal(acquisition):
    output, binary = acquisition
    report = module.run(output, binary)
    assert report["installed_matches_upstream_binary"]
    assert not report["reproducible_build_verified"]
    assert not report["executed_downloads"]
    assert json.loads((output / "receipt.json").read_text()) == report
    with pytest.raises(FileExistsError):
        module.run(output, binary)


def test_failed_download_receipt(acquisition, monkeypatch):
    output, binary = acquisition
    monkeypatch.setattr(module, "urlopen", lambda *a, **kw: Response(b"bad"))
    with pytest.raises(ValueError):
        module.run(output, binary)
    report = json.loads((output / "receipt.json").read_text())
    assert report["status"] == "acquisition_failed"
    assert "installed_matches_upstream_binary" not in report


def test_rejects_wrong_installed_binary(acquisition):
    output, binary = acquisition
    binary.write_bytes(b"other")
    with pytest.raises(ValueError, match="SHA-256"):
        module.run(output, binary)
    assert not output.exists()
