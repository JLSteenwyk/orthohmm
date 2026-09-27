import json
from pathlib import Path
from urllib.parse import urlsplit, unquote

import pytest

from benchmark_tools.stage_base_archives import identity, packages_from_receipt, stage


@pytest.fixture
def inputs(tmp_path):
    cache = tmp_path / "cache"
    cache.mkdir()
    archive = cache / "example-1.0-main.conda"
    archive.write_bytes(b"retained archive bytes")
    hashes = identity(archive)
    package = dict(name="example", fn=archive.name, version="1.0", build="main",
        subdir="linux-64", url="https://repo.anaconda.com/pkgs/main/linux-64/" + archive.name,
        sha256=hashes["sha256"], md5=hashes["md5"])
    data = dict(acquisition=dict(packages=[dict(package=package,
        archive=dict(path="/unavailable/old/archive", bytes=hashes["bytes"],
                     sha256=hashes["sha256"]))]))
    receipt = tmp_path / "receipt.json"
    receipt.write_text(json.dumps(data))
    return receipt, cache, data


def test_staging_and_relocation(inputs, tmp_path):
    receipt, cache, _ = inputs
    first = tmp_path / "first"
    result = stage(receipt, cache, first)
    assert result["status"] == "base_archives_staged"
    assert result["installation_performed"] is False
    second = tmp_path / "second with spaces"
    other = stage(receipt, first / "archives", second)
    assert result["packages"] == other["packages"]
    lines = (second / "explicit.txt").read_text().splitlines()
    assert lines[0] == "@EXPLICIT"
    url = urlsplit(lines[1])
    assert Path(unquote(url.path)).parent == second / "archives"
    assert url.fragment == result["packages"][0]["archive"]["md5"]
    assert json.loads((second / "staging.json").read_text()) == other
    with pytest.raises(FileExistsError):
        stage(receipt, cache, first)


@pytest.mark.parametrize("field,value", [
    ("fn", "../escape.conda"), ("fn", "bad\n.conda"),
    ("name", "../escape"), ("sha256", "0" * 64), ("md5", "nope"),
    ("subdir", "osx-arm64"), ("url", "http://repo.anaconda.com/file.conda"),
    ("url", "https://unknown.example/example-1.0-main.conda"),
    ("url", "https://user:password@repo.anaconda.com/example-1.0-main.conda"),
    ("url", "https://repo.anaconda.com/wrong.conda"),
])
def test_bad_metadata(inputs, field, value):
    _, _, data = inputs
    data["acquisition"]["packages"][0]["package"][field] = value
    with pytest.raises(ValueError):
        packages_from_receipt(data)


def test_duplicate_and_empty(inputs):
    _, _, data = inputs
    data["acquisition"]["packages"] *= 2
    with pytest.raises(ValueError, match="Duplicate"):
        packages_from_receipt(data)
    data["acquisition"]["packages"] = []
    with pytest.raises(ValueError, match="Empty"):
        packages_from_receipt(data)


@pytest.mark.parametrize("mode", ["changed", "missing", "symlink"])
def test_invalid_cache_leaves_no_output(inputs, tmp_path, mode):
    receipt, cache, _ = inputs
    archive = next(cache.iterdir())
    if mode == "changed":
        archive.write_bytes(b"corruption")
    else:
        archive.unlink()
        if mode == "symlink":
            archive.symlink_to(receipt)
    output = tmp_path / "out"
    with pytest.raises(ValueError):
        stage(receipt, cache, output)
    assert not output.exists()


def test_digest_mismatch_during_copy(inputs, tmp_path, monkeypatch):
    receipt, cache, _ = inputs
    import benchmark_tools.stage_base_archives as module
    monkeypatch.setattr(module.shutil, "copyfileobj", lambda src, dst, length: dst.write(b"bad"))
    output = tmp_path / "out"
    with pytest.raises(ValueError, match="Copied archive differs"):
        stage(receipt, cache, output)
    assert output.exists()
    assert not (output / "staging.json").exists()
    assert not (output / "explicit.txt").exists()


def test_receipt_changes_during_copy(inputs, tmp_path, monkeypatch):
    receipt, cache, _ = inputs
    import benchmark_tools.stage_base_archives as module
    original = module.shutil.copyfileobj

    def change_receipt(src, dst, length):
        original(src, dst, length)
        receipt.write_text("{}")

    monkeypatch.setattr(module.shutil, "copyfileobj", change_receipt)
    output = tmp_path / "out"
    with pytest.raises(ValueError, match="Receipt changed"):
        stage(receipt, cache, output)
    assert not (output / "staging.json").exists()


def test_broken_output_symlink(inputs, tmp_path):
    receipt, cache, _ = inputs
    output = tmp_path / "out"
    output.symlink_to(tmp_path / "missing")
    with pytest.raises(FileExistsError):
        stage(receipt, cache, output)


@pytest.mark.parametrize("size", [True, 0, -1, "20"])
def test_invalid_size(inputs, size):
    _, _, data = inputs
    data["acquisition"]["packages"][0]["archive"]["bytes"] = size
    with pytest.raises(ValueError, match="Inconsistent archive"):
        packages_from_receipt(data)
