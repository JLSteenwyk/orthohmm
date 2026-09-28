import copy
import io
import json
import tarfile

import pytest

from benchmark_tools import acquire_wheel_sources as mod


def release():
    wheel = dict(name="example_pkg", version="1.0", wheel=dict(path="/tmp/example.whl", bytes=12, sha256="a" * 64))
    source = dict(filename="example-1.0.tar.gz", packagetype="sdist", size=99,
                  url="https://files.pythonhosted.org/packages/example.tar.gz", digests=dict(sha256="b" * 64))
    meta = dict(info=dict(name="example-pkg", version="1.0"), urls=[
        dict(filename="example.whl", size=12, packagetype="bdist_wheel", digests=dict(sha256="a" * 64)), source])
    return wheel, meta


def test_select_exact_release():
    wheel, meta = release()
    assert mod.select_source(meta, wheel) == meta["urls"][1]


@pytest.mark.parametrize("target,key,value", [
    ("info", "name", "other"), ("info", "version", "2"),
    ("wheel", "filename", "other.whl"), ("wheel", "size", 11),
    ("wheel", "digests", {"sha256": "c" * 64}), ("wheel", "packagetype", "sdist"),
    ("source", "url", "http://files.pythonhosted.org/file"),
    ("source", "url", "https://files.pythonhosted.org.evil.test/file"),
    ("source", "url", "https://files.pythonhosted.org/file?x=1"),
    ("source", "filename", "../bad.tar.gz"), ("source", "filename", "dir/file.tar.gz"),
    ("source", "filename", "example.zip"), ("source", "size", True),
    ("source", "size", 100_000_001), ("source", "digests", {"sha256": "bad"}),
])
def test_select_rejects_mismatches(target, key, value):
    wheel, meta = release()
    selected = meta["info"] if target == "info" else meta["urls"][target == "source"]
    selected[key] = value
    with pytest.raises(ValueError):
        mod.select_source(meta, wheel)


@pytest.mark.parametrize("index", [0, 1])
def test_select_rejects_ambiguity(index):
    wheel, meta = release()
    meta["urls"].append(copy.deepcopy(meta["urls"][index]))
    with pytest.raises(ValueError):
        mod.select_source(meta, wheel)


def source_archive(tmp_path, extra=(), package=b"Name: example-pkg\nVersion: 1.0\n"):
    path = tmp_path / "source.tar.gz"
    with tarfile.open(path, "w:gz") as archive:
        for name, data, kind in [("root/PKG-INFO", package, tarfile.REGTYPE),
                                  ("root/LICENSE", b"notice", tarfile.REGTYPE), *extra]:
            item = tarfile.TarInfo(name)
            item.type = kind
            item.size = len(data)
            archive.addfile(item, io.BytesIO(data))
    return path


def test_inspection_no_extraction(tmp_path):
    path = source_archive(tmp_path, [("root/vendor/core/COPYING", b"core notice", tarfile.REGTYPE)])
    result = mod.inspect_source(path, "example_pkg", "1.0")
    assert len(result["files"]) == 3
    assert len(result["notice_candidates"]) == 2
    assert result["package_name"] == "example-pkg"
    assert list(tmp_path.iterdir()) == [path]


@pytest.mark.parametrize("name,kind", [
    ("../bad", tarfile.REGTYPE), ("/absolute", tarfile.REGTYPE),
    ("root/LICENSE", tarfile.REGTYPE), ("elsewhere/file", tarfile.REGTYPE),
    ("root/link", tarfile.SYMTYPE), ("root/link", tarfile.LNKTYPE),
    ("root/fifo", tarfile.FIFOTYPE),
])
def test_inspection_rejects_unsafe_members(tmp_path, name, kind):
    with pytest.raises(ValueError):
        mod.inspect_source(source_archive(tmp_path, [(name, b"", kind)]), "example-pkg", "1.0")


@pytest.mark.parametrize("package", [b"Name: other\nVersion: 1.0\n", b"Name: example-pkg\nVersion: 2\n"])
def test_source_metadata_mismatch(tmp_path, package):
    with pytest.raises(ValueError, match="identity"):
        mod.inspect_source(source_archive(tmp_path, package=package), "example-pkg", "1.0")


def test_download_limit_and_no_overwrite(tmp_path, monkeypatch):
    class Response(io.BytesIO):
        def geturl(self):
            return "https://test"
    monkeypatch.setattr(mod, "urlopen", lambda *a, **k: Response(b"abc"))
    path = tmp_path / "download"
    assert mod.download("https://test", path, 3)["bytes"] == 3
    with pytest.raises(FileExistsError):
        mod.download("https://test", path, 3)
    with pytest.raises(ValueError, match="limit"):
        mod.download("https://test", tmp_path / "limited", 2)
    with pytest.raises(ValueError, match="redirect"):
        mod.download("https://other", tmp_path / "redirect", 3)


@pytest.mark.parametrize("corrupt", [False, True])
def test_acquire_binds_artifacts_without_build(tmp_path, monkeypatch, corrupt):
    wheel, meta = release()
    wheel_path = tmp_path / "example.whl"
    wheel_path.write_bytes(b"wheel contents")
    wheel["wheel"] = mod.record(wheel_path)
    meta["urls"][0].update(size=wheel["wheel"]["bytes"], digests=dict(sha256=wheel["wheel"]["sha256"]))
    archive = source_archive(tmp_path)
    identity = mod.record(archive)
    meta["urls"][1].update(size=identity["bytes"], digests=dict(sha256=identity["sha256"]))
    inventory = tmp_path / "inventory.json"
    inventory.write_text(json.dumps(dict(status="selected_wheel_elf_inventory", wheels=[wheel])))
    calls = []
    def fake_download(url, path, maximum):
        calls.append(url)
        path.write_bytes(json.dumps(meta).encode() if url.endswith("/json")
                         else (b"bad" if corrupt else archive.read_bytes()))
        return mod.record(path)
    monkeypatch.setattr(mod, "download", fake_download)
    if corrupt:
        with pytest.raises(ValueError, match="download identity"):
            mod.acquire(inventory, ["example_pkg"], tmp_path / "output")
    else:
        result = mod.acquire(inventory, ["example_pkg"], tmp_path / "output")
        assert result["packages"][0]["inspection"]["archive"]["sha256"] == identity["sha256"]
        assert result["publication_ready"] is False
        assert result["redistribution_clearance"] is False
        assert result["packages"][0]["wheel"] == wheel["wheel"]
    assert len(calls) == 2
