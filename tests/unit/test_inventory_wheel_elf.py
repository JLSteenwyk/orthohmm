import subprocess
import stat
from pathlib import Path
import zipfile

import pytest

from benchmark_tools import inventory_wheel_elf as mod


def test_dynamic_tags():
    text = """Dynamic section at offset 0x123 contains 3 entries:
 0x0000000000000001 (NEEDED) Shared library: [libc.so.6]
 0x000000000000001d (RUNPATH) Library runpath: [$ORIGIN/../libs]
 0x000000000000000e (SONAME) Library soname: [example.so]
"""
    assert mod.dynamic_tags(text) == dict(NEEDED=["libc.so.6"], RUNPATH=["$ORIGIN/../libs"], SONAME=["example.so"], RPATH=[])


@pytest.mark.parametrize("text", ["(NEEDED) bad", " 0x1 (SONAME) name: [a]\n 0x2 (SONAME) name: [b]"])
def test_bad_dynamic_tags(text):
    with pytest.raises(ValueError):
        mod.dynamic_tags(text)


def test_static_object():
    assert all(not v for v in mod.dynamic_tags("There is no dynamic section in this file.").values())


def wheel(tmp_path, members):
    path = tmp_path / "test.whl"
    with zipfile.ZipFile(path, "w") as archive:
        for name, data in members:
            archive.writestr(name, data)
    return dict(name="test", version="1", wheel=mod.record(path), notice_candidates=[], native_members=["fake.so"])


def test_magic_scan_not_suffix(tmp_path, monkeypatch):
    item = wheel(tmp_path, [("native_without_suffix", b"\x7fELFfake"), ("fake.so", b"plain")])
    calls = []
    def run(command, **kwargs):
        calls.append(command)
        assert Path(command[-1]).read_bytes() == b"\x7fELFfake"
        return subprocess.CompletedProcess(command, 0, " 0x1 (NEEDED) Shared library: [libc.so.6]\n", "")
    monkeypatch.setattr(mod.subprocess, "run", run)
    result = mod.scan_wheel(item, Path("/usr/bin/readelf"))
    assert len(calls) == 1
    assert result["objects"][0]["member"] == "native_without_suffix"
    assert result["filename_candidates_not_elf"] == ["fake.so"]
    assert not Path(calls[0][-1]).exists()


def test_reject_diagnostics(tmp_path, monkeypatch):
    item = wheel(tmp_path, [("a", b"\x7fELFfake")])
    monkeypatch.setattr(mod.subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(a, 0, "", "warning"))
    with pytest.raises(ValueError, match="diagnostics"):
        mod.scan_wheel(item, Path("/usr/bin/readelf"))


@pytest.mark.parametrize("name", ["../escape", "/absolute", "a\\b"])
def test_reject_paths(tmp_path, name):
    with pytest.raises(ValueError, match="Unsafe"):
        mod.scan_wheel(wheel(tmp_path, [(name, b"x")]), Path("/usr/bin/readelf"))


def test_candidates_are_not_resolution():
    tags = dict(NEEDED=["alias.so", "absent.so"], SONAME=["alias.so"], RPATH=[], RUNPATH=[])
    wheels = [dict(wheel={"sha256": "a"}, objects=[dict(member="lib/object", tags=tags)])]
    edges = mod.candidate_edges(wheels)
    assert edges[0]["selected_wheel_candidates"] == [dict(wheel_sha256="a", member="lib/object")]
    assert edges[1]["selected_wheel_candidates"] == []


def test_reject_changed_wheel(tmp_path):
    item = wheel(tmp_path, [("a", b"x")])
    Path(item["wheel"]["path"]).write_bytes(b"changed")
    with pytest.raises(ValueError):
        mod.scan_wheel(item, Path("/usr/bin/readelf"))


def test_reject_duplicate_members(tmp_path):
    with pytest.warns(UserWarning, match="Duplicate"):
        item = wheel(tmp_path, [("a", b"x"), ("a", b"y")])
    with pytest.raises(ValueError, match="Duplicate"):
        mod.scan_wheel(item, Path("/usr/bin/readelf"))


def test_reject_symlink(tmp_path):
    info = zipfile.ZipInfo("link")
    info.external_attr = (stat.S_IFLNK | 0o777) << 16
    with pytest.raises(ValueError, match="Nonregular"):
        mod.scan_wheel(wheel(tmp_path, [(info, b"target")]), Path("/usr/bin/readelf"))
