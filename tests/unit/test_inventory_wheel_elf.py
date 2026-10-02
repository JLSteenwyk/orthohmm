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


ABI_HEADER = """ELF Header:
  Class: ELF64
  Data: 2's complement, little endian
  OS/ABI: UNIX - System V
  ABI Version: 0
  Type: DYN (Shared object file)
  Machine: Advanced Micro Devices X86-64
  Flags: 0x0
  Number of program headers: 1

Program Headers:
  INTERP 0x01 0x02
  [Requesting program interpreter: /lib64/ld-linux-x86-64.so.2]
"""
ABI_NEEDS = """
Version symbols section '.gnu.version' contains 1 entry:
  000: 2 (GLIBC_2.34)
Version needs section '.gnu.version_r' contains 1 entry:
 Addr: 0x01  Offset: 0x02  Link: 3 (.dynstr)
 000000: Version: 1  File: libc.so.6  Cnt: 2
 0x0010: Name: GLIBC_2.34  Flags: none  Version: 2
 0x0020: Name: GLIBC_2.2.5  Flags: WEAK  Version: 3
"""


def test_abi_requirements_keep_library_and_weak_flags():
    result = mod.abi_requirements(ABI_HEADER + ABI_NEEDS)
    assert result["interpreter"] == "/lib64/ld-linux-x86-64.so.2"
    assert result["header"]["Machine"] == "Advanced Micro Devices X86-64"
    assert result["version_requirements"] == [dict(library="libc.so.6", versions=[
        dict(name="GLIBC_2.34", flags="none", index=2), dict(name="GLIBC_2.2.5", flags="WEAK", index=3)])]


def test_static_abi_is_not_a_compatibility_pass():
    text = ABI_HEADER.replace("1\n\nProgram Headers:\n  INTERP 0x01 0x02\n  [Requesting program interpreter: /lib64/ld-linux-x86-64.so.2]", "0\n\nThere are no program headers in this file.")
    result = mod.abi_requirements(text + "No version information found in this file.\n")
    assert result["interpreter"] is None
    assert result["version_requirements"] == []


def test_version_definitions_are_not_requirements():
    text = ABI_HEADER + "Version definition section '.gnu.version_d' contains 1 entry:\n  Name: GLIBC_99.0\n"
    assert mod.abi_requirements(text)["version_requirements"] == []


@pytest.mark.parametrize("old,new", [
    ("  Class: ELF64\n", ""),
    ("  Class: ELF64\n", "  Class: ELF64\n  Class: ELF32\n"),
    ("Program Headers:", ""),
    ("Number of program headers: 1", "Number of program headers: 2"),
    ("  [Requesting program interpreter: /lib64/ld-linux-x86-64.so.2]\n", ""),
    ("  INTERP 0x01 0x02\n", ""),
    ("contains 1 entry:", "contains 2 entries:"),
    ("Cnt: 2", "Cnt: 3"),
    ("Cnt: 2", "Cnt: 1"),
    ("File: libc.so.6", "File: libc.so.6 unexpected"),
    ("Name: GLIBC_2.34", "Broken: GLIBC_2.34"),
    ("GLIBC_2.2.5", "GLIBC_2.34"),
    ("Version: 1  File:", "Version: 2  File:"),
])
def test_abi_fail_closed(old, new):
    with pytest.raises(ValueError):
        mod.abi_requirements((ABI_HEADER + ABI_NEEDS).replace(old, new))


def test_abi_missing_version_evidence():
    with pytest.raises(ValueError, match="version-information"):
        mod.abi_requirements(ABI_HEADER)


def test_optional_abi_scan_preserves_default_contract(tmp_path, monkeypatch):
    item = wheel(tmp_path, [("native", b"\x7fELFfake")])
    calls = []
    def run(command, **kwargs):
        calls.append(command)
        assert kwargs["env"]["LC_ALL"] == "C"
        assert kwargs["check"] is True and kwargs["timeout"] == 60
        output = ABI_HEADER + ABI_NEEDS if "--version-info" in command else " 0x1 (NEEDED) Shared library: [libc.so.6]\n"
        return subprocess.CompletedProcess(command, 0, output, "")
    monkeypatch.setattr(mod.subprocess, "run", run)
    default = mod.scan_wheel(item, Path("/usr/bin/readelf"))
    extended = mod.scan_wheel(item, Path("/usr/bin/readelf"), include_abi=True)
    abi = extended["objects"][0].pop("abi")
    assert default == extended
    assert abi["requirements"]["version_requirements"]
    assert len(calls) == 3


def test_abi_diagnostics_fail(tmp_path, monkeypatch):
    monkeypatch.setattr(mod.subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(a, 0, ABI_HEADER + ABI_NEEDS, "warning"))
    with pytest.raises(ValueError, match="ABI diagnostics"):
        mod.scan_abi(tmp_path / "binary", Path("/usr/bin/readelf"))
