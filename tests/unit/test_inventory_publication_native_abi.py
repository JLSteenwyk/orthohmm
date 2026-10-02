import json
from pathlib import Path

import pytest

from benchmark_tools import inventory_publication_native_abi as mod


def fixture(tmp_path, monkeypatch):
    root = tmp_path / "assembly"
    root.mkdir()
    files = [("assets/a.whl", b"wheel"), ("reader_wheels/a.whl", b"wheel"),
             ("assets/FastTree", b"\x7fELFnative"), ("assets/launcher", b"shell")]
    entries = []
    for relative, data in files:
        path = root / relative
        path.parent.mkdir(exist_ok=True, parents=True)
        path.write_bytes(data)
        entries.append(dict(kind="file", path=relative))
    entries.append(dict(kind="symlink", path="assets/alias"))
    (root / mod.assembly.INDEX).write_text(json.dumps(dict(entries=entries)))
    inspector = tmp_path / "readelf"
    inspector.write_bytes(b"inspector")
    checks = []
    def validate(*args):
        checks.append(args)
        return dict(status="verified", manifest=mod.record(root / mod.assembly.INDEX))
    monkeypatch.setattr(mod.assembly, "validate", validate)
    requirements = dict(header={"Machine": "test-machine"}, interpreter="/loader", version_requirements=[
        dict(library="libc.so.6", versions=[dict(name="GLIBC_2.34", flags="none", index=2)])])
    monkeypatch.setattr(mod, "scan_abi", lambda *args: dict(requirements=requirements, readelf_stdout="raw"))
    monkeypatch.setattr(mod, "scan_wheel", lambda w, r, **k: dict(objects=[dict(abi=dict(requirements=requirements))], wheel=w["wheel"]))
    return root, inspector, checks


def test_assembly_scan_deduplicates_exact_wheels_only(tmp_path, monkeypatch):
    root, inspector, checks = fixture(tmp_path, monkeypatch)
    result = mod.inspect(root, "anchor", inspector)
    assert len(checks) == 2
    assert result["unique_wheels"] == result["wheel_elf_objects"] == result["native_tool_objects"] == 1
    assert result["machines"] == ["test-machine"]
    assert result["required_library_versions"] == [dict(library="libc.so.6", name="GLIBC_2.34", flags="none")]
    assert all(result[key] is False for key in ("native_code_executed", "runtime_resolution_verified",
        "cross_host_compatibility_verified", "controlled_timing", "publication_ready", "security_clearance", "redistribution_clearance"))


def test_bad_anchor_precedes_inspection(tmp_path, monkeypatch):
    root, inspector, checks = fixture(tmp_path, monkeypatch)
    def reject(*args):
        raise ValueError("anchor")
    monkeypatch.setattr(mod.assembly, "validate", reject)
    monkeypatch.setattr(mod, "scan_wheel", lambda *a, **k: pytest.fail("must not inspect"))
    with pytest.raises(ValueError, match="anchor"):
        mod.inspect(root, "wrong", inspector)


def test_changed_inspector_rejected(tmp_path, monkeypatch):
    root, inspector, checks = fixture(tmp_path, monkeypatch)
    def mutate(*args):
        inspector.write_bytes(b"changed")
        return dict(requirements={})
    monkeypatch.setattr(mod, "scan_abi", mutate)
    with pytest.raises(ValueError):
        mod.inspect(root, "anchor", inspector)


def test_changed_assembly_rejected(tmp_path, monkeypatch):
    root, inspector, checks = fixture(tmp_path, monkeypatch)
    def mutate(*args):
        (root / mod.assembly.INDEX).write_text("changed")
        return dict(requirements={})
    monkeypatch.setattr(mod, "scan_abi", mutate)
    with pytest.raises(ValueError, match="Assembly changed"):
        mod.inspect(root, "anchor", inspector)
