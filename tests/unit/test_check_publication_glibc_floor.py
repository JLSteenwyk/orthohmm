import copy
import json

import pytest

from benchmark_tools import check_publication_glibc_floor as mod
from benchmark_tools import assemble_publication_runtime_assets as assembly


def inventory():
    return dict(schema="publication_native_abi_inventory_v1", status="declared_abi_inventory_completed",
        machines=["Advanced Micro Devices X86-64"], required_library_versions=[
            dict(library="libc.so.6", name="GLIBC_2.34", flags="none"),
            dict(library="libm.so.6", name="GLIBC_2.2.5", flags="none"),
            dict(library="libgomp.so.1", name="GOMP_4.5", flags="none")])


@pytest.mark.parametrize("host", ["glibc 2.34", "glibc 2.34.0", "glibc 2.39", "glibc 2.100"])
def test_numeric_floor_is_not_full_compatibility(host):
    result = mod.evaluate(inventory(), "Linux", "x86_64", host)
    assert result["declared_glibc_floor"] == "2.34"
    assert result["glibc_floor_check_passed"] is True
    assert all(result[k] is False for k in ("full_host_compatibility_verified", "runtime_resolution_verified",
        "reviewed_native_code_executed", "controlled_timing", "publication_ready"))


@pytest.mark.parametrize("host", ["glibc 2.9", "glibc 2.33", "glibc 2.3.4", "musl 1.2.5", None, "", "glibc unknown"])
def test_old_or_unknown_glibc_fails(host):
    with pytest.raises(ValueError):
        mod.evaluate(inventory(), "Linux", "x86_64", host)


@pytest.mark.parametrize("system,machine", [("Darwin", "x86_64"), ("Linux", "aarch64"), ("Windows", "AMD64")])
def test_unsupported_platform_fails(system, machine):
    with pytest.raises(ValueError):
        mod.evaluate(inventory(), system, machine, "glibc 2.39")


@pytest.mark.parametrize("requirement", ["GLIBC_PRIVATE", "GLIBC_ABI_DT_RELR", "GLIBC_2", "GLIBC_2.34.extra"])
def test_non_numeric_glibc_cannot_be_claimed_satisfied(requirement):
    item = inventory()
    item["required_library_versions"][0]["name"] = requirement
    with pytest.raises(ValueError):
        mod.evaluate(item, "Linux", "x86_64", "glibc 2.39")


def test_weak_glibc_requirements_are_conservatively_included():
    item = inventory()
    item["required_library_versions"][0].update(name="GLIBC_2.40", flags="WEAK")
    with pytest.raises(ValueError, match="below"):
        mod.evaluate(item, "Linux", "x86_64", "glibc 2.39")


@pytest.mark.parametrize("field,value", [("schema", "other"), ("status", "pending"), ("machines", ["AArch64"]),
    ("required_library_versions", [dict(library="libgomp", name="GOMP_4.5")])])
def test_incomplete_or_wrong_inventory_fails(field, value):
    item = inventory()
    item[field] = value
    with pytest.raises(ValueError):
        mod.evaluate(item, "Linux", "x86_64", "glibc 2.39")


def setup(tmp_path, monkeypatch):
    root = tmp_path / "bundle"
    root.mkdir()
    loader = tmp_path / "loader"
    loader.write_bytes(b"loader fixture")
    loader.chmod(0o755)
    expected = dict(status="verified", files=102, symlinks=18, payload_bytes=132139954,
        scientific_revision="scientific", manifest=dict(path="historical/index", sha256="a" * 64))
    item = inventory()
    item.update(assembly=copy.deepcopy(expected), interpreters=[str(loader)])
    path = tmp_path / "abi.json"
    path.write_text(json.dumps(item))
    monkeypatch.setattr(assembly, "validate", lambda *a: dict(expected, manifest=dict(path=str(root / "index"), sha256="a" * 64)))
    monkeypatch.setattr(mod.os, "confstr", lambda *a: "glibc 2.39")
    monkeypatch.setattr(mod.platform, "system", lambda: "Linux")
    monkeypatch.setattr(mod.platform, "machine", lambda: "x86_64")
    for key in ("LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"):
        monkeypatch.delenv(key, raising=False)
    return root, loader, path, item


def test_actual_file_binding_allows_relocated_assembly(tmp_path, monkeypatch):
    root, loader, path, item = setup(tmp_path, monkeypatch)
    result = mod.check_host(root, path, mod.record(path)["sha256"], "a" * 64)
    assert result["inventory"] == mod.record(path)
    assert result["loaders"] == [dict(declared_path=str(loader), identity=mod.record(loader))]
    assert result["assembly"]["manifest"]["path"] == str(root / "index")


@pytest.mark.parametrize("field", ["sha256", "payload_bytes", "files", "scientific_revision"])
def test_wrong_assembly_binding_fails(tmp_path, monkeypatch, field):
    root, loader, path, item = setup(tmp_path, monkeypatch)
    if field == "sha256":
        item["assembly"]["manifest"][field] = "b" * 64
    else:
        item["assembly"][field] = "wrong"
    path.write_text(json.dumps(item))
    with pytest.raises(ValueError, match="different assembly"):
        mod.check_host(root, path, mod.record(path)["sha256"], "a" * 64)


def test_bad_digest_precedes_host_observation(tmp_path, monkeypatch):
    root, loader, path, item = setup(tmp_path, monkeypatch)
    monkeypatch.setattr(mod.os, "confstr", lambda *a: pytest.fail("must not observe"))
    with pytest.raises(ValueError, match="external anchor"):
        mod.check_host(root, path, "b" * 64, "a" * 64)


@pytest.mark.parametrize("key", ["LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"])
def test_loader_overrides_rejected(tmp_path, monkeypatch, key):
    root, loader, path, item = setup(tmp_path, monkeypatch)
    monkeypatch.setenv(key, "override")
    with pytest.raises(ValueError, match="overrides"):
        mod.check_host(root, path, mod.record(path)["sha256"], "a" * 64)


def test_unavailable_confstr_fails(tmp_path, monkeypatch):
    root, loader, path, item = setup(tmp_path, monkeypatch)
    def missing(*a):
        raise ValueError("unsupported")
    monkeypatch.setattr(mod.os, "confstr", missing)
    with pytest.raises(ValueError, match="unavailable"):
        mod.check_host(root, path, mod.record(path)["sha256"], "a" * 64)


def test_missing_loader_fails(tmp_path, monkeypatch):
    root, loader, path, item = setup(tmp_path, monkeypatch)
    loader.unlink()
    with pytest.raises(ValueError, match="interpreter is unavailable"):
        mod.check_host(root, path, mod.record(path)["sha256"], "a" * 64)


def test_inventory_changed_during_observation_fails(tmp_path, monkeypatch):
    root, loader, path, item = setup(tmp_path, monkeypatch)
    def mutate(*a):
        path.write_text("changed")
        return "glibc 2.39"
    monkeypatch.setattr(mod.os, "confstr", mutate)
    with pytest.raises(ValueError, match="changed"):
        mod.check_host(root, path, mod.record(path)["sha256"], "a" * 64)
