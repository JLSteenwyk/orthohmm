import zipfile

import pytest

from benchmark_tools.audit_recovery_install import installed_payload


def test_exact_package_payload_and_explicit_exclusions(tmp_path):
    wheel = tmp_path / "example.whl"
    site = tmp_path / "site"
    site.mkdir()
    (site / "module.py").write_bytes(b"source")
    with zipfile.ZipFile(wheel, "w") as z:
        z.writestr("module.py", b"source")
        z.writestr("pkg-1.dist-info/RECORD", b"record")
        z.writestr("pkg-1.data/scripts/tool", b"script")
    result = installed_payload(wheel, site)
    assert result["matched_files"] == 1
    assert len(result["excluded"]) == 2


def test_vendored_record_is_payload_not_installer_generated(tmp_path):
    site = tmp_path / "site"
    path = site / "pkg/vendor/example.dist-info/RECORD"
    path.parent.mkdir(parents=True)
    path.write_bytes(b"changed")
    wheel = tmp_path / "vendor.whl"
    with zipfile.ZipFile(wheel, "w") as z:
        z.writestr("pkg/vendor/example.dist-info/RECORD", b"original")
    with pytest.raises(ValueError, match="payload differs"):
        installed_payload(wheel, site)


@pytest.mark.parametrize("change", ["content", "missing", "symlink", "traversal"])
def test_changed_or_unsafe_package_rejected(tmp_path, change):
    site = tmp_path / "site"
    site.mkdir()
    target = site / "module.py"
    if change == "content":
        target.write_bytes(b"other")
    elif change == "symlink":
        outside = tmp_path / "outside"
        outside.write_bytes(b"source")
        target.symlink_to(outside)
    wheel = tmp_path / "example.whl"
    with zipfile.ZipFile(wheel, "w") as z:
        z.writestr("../outside" if change == "traversal" else "module.py", b"source")
    with pytest.raises((ValueError, FileNotFoundError)):
        installed_payload(wheel, site)


@pytest.mark.parametrize("source_present", [True, False])
def test_regenerated_bytecode_requires_corresponding_verified_source(tmp_path, source_present):
    site = tmp_path / "site"
    (site / "pkg").mkdir(parents=True)
    (site / "pkg/module.py").write_bytes(b"source")
    wheel = tmp_path / "bytecode.whl"
    with zipfile.ZipFile(wheel, "w") as z:
        z.writestr("pkg/__pycache__/module.cpython-310.pyc", b"regenerated")
        if source_present:
            z.writestr("pkg/module.py", b"source")
    if source_present:
        result = installed_payload(wheel, site)
        assert result["matched_files"] == 1 and len(result["excluded"]) == 1
        (site / "pkg/module.py").write_bytes(b"altered")
        with pytest.raises(ValueError, match="payload differs"):
            installed_payload(wheel, site)
    else:
        with pytest.raises(ValueError, match="lacks auditable source"):
            installed_payload(wheel, site)
