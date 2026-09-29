import hashlib
import zipfile

import pytest

from benchmark_tools.build_private_timing_environment import (
    PACKAGES, build, hash_lock, selected_versions, wheel_metadata,
)


def wheel(path, name="Example", version="1", extra=False):
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("example.dist-info/METADATA", f"Name: {name}\nVersion: {version}\n")
        if extra:
            archive.writestr("other.dist-info/METADATA", "Name: Other\nVersion: 1\n")


def test_selected_scientific_versions_preserved():
    packages = {name: str(i) for i, name in enumerate(PACKAGES)}
    packages.update(pip="old", unrelated="9")
    result = selected_versions({"environments": {"orthohmm": {"packages": packages}}})
    assert {k: v for k, v in result.items() if k != "pip"} == {k: packages[k] for k in PACKAGES}
    assert result["pip"] == "26.2.1" and "unrelated" not in result


def test_missing_baseline_version_rejected():
    with pytest.raises(KeyError):
        selected_versions({"environments": {"orthohmm": {"packages": {}}}})


def test_patched_deployment_only_changes_declared_packages():
    packages = {name: "1" for name in PACKAGES}
    baseline = {"environments": {"orthohmm": {"packages": packages}}}
    original = selected_versions(baseline)
    patched = selected_versions(baseline, patched_deployment=True)
    assert {k: v for k, v in patched.items() if original[k] != v} == {
        "packaging": "26.1", "pillow": "12.3.0", "setuptools": "83.0.0"}
    assert all(v == "1" for v in packages.values())


def test_lock_binds_actual_wheel(tmp_path):
    path = tmp_path / "example.whl"
    wheel(path)
    text, rows = hash_lock(tmp_path, {"example": "1"})
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    assert text == f"Example==1 --hash=sha256:{digest}\n"
    assert rows[0]["sha256"] == digest


@pytest.mark.parametrize("selected", [{}, {"Example": "2"}, {"Example": "1", "Missing": "1"}])
def test_incomplete_or_unexpected_wheels_rejected(tmp_path, selected):
    wheel(tmp_path / "example.whl")
    with pytest.raises(ValueError):
        hash_lock(tmp_path, selected)


def test_duplicate_wheels_rejected(tmp_path):
    wheel(tmp_path / "one.whl")
    wheel(tmp_path / "two.whl")
    with pytest.raises(ValueError, match="duplicate"):
        hash_lock(tmp_path, {"Example": "1"})


def test_ambiguous_metadata_rejected(tmp_path):
    path = tmp_path / "example.whl"
    wheel(path, extra=True)
    with pytest.raises(ValueError, match="metadata"):
        wheel_metadata(path)


def test_no_overwrite(tmp_path):
    with pytest.raises(FileExistsError):
        build(tmp_path / "absent", "0" * 64, tmp_path)


def test_wrong_baseline_pin(tmp_path):
    path = tmp_path / "baseline.json"
    path.write_text("{}")
    output = tmp_path / "output"
    with pytest.raises(ValueError, match="checksum"):
        build(path, "0" * 64, output)
    assert not output.exists()
