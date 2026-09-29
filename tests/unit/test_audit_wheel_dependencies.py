import hashlib
import zipfile

import pytest
from packaging.markers import default_environment

from benchmark_tools.audit_wheel_dependencies import evaluate, metadata


def package(name="root", requires=None, python=None):
    return dict(name=name, version="1.0", requires_python=python, requires_dist=requires or [])


@pytest.mark.parametrize("requirement,status", [
    ("target>=1", "satisfied"), ("target>=2", "version_mismatch"),
    ("absent", "missing"), ("target[extra]", "unsupported_url_or_extras"),
    ("target @ https://example.org/a.whl", "unsupported_url_or_extras"),
    ('absent; extra == "test"', "inactive_marker")])
def test_declared_constraints(requirement, status):
    result = evaluate([package(requires=[requirement]), package("target")], default_environment())
    assert result["dependency_rows"][0]["status"] == status
    assert result["declared_runtime_dependencies_satisfied"] == (status in ("satisfied", "inactive_marker"))


def test_python_and_duplicate_rejection():
    result = evaluate([package(python=">=999")], default_environment())
    assert result["failures"][0]["reason"] == "requires_python"
    with pytest.raises(ValueError):
        evaluate([package(), package()], default_environment())


@pytest.mark.parametrize("case", ["valid", "hash", "name", "duplicate"])
def test_wheel_metadata_binding(tmp_path, case):
    path = tmp_path / "root-1.0-py3-none-any.whl"
    content = "Name: root\nVersion: 1.0\nRequires-Dist: target>=1\n\n"
    if case == "name":
        content = content.replace("Name: root", "Name: other")
    if case == "duplicate":
        content = "Name: root\n" + content
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("root-1.0.dist-info/METADATA", content)
        archive.writestr("root/_vendor/other-1.0.dist-info/METADATA", "Name: other\nVersion: 1.0\n")
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    pin = dict(name="root", version="1.0", hashes=["0" * 64 if case == "hash" else sha])
    if case == "valid":
        assert metadata(path, pin)["requires_dist"] == ["target>=1"]
    else:
        with pytest.raises(ValueError):
            metadata(path, pin)
