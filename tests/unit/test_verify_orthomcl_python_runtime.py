import copy

import pytest

from benchmark_tools.verify_orthomcl_python_runtime import compare, verify_runtime


def fixture():
    observed = {"python_version": "3.10.13", "prefix": "/env", "base_prefix": "/base",
                "executable": {"path": "/base/python", "bytes": 10, "sha256": "x"},
                "packages": {"biopython": "1.86", "numpy": "2.2.6"},
                "mapped_files": [{"path": "/lib/library.so", "bytes": 20, "sha256": "y"}]}
    expected = {"inspection": copy.deepcopy(observed),
                "runtime": {"records": [{"kind": "file", **observed["mapped_files"][0]}]}}
    return expected, observed


def test_exact_identity_and_mapped_subset():
    expected, observed = fixture()
    compare(expected, observed)
    observed["mapped_files"] = []
    compare(expected, observed)


@pytest.mark.parametrize("key", ["python_version", "prefix", "base_prefix", "executable", "packages"])
def test_changed_identity_rejected(key):
    expected, observed = fixture()
    observed[key] = "wrong"
    with pytest.raises(ValueError, match="identity"):
        compare(expected, observed)


@pytest.mark.parametrize("key,value", [("path", "/other/library.so"), ("bytes", 21), ("sha256", "z")])
def test_changed_or_unbound_mapped_library_rejected(key, value):
    expected, observed = fixture()
    observed["mapped_files"][0][key] = value
    with pytest.raises(ValueError, match="mapped"):
        compare(expected, observed)


def test_changed_manifest_rejected_before_runtime_inspection(tmp_path):
    directory = tmp_path / "benchmark_tools/results"
    directory.mkdir(parents=True)
    (directory / "orthomcl_python_runtime_20260918.json").write_text("{}")
    with pytest.raises(ValueError, match="manifest"):
        verify_runtime(tmp_path)
