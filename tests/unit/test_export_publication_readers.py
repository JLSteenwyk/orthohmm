import json

import pytest

from benchmark_tools.export_publication_readers import closure, dependencies, identity, verify


def test_dependency_syntax_and_function_local_imports():
    assert dependencies("from benchmark_tools import a as x\nimport os\ndef f():\n from benchmark_tools.b import y\n") == [
        "benchmark_tools.a", "benchmark_tools.b"]


def test_relative_import_rejected():
    with pytest.raises(ValueError):
        dependencies("from . import x")


def test_cycle_closure():
    sources = {"benchmark_tools/audit_publication_pipeline.py": b"from benchmark_tools import b",
               "benchmark_tools/b.py": b"from benchmark_tools import audit_publication_pipeline"}
    assert closure(sources.__getitem__) == sources


@pytest.mark.parametrize("problem", [None, "changed", "missing", "extra", "symlink"])
def test_export_integrity(tmp_path, problem):
    path = tmp_path / "a.py"
    path.write_bytes(b"x")
    (tmp_path / "manifest.json").write_text(json.dumps(dict(files={"a.py": identity(b"x")})))
    if problem == "changed":
        path.write_bytes(b"y")
    elif problem == "missing":
        path.unlink()
    elif problem == "extra":
        (tmp_path / "extra").write_text("x")
    elif problem == "symlink":
        path.unlink()
        path.symlink_to(tmp_path / "manifest.json")
    if problem is None:
        verify(tmp_path)
    else:
        with pytest.raises(ValueError):
            verify(tmp_path)
