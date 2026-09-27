import pytest

from benchmark_tools.probe_ob_leiden_overlay import allowed_distribution_path


@pytest.mark.parametrize("path", ["leidenalg/functions.py", "leidenalg/_c_leiden.abi3.so",
    "leidenalg.libs/libigraph.so", "leidenalg-0.11.0.dist-info/METADATA"])
def test_allowed_distribution_files(path):
    assert allowed_distribution_path(path)


@pytest.mark.parametrize("path", ["/leidenalg/a", "../leidenalg/a", "leidenalg/../../a",
    "numpy/__init__.py", "leidenalg-0.12.0.dist-info/METADATA", "."])
def test_reject_foreign_paths(path):
    with pytest.raises(ValueError):
        allowed_distribution_path(path)


def test_bytecode_not_copied():
    assert not allowed_distribution_path("leidenalg/__pycache__/functions.cpython-310.pyc")
