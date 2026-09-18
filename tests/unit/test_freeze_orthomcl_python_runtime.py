from pathlib import Path

import pytest

from benchmark_tools.freeze_orthomcl_python_runtime import minimal_roots, require_packages


@pytest.mark.parametrize("packages", [{}, {"biopython": "1.86"},
    {"biopython": "1.86", "numpy": "2.2.5"},
    {"biopython": "1.86", "numpy": "2.2.6", "extra": "1"}])
def test_wrong_packages_rejected(packages):
    with pytest.raises(ValueError):
        require_packages(packages)


def test_exact_packages():
    require_packages({"biopython": "1.86", "numpy": "2.2.6"})


def test_minimal_roots_preserve_uncovered_files(tmp_path):
    root = tmp_path / "env"
    root.mkdir()
    other = tmp_path / "library.so"
    assert minimal_roots([root / "bin", root, other, root]) == sorted([root, other], key=str)


def test_symlink_alias_does_not_duplicate_runtime_root(tmp_path):
    path = tmp_path / "real"
    path.mkdir()
    alias = tmp_path / "alias"
    alias.symlink_to(path, target_is_directory=True)
    assert minimal_roots([path, alias]) == [path]
