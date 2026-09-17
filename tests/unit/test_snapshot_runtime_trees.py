import pytest

from benchmark_tools.snapshot_runtime_trees import inventory, verify


def test_roundtrip_and_changed_file(tmp_path):
    source = tmp_path / "source.py"
    source.write_text("value = 1\n")
    frozen = inventory([tmp_path])
    assert verify(frozen)["status"] == "runtime_tree_identity_matches"
    source.write_text("value = 2\n")
    with pytest.raises(ValueError, match="changed"):
        verify(frozen)


def test_addition_deletion_and_mode(tmp_path):
    source = tmp_path / "source"
    source.write_text("data")
    frozen = inventory([tmp_path])
    extra = tmp_path / "extra"
    extra.write_text("new")
    with pytest.raises(ValueError, match="changed"):
        verify(frozen)
    extra.unlink()
    source.chmod(0o700)
    with pytest.raises(ValueError, match="changed"):
        verify(frozen)
    source.unlink()
    with pytest.raises(ValueError, match="changed"):
        verify(frozen)


def test_cache_excluded_but_source_included(tmp_path):
    source = tmp_path / "module.py"
    source.write_text("pass")
    frozen = inventory([tmp_path])
    cache = tmp_path / "__pycache__"
    cache.mkdir()
    (cache / "module.pyc").write_bytes(b"cache")
    assert verify(frozen)["records"] == 2


def test_external_directory_symlink_disclosed(tmp_path):
    root = tmp_path / "runtime"
    root.mkdir()
    external = tmp_path / "external"
    external.mkdir()
    (external / "file").write_text("outside")
    link = root / "link"
    link.symlink_to(external, target_is_directory=True)
    frozen = inventory([root])
    assert frozen["external_symlinks"] == [str(link)]
    assert len(frozen["records"]) == 2


def test_file_symlink_target_changes(tmp_path):
    root = tmp_path / "runtime"
    root.mkdir()
    target = tmp_path / "file"
    target.write_text("first")
    (root / "link").symlink_to(target)
    frozen = inventory([root])
    target.write_text("other")
    with pytest.raises(ValueError, match="changed"):
        verify(frozen)


def test_reject_empty_duplicate_overlap_and_missing(tmp_path):
    with pytest.raises(ValueError):
        inventory([])
    with pytest.raises(ValueError):
        inventory([tmp_path, tmp_path])
    child = tmp_path / "child"
    child.mkdir()
    with pytest.raises(ValueError, match="Overlapping"):
        inventory([tmp_path, child])
    with pytest.raises(FileNotFoundError):
        inventory([tmp_path / "missing"])
