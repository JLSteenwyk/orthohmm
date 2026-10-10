import os
from pathlib import Path

import numpy as np
import pytest

from benchmark_tools.fragment_trace_artifact_access import ArtifactAccess, identity


def mapped(tmp_path, names):
    rows = []
    for i, (logical, data) in enumerate(names.items()):
        target = tmp_path / str(i)
        target.write_bytes(data)
        rows.append(dict(logical_path=logical, target=str(i), **identity(target)))
    return ArtifactAccess(tmp_path, rows)


def test_reads_use_copy_without_changing_historical_label(tmp_path):
    access = mapped(tmp_path, {"/old/data.txt": b"copied"})
    path = access.path("/old/data.txt")
    assert str(path.resolve()) == "/old/data.txt"
    assert path.read_text() == "copied"
    assert os.fspath(path) == str(tmp_path / "0")
    assert str(access.path(Path(os.fspath(path)))) == str(path)


def test_unindexed_real_original_file_is_not_a_fallback(tmp_path):
    original = tmp_path / "original"
    original.write_text("must not read")
    access = mapped(tmp_path, {"/old/data": b"valid"})
    with pytest.raises(ValueError, match="Unindexed"):
        access.path(original).read_text()


@pytest.mark.parametrize("mode", ["w", "a", "x", "r+", "wb"])
def test_writes_are_refused(tmp_path, mode):
    access = mapped(tmp_path, {"/old/data": b"valid"})
    with pytest.raises(ValueError, match="read-only"):
        access.path("/old/data").open(mode)


def test_mutated_copied_bytes_are_refused(tmp_path):
    access = mapped(tmp_path, {"/old/data": b"valid"})
    (tmp_path / "0").write_text("changed")
    with pytest.raises(ValueError, match="differs"):
        access.path("/old/data").read_bytes()


def test_two_logical_names_cannot_share_one_copied_target(tmp_path):
    with pytest.raises(ValueError, match="share one target"):
        ArtifactAccess(tmp_path, [dict(logical_path=name, target="one", bytes=0, sha256="0" * 64)
                                  for name in ("/old/a", "/old/b")])


def test_child_paths_names_and_shallow_vs_recursive_globs(tmp_path):
    access = mapped(tmp_path, {"/old/input/a.fasta": b"a", "/old/input/deep/b.fasta": b"b", "/old/input/deep/table.tsv": b"t"})
    parent = access.path("/old/input")
    assert [p.name for p in parent.glob("*.fasta")] == ["a.fasta"]
    assert [p.name for p in parent.glob("**/*.fasta")] == ["b.fasta"]
    assert (parent / "a.fasta").stem == "a"
    assert str((parent / "a.fasta").with_name("z.fasta")) == "/old/input/z.fasta"


def test_numpy_pathlike_reader_uses_copied_array(tmp_path):
    file = tmp_path / "array.npy"
    np.save(file, np.array([1, 2, 3]))
    access = mapped(tmp_path, {"/old/array.npy": file.read_bytes()})
    assert np.load(access.path("/old/array.npy"), allow_pickle=False).tolist() == [1, 2, 3]


def test_symlinked_copy_or_parent_escape_is_refused(tmp_path):
    target = tmp_path / "real"
    target.write_text("valid")
    alias = tmp_path / "alias"
    alias.symlink_to(target)
    access = ArtifactAccess(tmp_path, [dict(logical_path="/old/data", target="alias", **identity(target))])
    with pytest.raises(ValueError, match="escaped"):
        access.verify()


@pytest.mark.parametrize("logical,target", [("relative", "data"), ("/old/../data", "data"), ("/old/data", "../escape"), ("/old/data", "/absolute")])
def test_unsafe_mapping_is_refused(tmp_path, logical, target):
    with pytest.raises(ValueError):
        ArtifactAccess(tmp_path, [dict(logical_path=logical, target=target, bytes=0, sha256="0" * 64)])
