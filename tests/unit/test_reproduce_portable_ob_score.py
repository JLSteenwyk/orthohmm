import json

import pytest

from benchmark_tools.reproduce_portable_ob_score import identity, validate_reference_ids, verify_files, worker


def staged(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    path = data / "input.fa"
    path.write_text(">gene\nACDE\n")
    return dict(files=[dict(relative_path="data/input.fa", role="fasta", **identity(path))])


def test_exact_staged_inventory(tmp_path):
    verify_files(tmp_path, staged(tmp_path))


def test_changed_input(tmp_path):
    manifest = staged(tmp_path)
    (tmp_path / "data/input.fa").write_text("changed")
    with pytest.raises(ValueError, match="Changed staged input"):
        verify_files(tmp_path, manifest)


def test_extra_input(tmp_path):
    manifest = staged(tmp_path)
    (tmp_path / "data/extra.fa").write_text("extra")
    with pytest.raises(ValueError, match="inventory differs"):
        verify_files(tmp_path, manifest)


@pytest.mark.parametrize("name", ["/etc/passwd", "data/../../escape"])
def test_unsafe_input(tmp_path, name):
    manifest = staged(tmp_path)
    manifest["files"][0]["relative_path"] = name
    with pytest.raises(ValueError, match="Unsafe or duplicate"):
        verify_files(tmp_path, manifest)


def test_duplicate_input(tmp_path):
    manifest = staged(tmp_path)
    manifest["files"] *= 2
    with pytest.raises(ValueError, match="Unsafe or duplicate"):
        verify_files(tmp_path, manifest)


def test_symlink_input(tmp_path):
    manifest = staged(tmp_path)
    path = tmp_path / "data/input.fa"
    target = tmp_path / "original.fa"
    path.rename(target)
    path.symlink_to(target)
    with pytest.raises(ValueError, match="Changed staged input"):
        verify_files(tmp_path, manifest)


def test_worker_preserves_existing_output(tmp_path):
    path = tmp_path / "score.json"
    path.write_text(json.dumps(dict(retained=True)))
    with pytest.raises(FileExistsError):
        worker(tmp_path, path)
    assert json.loads(path.read_text()) == dict(retained=True)


def test_blank_uncertainty_line_is_not_gene():
    uncertain = dict(group={"", "b"})
    validate_reference_ids(dict(group={"a", "b"}), uncertain, {"a", "b"})
    assert uncertain == dict(group={"", "b"})


def test_nonempty_unknown_reference_gene_rejected():
    with pytest.raises(ValueError, match="unknown nonempty gene"):
        validate_reference_ids(dict(group={"a"}), dict(group={"bad"}), {"a"})
