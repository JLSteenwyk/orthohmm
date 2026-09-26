import pytest

from benchmark_tools.probe_relocated_phylogeny_tools import copy_regular, inspect_trace


def test_copy_preserves_bytes_and_executable_mode(tmp_path):
    source = tmp_path / "original"
    source.write_bytes(b"fixture tool")
    source.chmod(0o755)
    target = tmp_path / "relocated/tool"
    result = copy_regular(source, target)
    assert target.read_bytes() == source.read_bytes()
    assert target.stat().st_mode & 0o777 == result["mode"] == 0o755
    assert result["original"]["sha256"] == result["relocated"]["sha256"]


def test_symlink_rejected(tmp_path):
    source = tmp_path / "original"
    source.write_text("fixture")
    link = tmp_path / "link"
    link.symlink_to(source)
    with pytest.raises(ValueError, match="nonsymlink"):
        copy_regular(link, tmp_path / "copy")


@pytest.mark.parametrize("trace,valid", [
    ('1 execve("/tmp/new/mafft", [], []) = 0\n', True),
    ('1 execve("/original/tools/mafft", [], []) = 0\n', False),
    ('1 execve("/tmp/new/mafft", [], []) = 0\n2 access("/original/tools/nope", F_OK) = -1 ENOENT\n', False),
    ('', False), ('1 openat(AT_FDCWD, "fixture", O_RDONLY) = 3\n', False),
])
def test_trace_rejects_original_access_even_if_failed(tmp_path, trace, valid):
    path = tmp_path / "trace"
    path.write_text(trace)
    if valid:
        assert inspect_trace(path, ["/original/tools"])["original_prefix_matches"] == 0
    else:
        with pytest.raises(ValueError):
            inspect_trace(path, ["/original/tools"])
