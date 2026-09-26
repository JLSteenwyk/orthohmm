from pathlib import Path

import pytest

from benchmark_tools.build_publication_fasttree import compile_command, run, validate_help


def test_fixed_baseline_compile_flags():
    command = compile_command()
    assert command[0] == "/usr/bin/gcc"
    assert "-march=x86-64" in command and "-mtune=generic" in command
    assert not any(flag in command for flag in ("-march=native", "-mavx2", "-DUSE_SINGLE", "-DOPENMP"))
    assert command[-2:] == ["FastTree.c", "-lm"]


@pytest.mark.parametrize("line", ["FastTree 2.2.0 Double precision:", "FastTree Version 2.2.0 Double precision"])
def test_supported_help_banner(line):
    validate_help(0, "", line + "\nhelp text\n")


@pytest.mark.parametrize("code,text", [(1, "FastTree 2.2.0 Double precision:"),
    (0, "FastTree 2.1.11 Double precision:"), (0, "FastTree 2.2.0 SSE3:"), (0, "")])
def test_wrong_help_banner(code, text):
    with pytest.raises(ValueError):
        validate_help(code, "", text)


def test_refuse_existing_directory(tmp_path):
    with pytest.raises(FileExistsError):
        run(tmp_path, tmp_path / "absent", tmp_path / "mafft", tmp_path)


def test_refuse_symlink_output(tmp_path):
    output = tmp_path / "output"
    output.symlink_to(tmp_path / "missing")
    with pytest.raises(FileExistsError):
        run(tmp_path, tmp_path / "absent", tmp_path / "mafft", output)


def test_changed_source_rejected_before_output_creation(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    (source / "FastTree.c").write_text("not the pinned source")
    output = tmp_path / "output"
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        run(tmp_path, source, Path("/unused"), output)
    assert not output.exists()
