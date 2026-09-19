from pathlib import Path
import runpy

import setuptools
from setuptools.command.build_py import build_py
import pytest


def test_ctypes_distribution_is_not_platform_independent(monkeypatch):
    captured = {}
    monkeypatch.setattr(setuptools, "setup", lambda **kwargs: captured.update(kwargs))
    runpy.run_path(str(Path(__file__).resolve().parents[2] / "setup.py"))
    distribution = captured["distclass"]()
    assert distribution.has_ext_modules()
    assert not distribution.is_pure()
    command = distribution.get_command_obj("bdist_wheel")
    command.ensure_finalized()
    assert command.root_is_pure is False
    assert command.get_tag()[2] != "any"
    assert captured["entry_points"]["console_scripts"] == [
        "orthohmm = orthohmm.orthohmm:main"
    ]


@pytest.mark.parametrize("compiler_available", [False, True])
def test_wheel_build_excludes_stale_native_files(tmp_path, monkeypatch, compiler_available):
    captured = {}
    monkeypatch.setattr(setuptools, "setup", lambda **kwargs: captured.update(kwargs))
    runpy.run_path(str(Path(__file__).resolve().parents[2] / "setup.py"))
    command_class = captured["cmdclass"]["build_py"]
    command = command_class(captured["distclass"]())
    command.build_lib = str(tmp_path / "build")
    source = tmp_path / "source"
    source.mkdir()
    (source / "hmm_viterbi.so").write_bytes(b"stale source binary")
    target = Path(command.build_lib) / "orthohmm" / "search" / "csrc"
    target.mkdir(parents=True)
    (target / "old_kernel.so").write_bytes(b"stale build binary")

    def copy_package(_command):
        (target / "hmm_viterbi.c").write_text("/* source */\n")
        (target / "hmm_viterbi.so").write_bytes(b"stale copied binary")

    calls = []

    def compile_cpu(directory):
        assert directory == target
        assert not list(directory.glob("*.so"))
        assert (directory / "hmm_viterbi.c").is_file()
        calls.append("cpu")
        if compiler_available:
            (directory / "hmm_viterbi.so").write_bytes(b"fresh binary")
        return ["hmm_viterbi.so"] if compiler_available else []

    def compile_cuda(directory):
        assert directory == target
        calls.append("cuda")
        return []

    monkeypatch.setattr(build_py, "run", copy_package)
    monkeypatch.setitem(command_class.run.__globals__, "CSRC", source)
    monkeypatch.setitem(command_class.run.__globals__, "build_cpu_kernels", compile_cpu)
    monkeypatch.setitem(command_class.run.__globals__, "build_cuda_kernels", compile_cuda)
    command.run()
    assert calls == ["cpu", "cuda"]
    assert (source / "hmm_viterbi.so").read_bytes() == b"stale source binary"
    assert not (target / "old_kernel.so").exists()
    if compiler_available:
        assert (target / "hmm_viterbi.so").read_bytes() == b"fresh binary"
    else:
        assert not list(target.glob("*.so"))
