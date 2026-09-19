from pathlib import Path
import runpy
import subprocess

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


@pytest.mark.parametrize("backend", ["cpu", "cuda"])
@pytest.mark.parametrize("outcome", ["success", "partial_failure", "no_output_failure"])
def test_kernel_build_discards_failed_outputs(tmp_path, monkeypatch, backend, outcome):
    monkeypatch.setattr(setuptools, "setup", lambda **kwargs: None)
    namespace = runpy.run_path(str(Path(__file__).resolve().parents[2] / "setup.py"))
    builder = namespace[f"build_{backend}_kernels"]
    globals_ = builder.__globals__
    monkeypatch.setitem(globals_, "_have_gcc", lambda: True)
    monkeypatch.setitem(globals_, "_have_nvcc", lambda: True)
    monkeypatch.setitem(globals_, "_gcc_supports", lambda *args, **kwargs: True)
    kernels = namespace[f"{backend.upper()}_KERNELS"]
    expected = []
    for kernel in kernels:
        (tmp_path / kernel[0]).write_text("/* fixture */\n")
        expected.append(Path(kernel[0]).with_suffix(".so").name)

    def compile_command(argv, cwd=None):
        target = Path(argv[argv.index("-o") + 1])
        if outcome != "no_output_failure":
            target.write_bytes(b"complete" if outcome == "success" else b"partial")
        if outcome != "success":
            raise subprocess.CalledProcessError(1, argv)

    monkeypatch.setattr(subprocess, "check_call", compile_command)
    built = builder(tmp_path)
    if outcome == "success":
        assert built == expected
        assert all((tmp_path / name).read_bytes() == b"complete" for name in expected)
    else:
        assert built == []
        assert not list(tmp_path.glob("*.so"))
    assert all((tmp_path / kernel[0]).is_file() for kernel in kernels)
