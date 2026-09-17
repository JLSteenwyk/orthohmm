from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import build_publication_arm_runtime as arm


@pytest.mark.parametrize("name", arm.KERNELS)
def test_arm_flags_and_exact_sources(name):
    command = arm.compile_command("/usr/bin/gcc", Path("/core/csrc"), name, "aarch64")
    assert command == ["/usr/bin/gcc", "-O3", "-fopenmp", "-shared", "-fPIC", "-march=armv8-a",
                       "-o", "/core/csrc/" + name + ".so", "/core/csrc/" + name + ".c"]


@pytest.mark.parametrize("machine,name", [("x86_64", "hmm_viterbi"), ("aarch64", "cuda")])
def test_wrong_architecture_or_kernel_rejected(machine, name):
    with pytest.raises(ValueError):
        arm.compile_command("gcc", Path("/core"), name, machine)


def test_no_overwrite(tmp_path):
    with pytest.raises(FileExistsError):
        arm.build(tmp_path, tmp_path)


def test_build_requires_arm(monkeypatch, tmp_path):
    monkeypatch.setattr(arm.platform, "machine", lambda: "x86_64")
    with pytest.raises(ValueError, match="AArch64"):
        arm.build(tmp_path, tmp_path / "report.json")
    assert not (tmp_path / "report.json").exists()


def test_missing_symbol_rejected(monkeypatch):
    monkeypatch.setattr(arm.ctypes, "CDLL", lambda path: SimpleNamespace())
    with pytest.raises(AttributeError):
        arm.inspect_library("/missing.so", "pair_align")


@pytest.mark.parametrize("avx2", [0, 1])
def test_hmm_backend_checked(monkeypatch, avx2):
    lib = SimpleNamespace(**{name: lambda: avx2 for name in arm.SYMBOLS["hmm_viterbi"]})
    monkeypatch.setattr(arm.ctypes, "CDLL", lambda path: lib)
    if avx2:
        with pytest.raises(ValueError, match="AVX2"):
            arm.inspect_library("/library.so", "hmm_viterbi")
    else:
        assert arm.inspect_library("/library.so", "hmm_viterbi")["hmm_have_avx2"] == 0
