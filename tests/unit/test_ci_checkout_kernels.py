from pathlib import Path
import ctypes
import os
import shutil
import subprocess
import sys
from types import SimpleNamespace

import pytest

from benchmark_tools import build_ci_checkout_kernels as module
from benchmark_tools.run_ci_cpu_wheel import stage


@pytest.fixture
def build_fixture(tmp_path, monkeypatch):
    root = tmp_path / "source"
    source = root / "orthohmm/search/csrc"
    source.mkdir(parents=True)
    for name in module.KERNELS:
        (source / name.replace(".so", ".c")).write_bytes(b"source fixture")
    output = tmp_path / "compiled"
    calls = []

    def run(argv, **kwargs):
        calls.append((argv, kwargs))
        assert kwargs["env"]["ORTHOHMM_CPU_TARGET"] == "baseline"
        assert kwargs["check"] is True and kwargs["timeout"] == 300
        compiled = output / "orthohmm/search/csrc"
        compiled.mkdir(parents=True)
        for name in module.KERNELS:
            (compiled / name).write_bytes(name.encode())

    class Scalar:
        restype = None
        def __call__(self):
            return 0

    monkeypatch.setattr(module.subprocess, "run", run)
    monkeypatch.setattr(module.ctypes, "CDLL", lambda _: SimpleNamespace(hmm_have_avx2=Scalar()))
    return root, source, output, calls


def test_verified_build_uses_fresh_baseline_and_exclusive_checkout_files(build_fixture):
    root, source, output, calls = build_fixture
    result = module.build(root, output)
    assert result["cpu_target"] == "baseline" and result["controlled_timing"] is False
    assert result["frozen_benchmark_runtime"] is False
    assert len(calls) == 1
    for name in module.KERNELS:
        assert (source / name).read_bytes() == name.encode()
        assert (source / name.replace(".so", ".c")).read_bytes() == b"source fixture"


@pytest.mark.parametrize("case", ["existing_output", "output_symlink", "inherited_library"])
def test_preflight_failure_never_builds(build_fixture, case):
    root, source, output, calls = build_fixture
    if case == "existing_output":
        output.mkdir()
    elif case == "output_symlink":
        output.symlink_to(root / "missing")
    else:
        (source / "old.so").write_bytes(b"inherited")
    with pytest.raises((FileExistsError, ValueError)):
        module.build(root, output)
    assert calls == []


@pytest.mark.parametrize("case", ["missing_library", "extra_library", "load_failure", "nonbaseline", "source_drift"])
def test_invalid_build_never_exposes_checkout_libraries(build_fixture, monkeypatch, case):
    root, source, output, _ = build_fixture
    original = module.subprocess.run
    def run(*args, **kwargs):
        original(*args, **kwargs)
        compiled = output / "orthohmm/search/csrc"
        if case == "missing_library":
            (compiled / module.KERNELS[0]).unlink()
        elif case == "extra_library":
            (compiled / "extra.so").write_bytes(b"extra")
        elif case == "source_drift":
            (source / "pair_align.c").write_bytes(b"changed")
    monkeypatch.setattr(module.subprocess, "run", run)
    if case == "load_failure":
        def failed_load(_):
            raise OSError("unloadable fixture")
        monkeypatch.setattr(module.ctypes, "CDLL", failed_load)
    elif case == "nonbaseline":
        class Nonbaseline:
            restype = None
            def __call__(self):
                return 1
        monkeypatch.setattr(module.ctypes, "CDLL", lambda _: SimpleNamespace(hmm_have_avx2=Nonbaseline()))
    with pytest.raises((ValueError, OSError)):
        module.build(root, output)
    assert not list(source.glob("*.so"))


def test_failed_compiler_propagates_without_checkout_libraries(build_fixture, monkeypatch):
    root, source, output, _ = build_fixture
    def fail(*args, **kwargs):
        raise subprocess.CalledProcessError(1, args[0])
    monkeypatch.setattr(module.subprocess, "run", fail)
    with pytest.raises(subprocess.CalledProcessError):
        module.build(root, output)
    assert not list(source.glob("*.so"))


def test_actual_build_from_staged_committed_sources(tmp_path, monkeypatch):
    if not shutil.which("gcc"):
        pytest.skip("Native GCC required")
    root = Path(__file__).resolve().parents[2]
    checkout = tmp_path / "source"
    original = stage(root, checkout)
    compiler_directory = str(Path(shutil.which("gcc")).parent)
    monkeypatch.setenv("PATH", compiler_directory + os.pathsep + os.defpath)
    result = module.build(checkout, tmp_path / "build")
    assert result["kernels"] == list(module.KERNELS)
    for name in module.KERNELS:
        ctypes.CDLL(str(checkout / "orthohmm/search/csrc" / name))
    env = {key: value for key, value in os.environ.items() if key not in {"PYTHONPATH", "PYTHONHOME"}}
    code = (
        "from pathlib import Path; import orthohmm; "
        "assert Path(orthohmm.__file__).resolve() == Path.cwd() / 'orthohmm/__init__.py'; "
        "from orthohmm.search.matrices import get_matrix; "
        "from orthohmm.search.msa_center_star import center_star_msa; "
        "sequence='MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSF'; "
        "assert center_star_msa([sequence]*3, get_matrix('BLOSUM62')) == [sequence]*3; "
        "print('fresh checkout native alignment verified')"
    )
    run = subprocess.run([sys.executable, "-c", code], cwd=checkout, env=env,
                         capture_output=True, text=True, check=True, timeout=30)
    assert run.stdout.strip() == "fresh checkout native alignment verified"
    for row in original["files"]:
        assert Path(row["original"]["path"]).read_bytes() == Path(row["staged"]["path"]).read_bytes()
