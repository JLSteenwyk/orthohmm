"""Compile current source so installed/stale binaries cannot mask regressions."""

import ctypes
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pytest

from orthohmm.search import viterbi


@pytest.fixture(scope="module")
def compiled_kernel(tmp_path_factory):
    compiler = shutil.which("gcc")
    if sys.platform != "linux" or compiler is None:
        pytest.skip("Native regression requires Linux GCC/OpenMP")
    source = Path(viterbi.__file__).parent / "csrc/hmm_viterbi.c"
    library = tmp_path_factory.mktemp("banding_kernel") / "hmm_viterbi.so"
    subprocess.run([compiler, "-O3", "-fopenmp", "-shared", "-fPIC", "-march=native",
                    str(source), "-o", str(library)], check=True, capture_output=True)
    lib = ctypes.CDLL(str(library))
    lib.hmm_have_avx2.argtypes = []
    lib.hmm_have_avx2.restype = ctypes.c_int32
    if not lib.hmm_have_avx2():
        pytest.skip("Native host has no compiled AVX2 backend")
    lib.hmm_set_num_threads.argtypes = [ctypes.c_int32]
    lib.hmm_set_num_threads.restype = None
    for name in ("batch_hmm_viterbi_c", "batch_hmm_viterbi_multipair_avx2_c"):
        function = getattr(lib, name)
        function.argtypes = [ctypes.c_void_p] * 9 + [ctypes.c_int32, ctypes.c_int32, ctypes.c_void_p]
        function.restype = None
    return lib


@pytest.fixture(scope="module")
def boundary_data():
    rng = np.random.Generator(np.random.PCG64(20260920))
    lengths = np.array([1, 31, 49, 50, 51, 64], dtype=np.int32)
    # Sorting still leaves a 51-residue target in the first eight-lane batch.
    target_lengths = np.array([0, 1, 8, 9, 31, 49, 50, 51, 52, 64, 129, 130, 200], dtype=np.int32)
    emissions = rng.integers(-5, 11, (int(lengths.sum()), 20), dtype=np.int8)
    insert = np.full(20, -1, dtype=np.int8)
    transitions = np.array([0, -12, -12, -2, -3, -2, -3], dtype=np.int32)
    offsets = np.concatenate(([0], np.cumsum(lengths, dtype=np.int64)[:-1]))
    target_offsets = np.concatenate(([0], np.cumsum(target_lengths, dtype=np.int64)[:-1]))
    targets = rng.integers(0, 21, int(target_lengths.sum()), dtype=np.uint8)
    pairs = np.array([(q, t) for q in range(len(lengths)) for t in range(len(target_lengths))], dtype=np.int32)
    return emissions, insert, transitions, offsets, lengths, targets, target_offsets, target_lengths, pairs


@pytest.mark.parametrize("band", [0, 1, 8, 64, 128])
@pytest.mark.parametrize("order", ["original", "reverse", "shuffle"])
@pytest.mark.parametrize("threads", [1, 4])
def test_per_pair_band_semantics(compiled_kernel, boundary_data, monkeypatch, band, order, threads):
    monkeypatch.setattr(viterbi, "_c_hmm_lib", compiled_kernel)
    monkeypatch.setattr(viterbi, "_c_hmm_available", True)
    pairs = boundary_data[-1]
    reference = viterbi.batch_viterbi_score(*boundary_data, band)
    scalar = viterbi.batch_viterbi_c(*boundary_data, band, n_threads=1)
    np.testing.assert_array_equal(scalar, reference)
    indices = np.arange(len(pairs))
    if order == "reverse":
        indices = indices[::-1].copy()
    elif order == "shuffle":
        np.random.Generator(np.random.PCG64(20260921)).shuffle(indices)
    actual = viterbi.batch_viterbi_multipair_c(*boundary_data[:-1], pairs[indices], band, n_threads=threads)
    np.testing.assert_array_equal(actual, reference[indices])
