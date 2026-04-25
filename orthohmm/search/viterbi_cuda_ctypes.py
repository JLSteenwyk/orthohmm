"""Python wrapper for the warp-collaborative CUDA Viterbi kernel.

This replaces the earlier naive numba.cuda kernel (viterbi_cuda.py) with
a proper warp-collaborative implementation in raw CUDA C. Design:

  - 1 pair per warp (32 threads per block)
  - Shared-memory workspace (6 row buffers per block)
  - Parallel M/D computation across threads
  - Parallel-scan I state (algebraic reformulation)

Empirical: 2.4-2.7x faster than the CPU multipair AVX2 kernel on
realistic batches. 100% parity with scalar CPU reference.

Build:
  nvcc -O3 -Xcompiler -fPIC -arch=sm_89 -shared \\
       -o csrc/hmm_viterbi_cuda.so csrc/hmm_viterbi_cuda.cu
"""

from __future__ import annotations

import ctypes
from pathlib import Path as _Path
from typing import Optional

import numpy as np


_lib = None
_available: Optional[bool] = None


def _load():
    global _lib, _available
    if _available is not None:
        return _lib, _available
    so_path = _Path(__file__).resolve().parent / "csrc" / "hmm_viterbi_cuda.so"
    if not so_path.exists():
        _available = False
        return None, False
    try:
        lib = ctypes.cdll.LoadLibrary(str(so_path))
        lib.hmm_viterbi_cuda_run.restype = ctypes.c_int
        lib.hmm_viterbi_cuda_run.argtypes = [
            ctypes.c_void_p, ctypes.c_int64,   # fme, fme_n
            ctypes.c_void_p,                    # ie
            ctypes.c_void_p,                    # tr
            ctypes.c_void_p, ctypes.c_int32,   # po, n_profs
            ctypes.c_void_p,                    # pl
            ctypes.c_void_p, ctypes.c_int64,   # tf, tf_n
            ctypes.c_void_p, ctypes.c_int32,   # to, n_tgts
            ctypes.c_void_p,                    # tl
            ctypes.c_void_p, ctypes.c_int32,   # pairs, num_pairs
            ctypes.c_int32, ctypes.c_int32,     # band_width, max_T
            ctypes.c_void_p,                    # out_scores
        ]
        lib.hmm_viterbi_cuda_device_count.restype = ctypes.c_int
        lib.hmm_viterbi_cuda_device_count.argtypes = []
        n_dev = lib.hmm_viterbi_cuda_device_count()
        if n_dev <= 0:
            _available = False
            return None, False
        _lib = lib
        _available = True
    except OSError:
        _lib = None
        _available = False
    return _lib, _available


def is_cuda_available() -> bool:
    _, ok = _load()
    return ok


def batch_viterbi_cuda(
    flat_match_emit, insert_emit, transitions,
    profile_offsets, profile_lengths,
    target_flat, target_offsets, target_lengths,
    pairs, band_width,
):
    """Warp-collaborative CUDA batch Viterbi.

    Returns (scores,) int32 array, 1 per pair, in the caller's order.
    Raises RuntimeError if CUDA isn't available or if target lengths
    exceed the shared-memory cap (~1998 residues).
    """
    lib, ok = _load()
    if not ok:
        raise RuntimeError("CUDA Viterbi not available")

    M = len(pairs)
    _me = np.ascontiguousarray(flat_match_emit, dtype=np.int8)
    _ie = np.ascontiguousarray(insert_emit, dtype=np.int8)
    _tr = np.ascontiguousarray(transitions, dtype=np.int32)
    _po = np.ascontiguousarray(profile_offsets, dtype=np.int64)
    _pl = np.ascontiguousarray(profile_lengths, dtype=np.int32)
    _tf = np.ascontiguousarray(target_flat, dtype=np.uint8)
    _to = np.ascontiguousarray(target_offsets, dtype=np.int64)
    _tl = np.ascontiguousarray(target_lengths, dtype=np.int32)
    _pr = np.ascontiguousarray(np.asarray(pairs, dtype=np.int32).flatten())
    out_scores = np.empty(M, dtype=np.int32)

    max_T = int(_tl.max()) if len(_tl) else 0
    rc = lib.hmm_viterbi_cuda_run(
        _me.ctypes.data, _me.size,
        _ie.ctypes.data, _tr.ctypes.data,
        _po.ctypes.data, len(_pl), _pl.ctypes.data,
        _tf.ctypes.data, _tf.size,
        _to.ctypes.data, len(_tl), _tl.ctypes.data,
        _pr.ctypes.data, M,
        int(band_width), max_T,
        out_scores.ctypes.data,
    )
    if rc == -1:
        raise RuntimeError(
            f"CUDA Viterbi: max target length {max_T} exceeds shared-memory cap"
        )
    if rc != 0:
        raise RuntimeError(f"CUDA Viterbi launch failed (rc={rc})")
    return out_scores
