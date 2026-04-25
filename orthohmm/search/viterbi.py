"""3-state banded Viterbi scoring for profile HMM search.

Implements the Plan7 HMM Viterbi algorithm with Match, Insert, and
Delete states. Uses banding and row-toggling for efficiency.

Two implementations:
  1. Numba JIT (always available, fallback)
  2. C/OpenMP (fast path, lazy-loaded from csrc/hmm_viterbi.so)

The Viterbi recurrence for profile position i, target position j:

  V_M(i,j) = emit_M(i, target[j]) + max(
      V_M(i-1, j-1) + t_MM,
      V_I(i-1, j-1) + t_IM,
      V_D(i-1, j-1) + t_DM
  )
  V_I(i,j) = emit_I(target[j]) + max(
      V_M(i, j-1) + t_MI,
      V_I(i, j-1) + t_II
  )
  V_D(i,j) = max(
      V_M(i-1, j) + t_MD,
      V_D(i-1, j) + t_DD
  )

This is a local alignment variant: V_M can start from 0 at any position,
and the maximum V_M over all (i,j) is the Viterbi score.
"""

import ctypes
import logging
import os
import numpy as np
from numba import njit, prange, int32, int8, float32
from pathlib import Path as _Path

from .profile import T_MM, T_MI, T_MD, T_IM, T_II, T_DM, T_DD

_logger = logging.getLogger(__name__)

# ──────────────────────────────────────────────────────────────────────
# C/OpenMP Viterbi extension (lazy-loaded)
# ──────────────────────────────────────────────────────────────────────

_c_hmm_lib = None
_c_hmm_available = None


def _load_c_hmm():
    """Load the C/OpenMP HMM Viterbi shared library."""
    global _c_hmm_lib, _c_hmm_available
    if _c_hmm_available is not None:
        return _c_hmm_lib, _c_hmm_available
    try:
        so_path = _Path(__file__).resolve().parent / "csrc" / "hmm_viterbi.so"
        lib = ctypes.cdll.LoadLibrary(str(so_path))
        lib.batch_hmm_viterbi_c.restype = None
        lib.batch_hmm_viterbi_c.argtypes = [
            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
            ctypes.c_void_p, ctypes.c_void_p,
            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
            ctypes.c_void_p, ctypes.c_int32, ctypes.c_int32,
            ctypes.c_void_p,
        ]
        if hasattr(lib, "batch_hmm_viterbi_xdrop_c"):
            lib.batch_hmm_viterbi_xdrop_c.restype = None
            lib.batch_hmm_viterbi_xdrop_c.argtypes = [
                ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
                ctypes.c_void_p, ctypes.c_void_p,
                ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
                ctypes.c_void_p, ctypes.c_int32, ctypes.c_int32,
                ctypes.c_int32, ctypes.c_void_p,
            ]
        if hasattr(lib, "batch_hmm_viterbi_avx2_c"):
            lib.batch_hmm_viterbi_avx2_c.restype = None
            lib.batch_hmm_viterbi_avx2_c.argtypes = [
                ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
                ctypes.c_void_p, ctypes.c_void_p,
                ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
                ctypes.c_void_p, ctypes.c_int32, ctypes.c_int32,
                ctypes.c_int32, ctypes.c_void_p,
            ]
        if hasattr(lib, "batch_hmm_viterbi_int16_c"):
            lib.batch_hmm_viterbi_int16_c.restype = None
            lib.batch_hmm_viterbi_int16_c.argtypes = [
                ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
                ctypes.c_void_p, ctypes.c_void_p,
                ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
                ctypes.c_void_p, ctypes.c_int32, ctypes.c_int32,
                ctypes.c_int32, ctypes.c_void_p, ctypes.c_void_p,
            ]
        if hasattr(lib, "batch_hmm_viterbi_multipair_avx2_c"):
            lib.batch_hmm_viterbi_multipair_avx2_c.restype = None
            lib.batch_hmm_viterbi_multipair_avx2_c.argtypes = [
                ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
                ctypes.c_void_p, ctypes.c_void_p,
                ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
                ctypes.c_void_p, ctypes.c_int32, ctypes.c_int32,
                ctypes.c_void_p,
            ]
        if hasattr(lib, "hmm_have_avx2"):
            lib.hmm_have_avx2.restype = ctypes.c_int32
            lib.hmm_have_avx2.argtypes = []
        lib.hmm_set_num_threads.restype = None
        lib.hmm_set_num_threads.argtypes = [ctypes.c_int32]
        _c_hmm_lib = lib
        _c_hmm_available = True
    except (OSError, AttributeError):
        _c_hmm_lib = None
        _c_hmm_available = False
    return _c_hmm_lib, _c_hmm_available


def batch_viterbi_multipair_c(
    flat_match_emit, insert_emit, transitions,
    profile_offsets, profile_lengths,
    target_flat, target_offsets, target_lengths,
    pairs, band_width, n_threads=4,
):
    """Inter-pair SIMD Viterbi (AVX2). Each AVX2 int32 lane runs an
    independent DP for one pair — 8 pairs are processed in lockstep.

    Pairs are auto-sorted by (query_idx, target_length) internally so
    8-lane batches share a query profile AND have similar target
    lengths, keeping the per-lane band union tight. Scores are returned
    in the CALLER'S original pair order.

    Expected speedup: ~2-3x over the scalar banded kernel on realistic
    data (many pairs per query — the typical prefilter output pattern).
    """
    lib, available = _load_c_hmm()
    if not available or not hasattr(lib, "batch_hmm_viterbi_multipair_avx2_c"):
        raise RuntimeError("Multipair AVX2 Viterbi not available")
    lib.hmm_set_num_threads(n_threads)

    M = len(pairs)
    pairs = np.ascontiguousarray(pairs, dtype=np.int32)

    # Sort by (query_idx, target_length) so 8-lane batches share a
    # profile and have similar target lengths.
    tlens_per_pair = np.ascontiguousarray(target_lengths, dtype=np.int32)[pairs[:, 1]]
    sort_key = pairs[:, 0].astype(np.int64) * (1 << 32) + tlens_per_pair.astype(np.int64)
    order = np.argsort(sort_key, kind="stable")
    pairs_sorted = pairs[order]

    _me = np.ascontiguousarray(flat_match_emit, dtype=np.int8)
    _ie = np.ascontiguousarray(insert_emit, dtype=np.int8)
    _tr = np.ascontiguousarray(transitions, dtype=np.int32)
    _po = np.ascontiguousarray(profile_offsets, dtype=np.int64)
    _pl = np.ascontiguousarray(profile_lengths, dtype=np.int32)
    _tf = np.ascontiguousarray(target_flat, dtype=np.uint8)
    _to = np.ascontiguousarray(target_offsets, dtype=np.int64)
    _tl = np.ascontiguousarray(target_lengths, dtype=np.int32)
    _pr = np.ascontiguousarray(pairs_sorted.flatten(), dtype=np.int32)
    out_scores_sorted = np.empty(M, dtype=np.int32)

    lib.batch_hmm_viterbi_multipair_avx2_c(
        _me.ctypes.data, _ie.ctypes.data, _tr.ctypes.data,
        _po.ctypes.data, _pl.ctypes.data,
        _tf.ctypes.data, _to.ctypes.data, _tl.ctypes.data,
        _pr.ctypes.data, M, band_width,
        out_scores_sorted.ctypes.data,
    )

    # Unsort back to caller's pair order
    out_scores = np.empty(M, dtype=np.int32)
    out_scores[order] = out_scores_sorted
    return out_scores


def batch_viterbi_int16_c(  # pragma: no cover — experimental reference
    flat_match_emit, insert_emit, transitions,
    profile_offsets, profile_lengths,
    target_flat, target_offsets, target_lengths,
    pairs, band_width, n_threads=4, xdrop=0,
):
    """AVX2 int16 saturating Viterbi with int32 fallback on overflow.

    Returns (scores, n_overflow). Pairs whose max score would saturate
    int16 are automatically rescored with the int32 kernel. In practice
    this is empirically SLOWER than the scalar int32 kernel on this CPU
    (saturating adds have a latency penalty and the scalar I-state pass
    dominates); the implementation is here for GPU reuse and for
    architectures where int16 helps more.
    """
    lib, available = _load_c_hmm()
    if not available or not hasattr(lib, "batch_hmm_viterbi_int16_c"):
        raise RuntimeError("Int16 Viterbi not available")
    lib.hmm_set_num_threads(n_threads)
    M = len(pairs)
    _me = np.ascontiguousarray(flat_match_emit, dtype=np.int8)
    _ie = np.ascontiguousarray(insert_emit, dtype=np.int8)
    _tr = np.ascontiguousarray(transitions, dtype=np.int32)
    _po = np.ascontiguousarray(profile_offsets, dtype=np.int64)
    _pl = np.ascontiguousarray(profile_lengths, dtype=np.int32)
    _tf = np.ascontiguousarray(target_flat, dtype=np.uint8)
    _to = np.ascontiguousarray(target_offsets, dtype=np.int64)
    _tl = np.ascontiguousarray(target_lengths, dtype=np.int32)
    _pr = np.ascontiguousarray(pairs.flatten(), dtype=np.int32)
    out_scores = np.empty(M, dtype=np.int32)
    n_overflow = ctypes.c_int32(0)
    lib.batch_hmm_viterbi_int16_c(
        _me.ctypes.data, _ie.ctypes.data, _tr.ctypes.data,
        _po.ctypes.data, _pl.ctypes.data,
        _tf.ctypes.data, _to.ctypes.data, _tl.ctypes.data,
        _pr.ctypes.data, M, band_width, int(xdrop),
        out_scores.ctypes.data, ctypes.byref(n_overflow),
    )
    return out_scores, n_overflow.value


def batch_viterbi_c(
    flat_match_emit, insert_emit, transitions,
    profile_offsets, profile_lengths,
    target_flat, target_offsets, target_lengths,
    pairs, band_width, n_threads=4, xdrop=0, use_avx2=False,
):
    """C/OpenMP batch Viterbi scoring.

    If xdrop > 0, uses the X-drop variant that terminates a pair early
    when the running max score hasn't improved for `xdrop` rows.
    If use_avx2=True, uses the AVX2 two-pass kernel (~10% faster than
    scalar; most of the scalar kernel is already auto-vectorized or
    ILP-bound — set by empirical microbenchmark).

    Returns scores array of shape (M,) int32.
    """
    lib, available = _load_c_hmm()
    if not available:
        raise RuntimeError("C HMM Viterbi extension not available")

    lib.hmm_set_num_threads(n_threads)

    M = len(pairs)
    _me = np.ascontiguousarray(flat_match_emit, dtype=np.int8)
    _ie = np.ascontiguousarray(insert_emit, dtype=np.int8)
    _tr = np.ascontiguousarray(transitions, dtype=np.int32)
    _po = np.ascontiguousarray(profile_offsets, dtype=np.int64)
    _pl = np.ascontiguousarray(profile_lengths, dtype=np.int32)
    _tf = np.ascontiguousarray(target_flat, dtype=np.uint8)
    _to = np.ascontiguousarray(target_offsets, dtype=np.int64)
    _tl = np.ascontiguousarray(target_lengths, dtype=np.int32)
    _pr = np.ascontiguousarray(pairs.flatten(), dtype=np.int32)
    out_scores = np.empty(M, dtype=np.int32)

    if use_avx2 and hasattr(lib, "batch_hmm_viterbi_avx2_c") and \
            hasattr(lib, "hmm_have_avx2") and lib.hmm_have_avx2():
        lib.batch_hmm_viterbi_avx2_c(
            _me.ctypes.data, _ie.ctypes.data, _tr.ctypes.data,
            _po.ctypes.data, _pl.ctypes.data,
            _tf.ctypes.data, _to.ctypes.data, _tl.ctypes.data,
            _pr.ctypes.data, M, band_width, int(xdrop),
            out_scores.ctypes.data,
        )
    elif xdrop > 0 and hasattr(lib, "batch_hmm_viterbi_xdrop_c"):
        lib.batch_hmm_viterbi_xdrop_c(
            _me.ctypes.data, _ie.ctypes.data, _tr.ctypes.data,
            _po.ctypes.data, _pl.ctypes.data,
            _tf.ctypes.data, _to.ctypes.data, _tl.ctypes.data,
            _pr.ctypes.data, M, band_width, int(xdrop),
            out_scores.ctypes.data,
        )
    else:
        lib.batch_hmm_viterbi_c(
            _me.ctypes.data, _ie.ctypes.data, _tr.ctypes.data,
            _po.ctypes.data, _pl.ctypes.data,
            _tf.ctypes.data, _to.ctypes.data, _tl.ctypes.data,
            _pr.ctypes.data, M, band_width,
            out_scores.ctypes.data,
        )
    return out_scores


def is_c_available():
    """Check if C/OpenMP Viterbi is available."""
    _, available = _load_c_hmm()
    return available


NEG_INF = int32(-1000000)


@njit(cache=True)
def _viterbi_score_one(
    match_emit,       # (L, 20) int8
    insert_emit,      # (20,) int8
    transitions,      # (7,) int32
    target_seq,       # uint8 array
    target_len,       # int
    profile_len,      # int
    band_width,       # int
):
    """Banded 3-state Viterbi for one (profile, target) pair.

    Returns the Viterbi score (int32).
    """
    L = profile_len
    T = target_len

    if L == 0 or T == 0:
        return int32(0)

    t_mm = transitions[0]
    t_mi = transitions[1]
    t_md = transitions[2]
    t_im = transitions[3]
    t_ii = transitions[4]
    t_dm = transitions[5]
    t_dd = transitions[6]

    bw = band_width
    max_dim = max(L, T)
    if bw <= 0 or max_dim <= 50:
        bw = max_dim

    cols = T + 1

    # Row-toggled DP: two rows per state (previous and current)
    # M = match, I = insert, D = delete
    M_prev = np.full(cols, NEG_INF, dtype=np.int32)
    M_curr = np.full(cols, NEG_INF, dtype=np.int32)
    I_prev = np.full(cols, NEG_INF, dtype=np.int32)
    I_curr = np.full(cols, NEG_INF, dtype=np.int32)
    D_prev = np.full(cols, NEG_INF, dtype=np.int32)
    D_curr = np.full(cols, NEG_INF, dtype=np.int32)

    # Local alignment: M can start from 0
    max_score = int32(0)

    for i in range(1, L + 1):
        # Band: target positions to consider
        j_center = int32((i * T) // L)  # proportional mapping
        j_lo = max(1, j_center - bw)
        j_hi = min(T, j_center + bw)

        # Reset current row borders
        for j in range(j_lo - 1, j_hi + 2):
            if 0 <= j < cols:
                M_curr[j] = NEG_INF
                I_curr[j] = NEG_INF
                D_curr[j] = NEG_INF

        emit_row = match_emit[i - 1]  # 0-indexed profile position

        for j in range(j_lo, j_hi + 1):
            t_aa = target_seq[j - 1]  # 0-indexed target position

            # Match emission score
            if t_aa < 20:
                m_emit = int32(emit_row[t_aa])
                i_emit = int32(insert_emit[t_aa])
            else:
                m_emit = int32(0)
                i_emit = int32(-1)

            # V_M(i,j) = emit + max(V_M(i-1,j-1)+t_mm, V_I(i-1,j-1)+t_im, V_D(i-1,j-1)+t_dm)
            best_m = int32(0)  # local: can start fresh
            val = M_prev[j - 1] + t_mm
            if val > best_m:
                best_m = val
            val = I_prev[j - 1] + t_im
            if val > best_m:
                best_m = val
            val = D_prev[j - 1] + t_dm
            if val > best_m:
                best_m = val
            M_curr[j] = best_m + m_emit

            # V_I(i,j) = emit + max(V_M(i,j-1)+t_mi, V_I(i,j-1)+t_ii)
            best_i = NEG_INF
            val = M_curr[j - 1] + t_mi
            if val > best_i:
                best_i = val
            val = I_curr[j - 1] + t_ii
            if val > best_i:
                best_i = val
            if best_i > NEG_INF // 2:
                I_curr[j] = best_i + i_emit
            else:
                I_curr[j] = NEG_INF

            # V_D(i,j) = max(V_M(i-1,j)+t_md, V_D(i-1,j)+t_dd)
            best_d = NEG_INF
            val = M_prev[j] + t_md
            if val > best_d:
                best_d = val
            val = D_prev[j] + t_dd
            if val > best_d:
                best_d = val
            D_curr[j] = best_d

            # Track max score (local alignment)
            if M_curr[j] > max_score:
                max_score = M_curr[j]

        # Swap rows
        M_prev, M_curr = M_curr, M_prev
        I_prev, I_curr = I_curr, I_prev
        D_prev, D_curr = D_curr, D_prev

    return max_score


@njit(parallel=True, cache=True)
def batch_viterbi_score(
    flat_match_emit,     # (total_profile_pos, 20) int8
    insert_emit,         # (20,) int8
    transitions,         # (7,) int32
    profile_offsets,     # (N_queries,) int64
    profile_lengths,     # (N_queries,) int32
    target_flat,         # uint8
    target_offsets,      # (N_targets,) int64
    target_lengths,      # (N_targets,) int32
    pairs,               # (M, 2) int32: (query_idx, target_idx)
    band_width,          # int
):
    """Parallel Viterbi scoring for M (query, target) pairs.

    Returns scores array of shape (M,) int32.
    """
    M = pairs.shape[0]
    scores = np.empty(M, dtype=np.int32)

    for idx in prange(M):
        qi = pairs[idx, 0]
        ti = pairs[idx, 1]

        # Extract profile for this query
        p_off = profile_offsets[qi]
        p_len = profile_lengths[qi]
        p_emit = flat_match_emit[p_off:p_off + p_len]

        # Extract target sequence
        t_off = target_offsets[ti]
        t_len = target_lengths[ti]
        t_seq = target_flat[t_off:t_off + t_len]

        scores[idx] = _viterbi_score_one(
            p_emit, insert_emit, transitions,
            t_seq, int(t_len), int(p_len),
            band_width,
        )

    return scores


@njit(cache=True)
def viterbi_all_targets(
    match_emit,          # (L, 20) int8
    insert_emit,         # (20,) int8
    transitions,         # (7,) int32
    target_flat,         # uint8
    target_offsets,      # (N_targets,) int64
    target_lengths,      # (N_targets,) int32
    target_indices,      # int32 array of target indices to score
    profile_len,         # int
    band_width,          # int
):
    """Score one profile against multiple targets (sequential).

    Returns scores array matching target_indices length.
    """
    n = len(target_indices)
    scores = np.empty(n, dtype=np.int32)

    for k in range(n):
        ti = target_indices[k]
        t_off = target_offsets[ti]
        t_len = target_lengths[ti]
        t_seq = target_flat[t_off:t_off + t_len]

        scores[k] = _viterbi_score_one(
            match_emit, insert_emit, transitions,
            t_seq, int(t_len), profile_len,
            band_width,
        )

    return scores
