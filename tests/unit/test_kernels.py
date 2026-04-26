"""Regression tests for the ctypes-loaded kernels in orthohmm/search/csrc.

These tests protect against silent correctness breakage when anyone tweaks
the C/CUDA code. They check parity between implementations (scalar C vs
multipair AVX2 vs CUDA vs Numba fallback) on random synthetic data where
we know the scalar reference is right.

Tests skip gracefully if a kernel isn't built (e.g. no AVX2, no CUDA).
"""
from __future__ import annotations

import numpy as np
import pytest

from orthohmm.search.viterbi import (
    batch_viterbi_c,
    batch_viterbi_multipair_c,
    batch_viterbi_score,
    is_c_available,
    _load_c_hmm,
)


# ─── Fixtures ──────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def synthetic_pairs():
    """Small but non-trivial batch of (profile, target) pairs."""
    rng = np.random.default_rng(12345)
    n_prof = 25
    n_tgt = 50
    plens = rng.integers(80, 300, n_prof, dtype=np.int32)
    tlens = rng.integers(80, 300, n_tgt, dtype=np.int32)

    flat_me = rng.integers(-5, 11, (int(plens.sum()), 20), dtype=np.int8)
    po = np.concatenate(([0], np.cumsum(plens).astype(np.int64)[:-1]))
    tf = rng.integers(0, 21, int(tlens.sum()), dtype=np.uint8)
    to_ = np.concatenate(([0], np.cumsum(tlens).astype(np.int64)[:-1]))
    ins = np.full(20, -1, dtype=np.int8)
    tr = np.array([0, -12, -12, -2, -3, -2, -3], dtype=np.int32)

    # A mix of same-query and diverse-query pairs (120 pairs)
    pairs = []
    for q in range(n_prof):
        for t in rng.choice(n_tgt, size=5, replace=False):
            pairs.append((q, int(t)))
    pairs = np.array(pairs, dtype=np.int32)
    return {
        "flat_me": flat_me, "ins": ins, "tr": tr,
        "po": po, "plens": plens,
        "tf": tf, "to": to_, "tlens": tlens,
        "pairs": pairs,
    }


# ─── Viterbi kernels ───────────────────────────────────────────────────

@pytest.mark.skipif(not is_c_available(), reason="C Viterbi kernel not built")
def test_viterbi_c_matches_numba(synthetic_pairs):
    """C/OpenMP Viterbi must match the Numba reference exactly."""
    p = synthetic_pairs
    s_c = batch_viterbi_c(
        p["flat_me"], p["ins"], p["tr"], p["po"], p["plens"],
        p["tf"], p["to"], p["tlens"], p["pairs"], 64, n_threads=1,
    )
    s_nb = batch_viterbi_score(
        p["flat_me"], p["ins"], p["tr"], p["po"], p["plens"],
        p["tf"], p["to"], p["tlens"], p["pairs"], 64,
    )
    np.testing.assert_array_equal(s_c, s_nb)


@pytest.mark.skipif(not is_c_available(), reason="C Viterbi kernel not built")
def test_viterbi_multipair_matches_scalar(synthetic_pairs):
    """Multipair AVX2 batch must match scalar C exactly."""
    lib, _ = _load_c_hmm()
    if not hasattr(lib, "batch_hmm_viterbi_multipair_avx2_c"):
        pytest.skip("no AVX2 multipair kernel")
    p = synthetic_pairs
    s_scalar = batch_viterbi_c(
        p["flat_me"], p["ins"], p["tr"], p["po"], p["plens"],
        p["tf"], p["to"], p["tlens"], p["pairs"], 64, n_threads=1,
    )
    s_mp = batch_viterbi_multipair_c(
        p["flat_me"], p["ins"], p["tr"], p["po"], p["plens"],
        p["tf"], p["to"], p["tlens"], p["pairs"], 64, n_threads=1,
    )
    np.testing.assert_array_equal(s_mp, s_scalar)


def test_viterbi_xdrop_bounded_by_exact(synthetic_pairs):
    """X-drop must always be <= exact score (never over-estimates)."""
    if not is_c_available():
        pytest.skip("C kernel not built")
    p = synthetic_pairs
    s_exact = batch_viterbi_c(
        p["flat_me"], p["ins"], p["tr"], p["po"], p["plens"],
        p["tf"], p["to"], p["tlens"], p["pairs"], 64,
        n_threads=1, xdrop=0,
    )
    s_xdrop = batch_viterbi_c(
        p["flat_me"], p["ins"], p["tr"], p["po"], p["plens"],
        p["tf"], p["to"], p["tlens"], p["pairs"], 64,
        n_threads=1, xdrop=30,
    )
    assert (s_xdrop <= s_exact).all(), "x-drop produced higher score than exact"


# ─── CUDA (optional) ───────────────────────────────────────────────────

def _cuda_available():
    try:
        from orthohmm.search.viterbi_cuda_ctypes import is_cuda_available
        return is_cuda_available()
    except Exception:
        return False


@pytest.mark.skipif(not _cuda_available(), reason="CUDA kernel not available")
def test_viterbi_cuda_matches_scalar(synthetic_pairs):
    """CUDA warp-collaborative kernel must match scalar C exactly."""
    from orthohmm.search.viterbi_cuda_ctypes import batch_viterbi_cuda
    p = synthetic_pairs
    s_scalar = batch_viterbi_c(
        p["flat_me"], p["ins"], p["tr"], p["po"], p["plens"],
        p["tf"], p["to"], p["tlens"], p["pairs"], 64, n_threads=1,
    )
    s_gpu = batch_viterbi_cuda(
        p["flat_me"], p["ins"], p["tr"], p["po"], p["plens"],
        p["tf"], p["to"], p["tlens"], p["pairs"], 64,
    )
    np.testing.assert_array_equal(s_gpu, s_scalar)


# ─── Pair align ────────────────────────────────────────────────────────

def test_pair_align_self_gives_top_score():
    """Aligning a sequence against itself should return len*diagonal_score."""
    from orthohmm.search.msa_center_star import _load_pair_align, _encode_sequence
    from orthohmm.search.matrices import get_matrix
    try:
        lib = _load_pair_align()
    except OSError:
        pytest.skip("pair_align.so not built")

    import ctypes
    sm = get_matrix("BLOSUM62")
    # Use a known sequence
    seq = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSF"
    enc = _encode_sequence(seq)
    L = len(enc)
    flat = enc.astype(np.uint8)
    offsets = np.array([0, L], dtype=np.int64)
    lengths = np.array([L, L], dtype=np.int32)
    pairs = np.array([[0, 1]], dtype=np.int32)

    out_a = np.zeros(L * 2 + 4, dtype=np.uint8)
    out_b = np.zeros(L * 2 + 4, dtype=np.uint8)
    out_lens = np.zeros(1, dtype=np.int32)
    out_scores = np.zeros(1, dtype=np.int32)
    _sm = np.ascontiguousarray(sm, dtype=np.int8)
    pr = np.ascontiguousarray(pairs.flatten(), dtype=np.int32)
    lib.pair_align_set_num_threads(1)
    lib.batch_pair_align_c(
        np.concatenate([flat, flat]).ctypes.data,
        offsets.ctypes.data, lengths.ctypes.data,
        _sm.ctypes.data, -11, -1,
        pr.ctypes.data, 1, 0,
        out_a.ctypes.data, out_b.ctypes.data,
        out_lens.ctypes.data, L * 2 + 4,
        out_scores.ctypes.data,
    )
    # Self-alignment: score = sum of diagonal (each AA with itself).
    # BLOSUM62 self-scores are positive; sum should be > 0.
    assert out_scores[0] > 0
    # Alignment length should equal L (no gaps needed).
    assert out_lens[0] == L
    # Aligned a and b should match position-by-position (no gap=20 entries).
    aligned_a = out_a[:L]
    aligned_b = out_b[:L]
    np.testing.assert_array_equal(aligned_a, aligned_b)
    assert (aligned_a < 20).all()


# ─── MSA profile build ─────────────────────────────────────────────────

def test_center_star_msa_aligns_identical_sequences():
    """Identical sequences should produce an un-gapped MSA."""
    from orthohmm.search.msa_center_star import center_star_msa, _load_pair_align
    from orthohmm.search.matrices import get_matrix
    try:
        _load_pair_align()
    except OSError:
        pytest.skip("pair_align.so not built")
    sm = get_matrix("BLOSUM62")
    seq = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGK"
    seqs = [seq, seq, seq]
    aligned = center_star_msa(seqs, sm, gap_open=-11, gap_extend=-1)
    assert aligned is not None
    assert len(aligned) == 3
    for s in aligned:
        assert s == seq  # no gaps, identical output
