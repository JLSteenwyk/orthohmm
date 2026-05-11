"""Unit tests for orthohmm/search/evalue.py (Karlin–Altschul math)."""
import numpy as np
import pytest

from orthohmm.search.evalue import batch_estimate_evalues, estimate_evalue
from orthohmm.search.matrices import get_ka_params


class TestEstimateEvalue:
    def test_known_blosum62_value(self):
        # E = K * m * n * exp(-lambda * S)
        lam, K = get_ka_params("BLOSUM62")
        m, n, S = 100, 10_000, 50
        e = estimate_evalue(S, m, n, lam, K)
        expected = K * m * n * np.exp(-lam * S)
        assert e == pytest.approx(expected, rel=1e-6)

    def test_zero_score_returns_large_sentinel(self):
        # Score <= 0 is the engine's "no hit" path; we return 1e10 to push
        # the E-value past any sane threshold.
        lam, K = get_ka_params("BLOSUM62")
        assert estimate_evalue(0, 100, 10_000, lam, K) == 1e10
        assert estimate_evalue(-5, 100, 10_000, lam, K) == 1e10

    def test_evalue_decreases_with_score(self):
        # Higher Viterbi score → smaller (better) E-value
        lam, K = get_ka_params("BLOSUM62")
        e_low  = estimate_evalue(10, 100, 10_000, lam, K)
        e_high = estimate_evalue(60, 100, 10_000, lam, K)
        assert e_high < e_low

    def test_evalue_increases_with_db_size(self):
        # Larger search space → more significant-by-chance hits → larger E
        lam, K = get_ka_params("BLOSUM62")
        e_small = estimate_evalue(40, 100,   1_000, lam, K)
        e_large = estimate_evalue(40, 100, 100_000, lam, K)
        assert e_large > e_small


class TestBatchEstimateEvalues:
    def test_matches_scalar_estimate_per_pair(self):
        lam, K = get_ka_params("BLOSUM62")
        scores       = np.array([10, 40,  80,  -1], dtype=np.int32)
        query_lens   = np.array([50, 100, 200, 100], dtype=np.int32)
        db_total     = 10_000

        batch = batch_estimate_evalues(scores, query_lens, db_total, lam, K)
        # Should match element-wise (within float32→64 round-trip noise)
        for i, s in enumerate(scores):
            expected = estimate_evalue(int(s), int(query_lens[i]), db_total, lam, K)
            assert batch[i] == pytest.approx(expected, rel=1e-6)

    def test_returns_correct_shape_and_dtype(self):
        lam, K = get_ka_params("BLOSUM62")
        scores = np.array([10, 40, 80], dtype=np.int32)
        ql     = np.array([50, 50, 50], dtype=np.int32)
        out = batch_estimate_evalues(scores, ql, 1000, lam, K)
        assert out.shape == (3,)
        assert out.dtype == np.float64

    def test_sentinel_for_non_positive_scores(self):
        lam, K = get_ka_params("BLOSUM62")
        scores = np.array([0, -1, 5], dtype=np.int32)
        ql     = np.array([100, 100, 100], dtype=np.int32)
        out = batch_estimate_evalues(scores, ql, 1000, lam, K)
        assert out[0] == 1e10
        assert out[1] == 1e10
        # Positive score gets a real value
        assert out[2] < 1e10
