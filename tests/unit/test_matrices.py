"""Unit tests for orthohmm/search/matrices.py (substitution matrix registry)."""
import numpy as np
import pytest

from orthohmm.search.matrices import (
    ALPHABET,
    ALPHABET_SIZE,
    KA_PARAMS,
    MATRICES,
    get_background_freqs,
    get_ka_params,
    get_matrix,
)


SUPPORTED = [
    "BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM90",
    "PAM30", "PAM70", "PAM120", "PAM240",
    "WAG", "LG",
]


class TestGetMatrix:
    @pytest.mark.parametrize("name", SUPPORTED)
    def test_is_20x20(self, name):
        m = get_matrix(name)
        assert m.shape == (ALPHABET_SIZE, ALPHABET_SIZE) == (20, 20)

    @pytest.mark.parametrize("name", SUPPORTED)
    def test_is_symmetric(self, name):
        m = get_matrix(name)
        assert np.array_equal(m, m.T), f"{name} matrix is not symmetric"

    def test_blosum62_diagonal_is_positive(self):
        # Each residue should align positively with itself.
        m = get_matrix("BLOSUM62")
        diag = np.diag(m)
        assert (diag > 0).all()

    def test_blosum62_known_entries(self):
        m = get_matrix("BLOSUM62")
        # BLOSUM62 cysteine-cysteine = 9 in the published matrix
        c = ALPHABET.index("C")
        assert m[c, c] == 9
        # Tryptophan-tryptophan = 11
        w = ALPHABET.index("W")
        assert m[w, w] == 11

    def test_unknown_matrix_raises(self):
        with pytest.raises((KeyError, ValueError)):
            get_matrix("NOT_A_REAL_MATRIX")

    def test_all_supported_names_in_registry(self):
        # Sanity: SUPPORTED list matches the MATRICES dict keys (case-insensitive
        # is OK; the registry keys are upper-case).
        registered = {k.upper() for k in MATRICES}
        for name in SUPPORTED:
            assert name in registered


class TestGetBackgroundFreqs:
    @pytest.mark.parametrize("name", SUPPORTED)
    def test_freqs_sum_to_one(self, name):
        f = get_background_freqs(name)
        assert f.shape == (ALPHABET_SIZE,)
        assert f.sum() == pytest.approx(1.0, abs=1e-3)

    @pytest.mark.parametrize("name", SUPPORTED)
    def test_freqs_all_positive(self, name):
        f = get_background_freqs(name)
        assert (f > 0).all(), f"{name} has non-positive background freqs"


class TestGetKaParams:
    @pytest.mark.parametrize("name", SUPPORTED)
    def test_returns_two_floats(self, name):
        lam, K = get_ka_params(name)
        assert isinstance(lam, float) and lam > 0
        assert isinstance(K, float) and K > 0

    def test_known_blosum62_values(self):
        # Karlin–Altschul ungapped BLOSUM62 published parameters
        lam, K = get_ka_params("BLOSUM62")
        # lambda ≈ 0.32, K ≈ 0.14 (BLAST defaults). Tolerance ±20% so any
        # standard literature value passes.
        assert 0.25 < lam < 0.40
        assert 0.05 < K < 0.30

    def test_registry_has_entry_for_every_supported_matrix(self):
        for name in SUPPORTED:
            assert name in KA_PARAMS, f"KA params missing for {name}"
