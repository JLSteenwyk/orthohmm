"""Unit tests for orthohmm/search/prefilter.py (k-mer index + candidate gen)."""
import numpy as np
import pytest

from orthohmm.search.prefilter import (
    REDUCED_ALPHA,
    REDUCED_ALPHA_SIZE,
    _remap_flat,
    build_kmer_index,
    compute_freq_threshold,
    is_c_prefilter_available,
    prefilter_candidates,
    prepare_kmer_index,
)
from orthohmm.search.sequences import SpeciesSequences


def _make_species(tmp_path, records):
    p = tmp_path / "x.fa"
    p.write_text("".join(f">{h}\n{s}\n" for h, s in records))
    return SpeciesSequences.from_fasta(str(p), "x.fa")


# ─── C prefilter availability ───────────────────────────────────────────


class TestCPrefilterAvailability:
    def test_returns_bool(self):
        # Just exercise the lazy-load path; result depends on build env
        result = is_c_prefilter_available()
        assert isinstance(result, bool)


# ─── Reduced alphabet ───────────────────────────────────────────────────


class TestReducedAlphabet:
    def test_size_consistent(self):
        # 9 distinct groups
        assert REDUCED_ALPHA_SIZE == 9
        assert set(REDUCED_ALPHA.tolist()).issubset(set(range(REDUCED_ALPHA_SIZE)))

    def test_covers_20_amino_acids(self):
        assert REDUCED_ALPHA.shape == (20,)


# ─── _remap_flat ────────────────────────────────────────────────────────


class TestRemapFlat:
    def test_each_residue_mapped_through_table(self):
        # Encode 'AACE' = [0, 0, 1, 3]
        src = np.array([0, 0, 1, 3], dtype=np.uint8)
        out = _remap_flat(src, REDUCED_ALPHA, len(src))
        # Apply the mapping table by hand
        expected = REDUCED_ALPHA[src]
        assert np.array_equal(out, expected)

    def test_returns_uint8_array(self):
        src = np.array([0, 1, 2], dtype=np.uint8)
        out = _remap_flat(src, REDUCED_ALPHA, len(src))
        assert out.dtype == np.uint8
        assert len(out) == len(src)


# ─── build_kmer_index ───────────────────────────────────────────────────


class TestBuildKmerIndex:
    def test_csr_offsets_monotonic(self, tmp_path):
        sp = _make_species(tmp_path, [
            ("g1", "ACDEFGHIKLMN"),
            ("g2", "ACDEFGHIKL"),
        ])
        kmer_offsets, kmer_entries, kmer_freqs = build_kmer_index(
            sp.flat_sequences, sp.offsets, sp.lengths,
            k=3, alpha_size=20,
        )
        # CSR offsets are monotonically non-decreasing
        assert np.all(np.diff(kmer_offsets) >= 0)
        # Index size = 20^3 + 1
        assert kmer_offsets.shape == (20**3 + 1,)

    def test_kmer_freqs_sums_match_total_kmer_count(self, tmp_path):
        sp = _make_species(tmp_path, [
            ("g1", "ACDEFGHIK"),     # 9 residues → 7 3-mers
            ("g2", "MNPQRSTVW"),     # 9 residues → 7 3-mers
        ])
        _, kmer_entries, kmer_freqs = build_kmer_index(
            sp.flat_sequences, sp.offsets, sp.lengths,
            k=3, alpha_size=20,
        )
        # 9 + 9 residues → (9-3+1) + (9-3+1) = 7 + 7 = 14 total k-mers
        assert int(kmer_freqs.sum()) == 14
        # Posting list total length matches
        assert len(kmer_entries) == 14

    def test_freqs_count_unique_kmers_correctly(self, tmp_path):
        # All "AAA" k-mers from a single gene of length 5 ("AAAAA") = 3 instances
        sp = _make_species(tmp_path, [("g1", "AAAAA")])
        _, _, kmer_freqs = build_kmer_index(
            sp.flat_sequences, sp.offsets, sp.lengths,
            k=3, alpha_size=20,
        )
        # k-mer index for "AAA" = 0 * 20^2 + 0 * 20 + 0 = 0
        assert int(kmer_freqs[0]) == 3
        # All other 20^3 - 1 k-mers absent
        assert int(kmer_freqs.sum()) == 3


# ─── compute_freq_threshold ─────────────────────────────────────────────


class TestComputeFreqThreshold:
    def test_empty_kmer_freqs_returns_db_size(self):
        empty = np.zeros(100, dtype=np.int32)
        thresh = compute_freq_threshold(empty, num_sequences=500)
        assert int(thresh) == 500

    def test_threshold_floored_at_100_or_n_over_200(self):
        # All k-mers occur exactly once; 99.5th percentile = 1, but floor kicks in
        freqs = np.ones(1000, dtype=np.int32)
        thresh = compute_freq_threshold(freqs, num_sequences=10_000)
        # max(percentile=1, max(100, 10000/200=50)) = max(1, 100) = 100
        assert int(thresh) == 100

    def test_high_frequency_kmers_set_threshold(self):
        # A few k-mers are very common; 99.5th percentile picks up high values
        freqs = np.zeros(1000, dtype=np.int32)
        freqs[:990] = 1
        freqs[990:1000] = 500  # top 1% are very frequent
        thresh = compute_freq_threshold(freqs, num_sequences=200)
        # max(100, 200/200=1) = floor of 100; percentile catches the 500s
        # so threshold should be at least 100
        assert int(thresh) >= 100


class TestPreparedKmerIndex:
    def test_reuse_preserves_candidates(self, tmp_path):
        species = _make_species(tmp_path, [
            ("g1", "ACDEFGHIKLMNPQRSTVWY"),
            ("g2", "ACDEFGHIKLMNPQAAAAAA"),
        ])
        kwargs = dict(
            k=3,
            use_reduced_alphabet=False,
            min_total_hits=1,
            min_diag_hits=1,
            max_candidates_per_query=10,
        )
        direct = prefilter_candidates(species, species, **kwargs)
        prepared = prepare_kmer_index(
            species, k=3, use_reduced_alphabet=False
        )
        reused = prefilter_candidates(
            species, species, prepared_index=prepared, **kwargs
        )
        assert np.array_equal(reused[0], direct[0])
        assert np.array_equal(reused[1], direct[1])

    def test_rejects_mismatched_settings(self, tmp_path):
        species = _make_species(tmp_path, [("g1", "ACDEFGHIK")])
        prepared = prepare_kmer_index(
            species, k=3, use_reduced_alphabet=False
        )
        with pytest.raises(ValueError, match="does not match"):
            prefilter_candidates(
                species,
                species,
                k=4,
                use_reduced_alphabet=False,
                prepared_index=prepared,
            )
