"""Unit tests for orthohmm/search/engine.py (built-in search orchestrator)."""
import numpy as np
import pytest

from orthohmm.search.engine import (
    SpeciesPairResults,
    _filter_significant,
    execute_builtin_search,
    results_to_phmmer_format,
    search_species_pair,
)
from orthohmm.search.sequences import SpeciesSequences


# ─── SpeciesPairResults dataclass ───────────────────────────────────────


class TestSpeciesPairResults:
    def _make(self, n=3):
        return SpeciesPairResults(
            query_species="q.fa",
            target_species="t.fa",
            target_names=np.array([f"t{i}" for i in range(n)], dtype="U50"),
            query_names=np.array([f"q{i}" for i in range(n)], dtype="U50"),
            evalues=np.full(n, 1e-10, dtype=np.float64),
            scores=np.full(n, 100.0, dtype=np.float64),
        )

    def test_field_shapes_match(self):
        r = self._make(5)
        for arr in (r.target_names, r.query_names, r.evalues, r.scores):
            assert len(arr) == 5

    def test_holds_string_names(self):
        r = self._make(2)
        assert r.target_names[0] == "t0"
        assert r.query_names[1] == "q1"


# ─── _filter_significant ────────────────────────────────────────────────


class TestFilterSignificant:
    def test_returns_object_unchanged_when_all_below_threshold(self):
        r = SpeciesPairResults(
            query_species="q.fa", target_species="t.fa",
            target_names=np.array(["t0", "t1"], dtype="U50"),
            query_names=np.array(["q0", "q1"], dtype="U50"),
            evalues=np.array([1e-10, 1e-50]),
            scores=np.array([100.0, 200.0]),
        )
        out = _filter_significant(r, 1e-4)
        # When all pass, the function short-circuits and returns the input
        assert out is r

    def test_drops_hits_above_threshold(self):
        r = SpeciesPairResults(
            query_species="q.fa", target_species="t.fa",
            target_names=np.array(["t0", "t1", "t2"], dtype="U50"),
            query_names=np.array(["q0", "q1", "q2"], dtype="U50"),
            evalues=np.array([1e-10, 1e-3, 1e-50]),
            scores=np.array([100.0, 30.0, 200.0]),
        )
        out = _filter_significant(r, 1e-4)
        # 1e-3 fails (not < 1e-4); other two pass
        assert len(out.evalues) == 2
        assert list(out.target_names) == ["t0", "t2"]
        # New object, not the original
        assert out is not r

    def test_zero_hits_returns_empty_arrays(self):
        r = SpeciesPairResults(
            query_species="q.fa", target_species="t.fa",
            target_names=np.array(["t0"], dtype="U50"),
            query_names=np.array(["q0"], dtype="U50"),
            evalues=np.array([1.0]),       # well above threshold
            scores=np.array([5.0]),
        )
        out = _filter_significant(r, 1e-4)
        assert len(out.evalues) == 0
        assert len(out.target_names) == 0


# ─── results_to_phmmer_format ───────────────────────────────────────────


class TestResultsToPhmmerFormat:
    def test_empty_results_returns_empty_structured_array(self):
        r = SpeciesPairResults(
            query_species="q.fa", target_species="t.fa",
            target_names=np.empty(0, dtype="U50"),
            query_names=np.empty(0, dtype="U50"),
            evalues=np.empty(0),
            scores=np.empty(0),
        )
        out = results_to_phmmer_format(r, 1e-4)
        assert out.shape == (0,)
        assert out.dtype.names == ("target_name", "query_name", "evalue", "score")

    def test_filters_by_evalue_threshold(self):
        r = SpeciesPairResults(
            query_species="q.fa", target_species="t.fa",
            target_names=np.array(["t0", "t1", "t2"], dtype="U50"),
            query_names=np.array(["q0", "q1", "q2"], dtype="U50"),
            evalues=np.array([1e-10, 1e-3, 1e-50]),
            scores=np.array([100.0, 30.0, 200.0]),
        )
        out = results_to_phmmer_format(r, 1e-4)
        assert out.shape == (2,)
        names = sorted(out["target_name"].tolist())
        assert names == ["t0", "t2"]

    def test_field_dtypes(self):
        r = SpeciesPairResults(
            query_species="q.fa", target_species="t.fa",
            target_names=np.array(["t0"], dtype="U50"),
            query_names=np.array(["q0"], dtype="U50"),
            evalues=np.array([1e-10]),
            scores=np.array([100.0]),
        )
        out = results_to_phmmer_format(r, 1e-4)
        # Names are unicode strings up to 50 chars; numbers are float
        assert out.dtype["target_name"].kind == "U"
        assert out.dtype["query_name"].kind == "U"
        assert out.dtype["evalue"].kind == "f"
        assert out.dtype["score"].kind == "f"


# ─── search_species_pair (synthetic tiny inputs) ────────────────────────


def _write_fasta(path, records):
    path.write_text("".join(f">{h}\n{s}\n" for h, s in records))


class TestSearchSpeciesPair:
    @pytest.fixture
    def tiny_species(self, tmp_path):
        # Three artificially-similar proteins (long enough to clear default
        # band width / k-mer thresholds), with one distinct one mixed in
        seq_a = "MKVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"
        seq_b = "MKVTSPADKTNVKAAWAKVGGHAEEYAAEALERMFISFPTTKTYFPHFDLSHGSAQIKGHGKKVADGLTLAVAHLDDLPGALSGLSDLHAYKLRVDPVNFKLLSHCLLVTLAAHHPDDFTPAVHASLDKFLAQVSAVLTSKYR"
        seq_c = "VHLTPEEKAAVTGLWGKVNVEEVGGEALGRLLVVYPWTQRYFDSFGDLSSPDAVMNNPKVKAHGKKVLGAFSDGLAHLDNLKGTFAQLSDLHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
        _write_fasta(tmp_path / "q.fa", [("q1", seq_a), ("q2", seq_b)])
        _write_fasta(tmp_path / "t.fa", [("t1", seq_a), ("t2", seq_c)])
        q = SpeciesSequences.from_fasta(str(tmp_path / "q.fa"), "q.fa")
        t = SpeciesSequences.from_fasta(str(tmp_path / "t.fa"), "t.fa")
        return q, t

    def test_returns_species_pair_results(self, tiny_species):
        q, t = tiny_species
        r = search_species_pair(
            q, t, matrix_name="BLOSUM62",
            kmer_k=4, max_candidates_per_query=20,
        )
        assert isinstance(r, SpeciesPairResults)
        assert r.query_species == "q.fa"
        assert r.target_species == "t.fa"

    def test_zero_query_sequences_short_circuits(self, tmp_path):
        # Build an empty query, non-empty target
        _write_fasta(tmp_path / "empty.fa", [])
        _write_fasta(tmp_path / "t.fa", [("t1", "ACDEFGHIKLMNPQRSTVWY")])
        q = SpeciesSequences.from_fasta(str(tmp_path / "empty.fa"), "empty.fa")
        t = SpeciesSequences.from_fasta(str(tmp_path / "t.fa"), "t.fa")
        r = search_species_pair(q, t, "BLOSUM62")
        assert len(r.scores) == 0
        assert len(r.target_names) == 0

    def test_score_count_matches_hit_count(self, tiny_species):
        q, t = tiny_species
        r = search_species_pair(
            q, t, matrix_name="BLOSUM62",
            kmer_k=4, max_candidates_per_query=20,
        )
        # Same number of entries across every per-pair array
        n = len(r.scores)
        assert len(r.evalues) == n
        assert len(r.query_names) == n
        assert len(r.target_names) == n


# Note: execute_builtin_search uses multiprocessing.Pool which is flaky
# inside pytest (interactions with numba JIT cache + child-process spawn).
# Its behavior is exercised by the integration tests in tests/integration/.
