from types import SimpleNamespace

import numpy as np
import pytest

from orthohmm.search import engine, prefilter


def species(n, name="Q"):
    return SimpleNamespace(num_sequences=n, species_file=name, flat_sequences=np.arange(n, dtype=np.int32),
                           offsets=np.arange(n, dtype=np.int64), lengths=np.ones(n, dtype=np.int32))


def index():
    return SimpleNamespace(alpha_size=20, kmer_offsets=np.array([]), kmer_entries=np.array([]),
                           kmer_freqs=np.array([]), freq_threshold=1)


def probe(monkeypatch, n, successful):
    observed = []
    def score(sequence, *args):
        query = int(sequence[0])
        observed.append(query)
        ids = np.array([0] if query in successful else [], dtype=np.int32)
        return ids, np.zeros(len(ids))
    monkeypatch.setattr(prefilter, "_score_query_two_stage", score)
    result = engine._adaptive_max_candidates(species(n), species(1, "T"), 5, 4, False, prepared_index=index())
    return result, observed


@pytest.mark.parametrize("n", [1, 49, 499, 500, 501, 999, 1000, 1001, 1499, 1501, 10001])
def test_coverage_uses_actual_sample_count_and_preserves_sampling_positions(monkeypatch, n):
    expected = list(range(0, n, max(1, n // min(n, 500))))
    result, observed = probe(monkeypatch, n, set(expected))
    assert observed == expected
    assert result == (20, 1.0)


def test_denominator_bias_previously_selected_wrong_candidate_tier(monkeypatch):
    result, observed = probe(monkeypatch, 999, set(range(420)))
    assert len(observed) == 999
    assert result == (100, 420 / 999)
    assert 420 / 500 > .8


@pytest.mark.parametrize("successes,cap", [(0, 200), (100, 200), (101, 100), (250, 100),
                                         (251, 60), (400, 60), (401, 20), (500, 20)])
def test_existing_strict_threshold_boundaries_are_unchanged(monkeypatch, successes, cap):
    result, observed = probe(monkeypatch, 500, set(range(successes)))
    assert observed == list(range(500))
    assert result == (cap, successes / 500)


def test_no_hits_and_partial_hits_use_the_same_actual_denominator(monkeypatch):
    expected = list(range(0, 1501, 3))
    result, observed = probe(monkeypatch, 1501, set(expected[:150]))
    assert observed == expected
    assert result == (100, 150 / len(expected))
    result, observed = probe(monkeypatch, 1501, set())
    assert result == (200, 0.0)


def test_empty_queries_skip_probe_and_index_preparation(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("Empty query collection must not prepare or score a probe")
    monkeypatch.setattr(engine, "prepare_kmer_index", forbidden)
    monkeypatch.setattr(prefilter, "_score_query_two_stage", forbidden)
    assert engine._adaptive_max_candidates(species(0), species(1), 5, 4, True) == (200, 0.0)


@pytest.mark.parametrize("cap", [20, 100, 400])
def test_fixed_caps_never_call_adaptive_probe(monkeypatch, cap):
    def forbidden(*args, **kwargs):
        pytest.fail("Fixed candidate cap called adaptive sampling")
    monkeypatch.setattr(engine, "_adaptive_max_candidates", forbidden)
    observed = {}
    def candidates(*args, **kwargs):
        observed.update(kwargs)
        return np.array([], dtype=np.int32), np.array([0, 0], dtype=np.int64)
    monkeypatch.setattr(engine, "prefilter_candidates", candidates)
    result = engine.search_species_pair_indexed(species(1, "Q"), species(1, "T"), "BLOSUM62",
        max_candidates_per_query=cap, prepared_index=index())
    assert observed["max_candidates_per_query"] == cap
    assert result.candidate_count == 0 and len(result.scores) == 0
