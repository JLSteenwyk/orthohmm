from types import SimpleNamespace

import numpy as np
import pytest

from benchmark_tools.search_decision_trace import classify_search_result


def result():
    return SimpleNamespace(query_indices=np.array([0, 1]),
                           target_indices=np.array([1, 0]),
                           scores=np.array([12., -2.]),
                           evalues=np.array([0., 1e-4]), candidate_count=2)


def classify(value, pairs=(("a", "y"), ("b", "x"), ("a", "x"))):
    return classify_search_result(value, ["a", "b"], ["x", "y"], pairs, 1e-4)


def test_three_decisions_and_strict_threshold():
    rows = classify(result())
    assert [r["decision"] for r in rows] == [
        "accepted", "scored_not_significant", "not_selected_by_prefilter"]
    assert rows[1]["score"] == -2.
    assert rows[2]["score"] is None and rows[2]["evalue"] is None


@pytest.mark.parametrize("attribute,value", [
    ("candidate_count", 3), ("scores", np.array([1.])),
    ("query_indices", np.array([-1, 1])),
    ("target_indices", np.array([2, 0])),
    ("query_indices", np.array([.5, 1.])),
    ("scores", np.array([np.nan, 1.])),
    ("scores", np.array([np.inf, 1.])),
    ("evalues", np.array([-1., 0.])),
    ("evalues", np.array([np.nan, 0.])),
    ("evalues", np.array([np.inf, 0.])),
])
def test_invalid_evidence_is_not_a_biological_rejection(attribute, value):
    data = result()
    setattr(data, attribute, value)
    with pytest.raises(ValueError):
        classify(data)


def test_duplicate_candidate():
    data = result()
    data.query_indices[:] = 0
    data.target_indices[:] = 1
    with pytest.raises(ValueError, match="Duplicate directed"):
        classify(data)


@pytest.mark.parametrize("pairs", [[("a", "y"), ("a", "y")], [("missing", "x")]])
def test_invalid_watched_pairs(pairs):
    with pytest.raises(ValueError):
        classify(result(), pairs)


@pytest.mark.parametrize("threshold", [0., -1., np.nan, np.inf])
def test_invalid_threshold(threshold):
    with pytest.raises(ValueError):
        classify_search_result(result(), ["a", "b"], ["x", "y"], [], threshold)


def test_duplicate_ids():
    with pytest.raises(ValueError, match="Duplicate sequence"):
        classify_search_result(result(), ["a", "a"], ["x", "y"], [], 1e-4)


def test_empty_prefilter_result():
    data = result()
    for name in ("query_indices", "target_indices", "scores", "evalues"):
        setattr(data, name, getattr(data, name)[:0])
    data.candidate_count = 0
    assert all(row["decision"] == "not_selected_by_prefilter" for row in classify(data))


def test_actual_cpu_search_agrees_with_engine_filter(tmp_path, monkeypatch):
    from orthohmm.search import engine
    from orthohmm.search.sequences import SpeciesSequences

    monkeypatch.setattr(engine, "is_cuda_available", lambda: False)
    sequence = "ACDEFGHIKLMNPQRSTVWY" * 4
    query_path, target_path = tmp_path / "q.fa", tmp_path / "t.fa"
    query_path.write_text(f">q\n{sequence}\n")
    target_path.write_text(f">t\n{sequence}\n>other\n{'A' * 80}\n")
    query = SpeciesSequences.from_fasta(str(query_path), "q.fa")
    target = SpeciesSequences.from_fasta(str(target_path), "t.fa")
    data = engine.search_species_pair_indexed(
        query, target, "BLOSUM62", max_candidates_per_query=100, n_threads=1)
    assert data.candidate_count > 0
    rows = classify_search_result(data, query.ids, target.ids,
                                  [("q", "t"), ("q", "other")], 1e-4)
    filtered = engine._filter_significant_indexed(data, 1e-4)
    accepted = {(query.ids[q], target.ids[t])
                for q, t in zip(filtered.query_indices, filtered.target_indices)}
    assert {(r["query"], r["target"]) for r in rows if r["decision"] == "accepted"} == accepted
