import json
from pathlib import Path

import numpy as np
import pytest

from benchmark_tools.run_search_recovery import (
    assert_control_fixture, link_changes, resolve_settings,
)
from orthohmm.search import engine, prefilter, sequences


def test_prespecified_arms_change_one_mechanism():
    plan = json.loads((Path(__file__).parents[2] /
        "benchmark_tools/hillclimb_search_recovery_20261012.json").read_text())
    assert resolve_settings(plan, "control") == plan["base_settings"]
    for arm in plan["arms"][1:]:
        settings = resolve_settings(plan, arm["label"])
        assert sum(settings[key] != value for key, value in plan["base_settings"].items()) == 1
    with pytest.raises(ValueError, match="Unknown"):
        resolve_settings(plan, "missing")


def test_duplicate_arms_and_unrecognized_overrides_rejected():
    with pytest.raises(ValueError, match="duplicate"):
        resolve_settings({"arms": [{"label": "a"}, {"label": "a"}]}, "a")
    with pytest.raises(ValueError, match="parameter"):
        resolve_settings({"base_settings": {}, "arms": [{"label": "a", "overrides": {"bad": 1}}]}, "a")


def rows(decisions):
    return [{"query": "q", "target": str(i), "decision": value}
            for i, value in enumerate(decisions)]


def test_link_gains_do_not_hide_control_losses():
    control = rows(["accepted", "scored_not_significant", "not_selected_by_prefilter", "accepted"])
    observed = rows(["not_selected_by_prefilter", "accepted", "not_selected_by_prefilter", "accepted"])
    assert link_changes(control, observed) == {"gained": 1, "lost": 1, "retained": 1, "absent_both": 1}


def test_link_comparison_rejects_missing_duplicate_or_unknown_decisions():
    control = rows(["accepted"])
    for observed in ([], control * 2, rows(["error"])):
        with pytest.raises(ValueError):
            link_changes(control, observed)


@pytest.mark.parametrize("value", [0, -1, 2**31, 1.2, "20", True, np.bool_(True)])
def test_invalid_shortlist_fails_before_native_allocation(value):
    with pytest.raises(ValueError, match="positive int32"):
        prefilter.prefilter_candidates(None, None, phase_a_topk=value)


@pytest.mark.parametrize("available,matrix", [(False, np.zeros((20, 20))), (True, None)])
def test_shortlist_change_cannot_silently_use_different_ranking(monkeypatch, available, matrix):
    monkeypatch.setattr(prefilter, "_load_c_prefilter", lambda: (None, available))
    with pytest.raises(ValueError, match="C extension ranking"):
        prefilter.prefilter_candidates(None, None, sub_matrix=matrix, phase_a_topk=20000)


def test_default_preserves_real_cpu_scoring_and_order(monkeypatch):
    monkeypatch.setattr(engine, "is_cuda_available", lambda: False)
    assert assert_control_fixture(engine, sequences)


def test_shortlist_reaches_native_ranker_and_changes_candidates(tmp_path):
    if not prefilter.is_c_prefilter_available():
        pytest.skip("C extension unavailable")
    path = tmp_path / "seq.fa"
    sequence = "ACDEFGHIKLMNPQRSTVWY" * 2
    path.write_text("".join(f">g{i}\n{sequence}\n" for i in range(3)))
    species = sequences.SpeciesSequences.from_fasta(str(path), "seq.fa")
    from orthohmm.search.matrices import get_matrix
    kwargs = dict(k=4, use_reduced_alphabet=False, min_total_hits=1,
                  min_diag_hits=1, max_candidates_per_query=3, n_threads=1,
                  sub_matrix=get_matrix("BLOSUM62"))
    narrow = prefilter.prefilter_candidates(species, species, phase_a_topk=1, **kwargs)
    broad = prefilter.prefilter_candidates(species, species, phase_a_topk=3, **kwargs)
    np.testing.assert_array_equal(np.diff(narrow[1]), [1, 1, 1])
    np.testing.assert_array_equal(np.diff(broad[1]), [3, 3, 3])
    for q in range(3):
        assert set(narrow[0][narrow[1][q]:narrow[1][q + 1]]) <= set(broad[0][broad[1][q]:broad[1][q + 1]])


def test_engine_passes_shortlist_setting(monkeypatch, tmp_path):
    path = tmp_path / "seq.fa"
    path.write_text(">g\nACDEFGHIKLMNPQRSTVWY\n")
    species = sequences.SpeciesSequences.from_fasta(str(path), "seq.fa")
    calls = []
    def capture(*args, **kwargs):
        calls.append(kwargs["phase_a_topk"])
        return np.empty(0, dtype=np.int32), np.zeros(2, dtype=np.int64)
    monkeypatch.setattr(engine, "prefilter_candidates", capture)
    result = engine.search_species_pair_indexed(species, species, "BLOSUM62",
        max_candidates_per_query=100, prefilter_phase_a_topk=20000)
    assert calls == [20000]
    assert result.candidate_count == 0
