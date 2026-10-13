from types import SimpleNamespace

import numpy as np
import pytest

from benchmark_tools.probe_matched_live_k4_search import search_pair


@pytest.mark.parametrize("instrumented", [False, True])
@pytest.mark.parametrize("fails", [False, True])
def test_matched_diagnostic_uses_actual_k4_and_preserves_calls_state_and_failures(instrumented, fails):
    query, target = object(), object()
    result = SimpleNamespace(evalues=np.array([1e-6, 1e-2]))
    def native(q, t, matrix, **kwargs):
        assert q is query and t is target and matrix == "BLOSUM62"
        assert kwargs == dict(band_width=64, kmer_k=4, max_candidates_per_query=400,
                              n_threads=4, prefilter_phase_a_topk=5000)
        if fails:
            raise ValueError("native failed")
        return result
    engine = SimpleNamespace(search_species_pair_indexed=native)
    def wrapper(function, costs, name):
        def wrapped(*args, **kwargs):
            costs[name] = {"calls": 1}
            return function(*args, **kwargs)
        return wrapped
    helpers = SimpleNamespace(STAGES=("search_species_pair_indexed",), timed_wrapper=wrapper,
                              fingerprint_result=lambda value: value is result)
    if fails:
        with pytest.raises(ValueError, match="native failed"):
            search_pair(engine, helpers, query, target, 400, instrumented)
    else:
        actual = search_pair(engine, helpers, query, target, 400, instrumented)
        assert actual["significant_hits"] == 1 and actual["fingerprint"]
        assert bool(actual["inclusive_stages"]) == instrumented
    assert engine.search_species_pair_indexed is native
