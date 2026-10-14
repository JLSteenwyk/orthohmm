import itertools

import numpy as np
import pytest

from benchmark_tools import benchmark_viterbi_boundary_clear as benchmark


def test_balanced_arms_use_division_only_and_preserve_arrays(monkeypatch):
    libraries = {name: object() for name in benchmark.ARMS}
    calls = []
    def score(data, library, mode):
        assert mode == "division"
        calls.append(library)
        return np.array([4], dtype=np.int32)
    monkeypatch.setattr(benchmark, "score", score)
    data = (np.array([1, 2]),)
    result = benchmark.paired_times(data, libraries)
    assert [r["order"] for r in result["repeats"]] == list(itertools.permutations(benchmark.ARMS))
    assert len(calls) == 21
    assert all(calls.count(lib) == 7 for lib in libraries.values())
    assert result["inputs_unchanged"] and result["exact_scores_match_original"]
    np.testing.assert_array_equal(data[0], [1, 2])


@pytest.mark.parametrize("bad_call", [2, 3, 4, 21])
def test_warm_and_timed_score_mismatches_are_rejected(monkeypatch, bad_call):
    count = 0
    def score(*args):
        nonlocal count
        count += 1
        return np.array([count == bad_call], dtype=np.int32)
    monkeypatch.setattr(benchmark, "score", score)
    with pytest.raises(ValueError, match="scores"):
        benchmark.paired_times((np.array([1]),), {name: object() for name in benchmark.ARMS})


@pytest.mark.parametrize("bad_call", [1, 21])
def test_input_mutation_is_rejected_including_original_warmup(monkeypatch, bad_call):
    count = 0
    def score(data, *args):
        nonlocal count
        count += 1
        if count == bad_call:
            data[0][0] += 1
        return np.array([4], dtype=np.int32)
    monkeypatch.setattr(benchmark, "score", score)
    with pytest.raises(ValueError, match="input arrays"):
        benchmark.paired_times((np.array([1]),), {name: object() for name in benchmark.ARMS})


def test_missing_or_reordered_library_inventory_cannot_change_design(monkeypatch):
    monkeypatch.setattr(benchmark, "score", lambda *args: pytest.fail("Invalid inventory must fail before scoring"))
    with pytest.raises(ValueError, match="prespecified"):
        benchmark.paired_times((), {name: object() for name in reversed(benchmark.ARMS)})
