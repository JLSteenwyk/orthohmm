import itertools

import numpy as np
import pytest

from benchmark_tools import benchmark_viterbi_center_recurrence as benchmark


@pytest.mark.parametrize("fails", [False, True])
def test_score_restores_native_loader_state(monkeypatch, fails):
    previous, library = object(), object()
    monkeypatch.setattr(benchmark.viterbi, "_c_hmm_lib", previous)
    monkeypatch.setattr(benchmark.viterbi, "_c_hmm_available", None)

    def native(*args, **kwargs):
        assert benchmark.viterbi._c_hmm_lib is library
        assert benchmark.viterbi._c_hmm_available is True
        assert args == (64,)
        assert kwargs == {"n_threads": 4, "band_center_mode": "recurrence"}
        if fails:
            raise RuntimeError("native failed")
        return np.array([4], dtype=np.int32)

    monkeypatch.setattr(benchmark.viterbi, "batch_viterbi_multipair_c", native)
    if fails:
        with pytest.raises(RuntimeError, match="native failed"):
            benchmark.score((), library, "recurrence")
    else:
        np.testing.assert_array_equal(benchmark.score((), library, "recurrence"), [4])
    assert benchmark.viterbi._c_hmm_lib is previous
    assert benchmark.viterbi._c_hmm_available is None


def test_three_arm_schedule_is_balanced_and_original_kernel_is_measured(monkeypatch):
    current, original = object(), object()
    calls = []

    def score(data, library, mode):
        calls.append((library, mode))
        return np.array([7, 9], dtype=np.int32)

    monkeypatch.setattr(benchmark, "score", score)
    result = benchmark.paired_times((np.array([1, 2]),), current, original)
    expected = list(itertools.permutations(("original_division", "division", "recurrence")))
    assert [row["order"] for row in result["repeats"]] == expected
    assert len(calls) == 21
    assert sum(library is original for library, _ in calls) == 7
    assert result["exact_scores_match_live_kernel"] and result["inputs_unchanged"]
    assert result["median_recurrence_over_division_wall"] > 0
    assert result["median_recurrence_over_original_wall"] > 0


@pytest.mark.parametrize("bad_call", [2, 3, 4, 21])
def test_score_mismatches_fail_during_warmup_or_timed_runs(monkeypatch, bad_call):
    count = 0

    def score(*args):
        nonlocal count
        count += 1
        return np.array([2 if count == bad_call else 1], dtype=np.int32)

    monkeypatch.setattr(benchmark, "score", score)
    with pytest.raises(ValueError, match="scores|integer"):
        benchmark.paired_times((np.array([1]),), object(), object())


@pytest.mark.parametrize("bad_call", [1, 21])
def test_input_mutation_including_original_warmup_fails(monkeypatch, bad_call):
    count = 0

    def score(data, *args):
        nonlocal count
        count += 1
        if count == bad_call:
            data[0][0] += 1
        return np.array([1], dtype=np.int32)

    monkeypatch.setattr(benchmark, "score", score)
    with pytest.raises(ValueError, match="modified input"):
        benchmark.paired_times((np.array([1]),), object(), object())


def test_input_array_fingerprints_capture_layout_and_values():
    first = np.array([[1, 2]], dtype=np.int32)
    same_bytes = first.reshape(2)
    assert benchmark.array_records((first,))[0]["sha256"] == benchmark.array_records((same_bytes,))[0]["sha256"]
    assert benchmark.array_records((first,)) != benchmark.array_records((same_bytes,))
    assert benchmark.array_records((first,))[0]["dtype"] == "int32"
