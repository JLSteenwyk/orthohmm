import numpy as np
import pytest

from orthohmm.search import engine
from orthohmm.search.sequences import SpeciesSequences


def species(name, lengths):
    return SpeciesSequences(
        name, [f"{name}_{i}" for i in range(len(lengths))],
        np.zeros(sum(lengths), dtype=np.uint8),
        np.array([0, *np.cumsum(lengths)[:-1]], dtype=np.int64),
        np.array(lengths, dtype=np.int32),
    )


@pytest.mark.parametrize("lengths", [[1998], [1999], [1998, 1999], [1999, 1998]])
@pytest.mark.parametrize("cuda", ["unavailable", "available", "failure"])
@pytest.mark.parametrize("cpu_backend", ["multipair", "scalar", "numba"])
def test_every_candidate_scored_once_by_successful_backend(monkeypatch, lengths, cuda, cpu_backend):
    query, target = species("q", [10]), species("t", lengths)
    count = len(lengths)
    monkeypatch.setattr(engine, "is_cuda_available", lambda: cuda != "unavailable")
    monkeypatch.setattr(engine, "is_c_available", lambda: cpu_backend != "numba")
    monkeypatch.setattr(engine, "prefilter_candidates", lambda *args, **kwargs: (
        np.arange(count, dtype=np.int32), np.array([0, count], dtype=np.int64)))
    monkeypatch.setattr(engine, "build_profiles_batch", lambda *args: (None,) * 5)
    successful = []
    attempted_gpu = []
    attempted_cpu = []

    def cpu(*args, **kwargs):
        pairs = args[8]
        attempted_cpu.extend(pairs[:, 1].tolist())
        successful.extend(pairs[:, 1].tolist())
        return 100 + pairs[:, 1]

    def gpu(*args, **kwargs):
        pairs = args[8]
        attempted_gpu.extend(pairs[:, 1].tolist())
        assert all(lengths[i] <= 1998 for i in pairs[:, 1])
        if cuda == "failure":
            raise RuntimeError("fixture CUDA failure")
        successful.extend(pairs[:, 1].tolist())
        return 100 + pairs[:, 1]

    def evalues(scores, *args):
        np.testing.assert_array_equal(scores, 100 + np.arange(count))
        return np.full(count, 1e-8)

    def multipair(*args, **kwargs):
        if cpu_backend == "scalar":
            raise RuntimeError("fixture multipair unavailable")
        assert cpu_backend == "multipair"
        return cpu(*args, **kwargs)

    monkeypatch.setattr(engine, "batch_viterbi_multipair_c", multipair)
    monkeypatch.setattr(engine, "batch_viterbi_c", cpu)
    monkeypatch.setattr(engine, "batch_viterbi_score", cpu)
    monkeypatch.setattr(engine, "batch_viterbi_cuda", gpu)
    monkeypatch.setattr(engine, "batch_estimate_evalues", evalues)
    result = engine.search_species_pair_indexed(
        query, target, "BLOSUM62", max_candidates_per_query=100, prepared_index=object())
    assert sorted(successful) == list(range(count))
    eligible = [i for i, length in enumerate(lengths) if length <= 1998]
    assert attempted_gpu == (eligible if cuda != "unavailable" else [])
    assert attempted_cpu == (list(range(count)) if cuda in ("unavailable", "failure")
                             else [i for i in range(count) if i not in eligible])
    np.testing.assert_allclose(result.scores, (100 + np.arange(count)) / np.sqrt(10 * np.array(lengths)))
    np.testing.assert_array_equal(result.target_indices, np.arange(count))


def test_long_targets_match_real_cpu_scoring_when_cuda_is_available(monkeypatch):
    query, target = species("q", [30]), species("t", [1999, 2001])
    monkeypatch.setattr(engine, "prefilter_candidates", lambda *args, **kwargs: (
        np.array([0, 1], dtype=np.int32), np.array([0, 2], dtype=np.int64)))

    def unexpected_gpu(*args, **kwargs):
        raise AssertionError("Long-only batch must not call CUDA")

    monkeypatch.setattr(engine, "batch_viterbi_cuda", unexpected_gpu)
    monkeypatch.setattr(engine, "is_cuda_available", lambda: False)
    cpu = engine.search_species_pair_indexed(
        query, target, "BLOSUM62", max_candidates_per_query=100, prepared_index=object())
    monkeypatch.setattr(engine, "is_cuda_available", lambda: True)
    routed = engine.search_species_pair_indexed(
        query, target, "BLOSUM62", max_candidates_per_query=100, prepared_index=object())
    assert np.isfinite(cpu.scores).all() and (cpu.scores > 0).all()
    np.testing.assert_array_equal(cpu.scores, routed.scores)
    np.testing.assert_array_equal(cpu.evalues, routed.evalues)
