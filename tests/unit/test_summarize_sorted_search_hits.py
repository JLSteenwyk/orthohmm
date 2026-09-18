import json

import numpy as np
import pytest

from benchmark_tools import compare_search_hit_coverage as original
from benchmark_tools import summarize_sorted_search_hits as module
from orthohmm.accuracy import load_accuracy_checkpoint, write_accuracy_checkpoint


def arrays(codes, n):
    codes = np.asarray(codes, dtype=np.int64)
    return codes // n, codes % n, np.arange(1, len(codes) + 1, dtype=float)


@pytest.mark.parametrize("seed", range(5))
@pytest.mark.parametrize("chunk", [1, 7, 1000])
def test_exact_overlap_and_coverage_match_existing_implementation(seed, chunk):
    rng = np.random.default_rng(seed)
    n = 41
    species = rng.integers(0, 5, n) * 3
    a = np.sort(rng.choice(n * n, 513, replace=False))
    b = np.sort(rng.choice(n * n, 701, replace=False))
    assert module.overlap(arrays(a, n), arrays(b, n), n, chunk) == original.overlap(a, b, n)
    report = module.summarize(*arrays(a, n), species, chunk)
    expected, _ = original.summarize(*arrays(a, n), species)
    for key in expected:
        if key in report:
            assert report[key] == expected[key], key
    assert sum(report["normalized_score_histogram"]["counts"]) == len(a)
    assert report["normalized_score_min"] == 1
    assert report["normalized_score_max"] == len(a)
    json.dumps(report, allow_nan=False)


@pytest.mark.parametrize("a,b", [([], []), ([], [1, 2]), ([0, 1], []),
    ([0, 1], [2, 3]), ([0, 3], [0, 3])])
def test_empty_disjoint_and_self_hits(a, b):
    assert module.overlap(arrays(a, 2), arrays(b, 2), 2, 1) == original.overlap(
        np.asarray(a, dtype=np.int64), np.asarray(b, dtype=np.int64), 2)


@pytest.mark.parametrize("q,t,s", [([0, 0], [1, 1], [1., 2.]), ([1, 0], [0, 1], [1., 2.]),
    ([2], [0], [1.]), ([-1], [0], [1.]), ([.5], [0], [1.]),
    ([0], [0], [float("nan")]), ([0], [0], [0.]), ([0], [0, 1], [1.])])
def test_invalid_hits_rejected_even_after_other_stream_exhausted(q, t, s):
    invalid = tuple(np.asarray(x) for x in (q, t, s))
    with pytest.raises(ValueError):
        module.overlap(arrays([], 2), invalid, 2, 1)


def test_score_bin_edges_and_empty_summary():
    q = np.zeros(4, dtype=int)
    t = np.arange(4)
    report = module.summarize(q, t, np.array([.01, .03, 100., 200.]), np.arange(4), 2)
    assert report["normalized_score_histogram"]["counts"] == [0, 1, 1, 0, 0, 0, 0, 0, 0, 2]
    empty = module.summarize(*arrays([], 2), np.array([0, 1]))
    assert empty["normalized_score_min"] is None
    assert empty["queries_without_hits"] == 2


@pytest.mark.parametrize("chunk", [0, -1, True, 1.5])
def test_invalid_chunk_size(chunk):
    with pytest.raises(ValueError):
        module.overlap(arrays([], 2), arrays([], 2), 2, chunk)


def test_real_readonly_checkpoint_arrays(tmp_path):
    values = arrays([0, 1, 3, 5, 8], 3)
    path = write_accuracy_checkpoint(str(tmp_path), ["a", "b", "c"], [0, 0, 1], *values)
    names, species, q, t, s = load_accuracy_checkpoint(path, verify=True)
    assert names == ["a", "b", "c"]
    assert all(isinstance(a, np.memmap) and not a.flags.writeable for a in (species, q, t, s))
    report = module.summarize(q, t, s, species, chunk_size=2)
    assert report["directed_hits"] == 5
    assert report["self_hits"] == 2
    assert report["cross_species_hits"] == 1
    assert module.overlap((q, t, s), values, 3, chunk_size=2)["all"]["jaccard"] == 1


@pytest.mark.parametrize("species", [[], [-1, 0], [.5, 1.], [[0, 1]]])
def test_invalid_species_universe(species):
    with pytest.raises(ValueError):
        module.summarize(*arrays([], 2), np.asarray(species))
