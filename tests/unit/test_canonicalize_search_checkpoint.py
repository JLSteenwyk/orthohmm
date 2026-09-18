import json
import sqlite3

import numpy as np
import pytest

from benchmark_tools import canonicalize_search_checkpoint as module
from orthohmm.accuracy import load_accuracy_checkpoint, write_accuracy_checkpoint


@pytest.mark.parametrize("seed", range(5))
def test_unsorted_checkpoint_exact_copy(tmp_path, seed):
    rng = np.random.default_rng(seed)
    n = 31
    codes = rng.choice(n * n, 311, replace=False)
    q, t = (codes // n).astype(np.int32), (codes % n).astype(np.int32)
    scores = rng.uniform(.001, 100., len(q))
    names = [f"gene{i:03d}" for i in range(n)]
    species = np.arange(n, dtype=np.int32) % 3
    source = write_accuracy_checkpoint(str(tmp_path / "source"), names, species, q, t, scores)
    original = {p.name: module.record(p) for p in source.iterdir()}
    report = module.canonicalize(source, original["manifest.json"]["sha256"], tmp_path / "copy")
    assert report["status"] == "canonical_search_checkpoint_written"
    result_names, result_species, rq, rt, rs = load_accuracy_checkpoint(report["checkpoint"])
    order = np.argsort(codes)
    assert result_names == names
    assert np.array_equal(result_species, species)
    assert np.array_equal(rq, q[order])
    assert np.array_equal(rt, t[order])
    assert np.array_equal(rs, scores[order])
    assert {p.name: module.record(p) for p in source.iterdir()} == original


def test_duplicate_fails_without_collapsing_or_mutating_source(tmp_path):
    source = write_accuracy_checkpoint(str(tmp_path / "source"), ["a", "b"], [0, 1], [0, 0], [1, 1], [1., 2.])
    before = module.record(source / "manifest.json")
    with pytest.raises(sqlite3.IntegrityError):
        module.canonicalize(source, before["sha256"], tmp_path / "copy")
    report = json.loads((tmp_path / "copy/manifest.json").read_text())
    assert report["status"] == "failed"
    assert not (tmp_path / "copy/orthohmm_working_res").exists()
    assert module.record(source / "manifest.json") == before


@pytest.mark.parametrize("q,t,s", [([-1], [0], [1.]), ([2], [0], [1.]),
    ([0], [1], [float("nan")]), ([0], [1], [0.])])
def test_invalid_tuple(q, t, s):
    with sqlite3.connect(":memory:") as database, pytest.raises(ValueError):
        module.ingest(database, np.asarray(q, dtype=np.int32), np.asarray(t, dtype=np.int32),
                      np.asarray(s, dtype=np.float64), 2, chunk_size=1)


def test_empty_checkpoint(tmp_path):
    source = write_accuracy_checkpoint(str(tmp_path / "source"), ["a"], [0], [], [], [])
    report = module.canonicalize(source, module.record(source / "manifest.json")["sha256"], tmp_path / "copy")
    assert report["hits"] == 0


def test_wrong_hash_creates_no_output(tmp_path):
    source = write_accuracy_checkpoint(str(tmp_path / "source"), ["a"], [0], [0], [0], [1.])
    with pytest.raises(ValueError):
        module.canonicalize(source, "0" * 64, tmp_path / "copy")
    assert not (tmp_path / "copy").exists()
