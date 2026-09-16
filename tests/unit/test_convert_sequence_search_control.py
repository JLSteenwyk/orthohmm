import sqlite3

import numpy as np
import pytest

from benchmark_tools.convert_sequence_search_control import ingest, initialize, parse_hits, selected_sql, write_variant


@pytest.fixture
def metadata():
    return {"a": {"length": 100, "species": "s1"}, "b": {"length": 400, "species": "s2"}}


def test_raw_score_normalized_once_and_direction_preserved(tmp_path, metadata):
    path = tmp_path / "hits"
    path.write_text("a\tb\t100\t400\t200\t80\t1e-10\n")
    assert list(parse_hits(path, metadata, {"a": 0, "b": 1}, "s2", 1)) == [(0, 1, 1, 200, 1)]


@pytest.mark.parametrize("line", [
    "a\tb\t99\t400\t200\t80\t1e-10", "a\tb\t100\t400\tnan\t80\t1e-10",
    "a\tb\t100\t400\t200\t80\t0.1", "x\tb\t100\t400\t200\t80\t1e-10",
    "a\tb\t100\t400\t-1\t80\t1e-10", "a\tb\t100\t400\t200\t80",
])
def test_invalid_hits_rejected(tmp_path, metadata, line):
    path = tmp_path / "hits"
    path.write_text(line + "\n")
    with pytest.raises(ValueError):
        list(parse_hits(path, metadata, {"a": 0, "b": 1}, "s2", 1))


def test_wrong_species_rejected(tmp_path, metadata):
    path = tmp_path / "hits"
    path.write_text("a\tb\t100\t400\t200\t80\t1e-10\n")
    with pytest.raises(ValueError, match="species"):
        list(parse_hits(path, metadata, {"a": 0, "b": 1}, "s1", 0))


def test_duplicate_pairs_fail_instead_of_silent_deduplication():
    with sqlite3.connect(":memory:") as db:
        initialize(db)
        with pytest.raises(sqlite3.IntegrityError):
            ingest(db, [(0, 1, 1, 200, 1), (0, 1, 1, 100, .5)])


def test_top100_is_per_query_target_species_with_stable_ties():
    with sqlite3.connect(":memory:") as db:
        initialize(db)
        ingest(db, [(0, i, 0, 200, 1) for i in reversed(range(101))] + [(0, 101, 1, 200, 1), (1, 102, 0, 200, 1)])
        rows = list(db.execute(selected_sql(100)))
        assert [r[1] for r in rows] == list(range(100)) + [101, 102]


def test_written_checkpoint_is_numeric_audited_and_retains_self_hits(tmp_path):
    with sqlite3.connect(":memory:") as db:
        initialize(db)
        ingest(db, [(0, 0, 0, 200, 2), (0, 1, 1, 100, .5)])
        result = write_variant(db, tmp_path / "variant", ["a", "b"], np.array([0, 1]), None)
    assert result["audit"]["summary"]["hits"] == 2
    assert result["audit"]["summary"]["self_hits"] == 1
    assert result["audit"]["summary"]["genes"] == 2
