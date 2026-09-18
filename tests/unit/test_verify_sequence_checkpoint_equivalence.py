import sqlite3
import random

import numpy as np
import pytest

from benchmark_tools import verify_sequence_checkpoint_equivalence as module
from benchmark_tools.convert_qfo_sequence_search_control import convert_hits


@pytest.mark.parametrize("seed", range(5))
def test_streamed_top100_matches_independent_group_selection(seed):
    rng = random.Random(seed)
    rows = [(q, 125 * s + t, str(s), float(rng.randrange(1, 8)), .5)
            for q in range(4) for s in range(3) for t in range(125)]
    rng.shuffle(rows)
    expected = []
    for query in range(4):
        for species in range(3):
            group = [r for r in rows if r[0] == query and r[2] == str(species)]
            selected = sorted(group, key=lambda r: (-r[3], r[1]))[:100]
            expected.extend((r[0], r[1], r[4]) for r in selected)
    with sqlite3.connect(":memory:") as db:
        db.execute("CREATE TABLE expected(q INTEGER,t INTEGER,species TEXT,raw REAL,score REAL)")
        db.executemany("INSERT INTO expected VALUES(?,?,?,?,?)", rows)
        assert list(module.expected_top100(db)) == sorted(expected)


def test_independent_reconstruction_matches_producer_variants(tmp_path):
    metadata = {"q": {"length": 10, "species": "a.fa"},
                **{f"t{i:03d}": {"length": 40, "species": "b.fa"} for i in range(105)}}
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    a.write_text("q\tq\t10\t10\t20\t10\t0\n")
    b.write_text("".join(f"q\tt{i:03d}\t10\t40\t40\t10\t1e-10\n" for i in reversed(range(105))))
    plan = {"searches": [{"output": str(a), "target_fasta": {"path": "/a.fa"}},
                         {"output": str(b), "target_fasta": {"path": "/b.fa"}}]}
    output = tmp_path / "converted"
    output.mkdir()
    variants, _ = convert_hits(plan, metadata, output)
    with sqlite3.connect(":memory:") as db:
        assert module.reconstruct(db, [(a, "a.fa"), (b, "b.fa")], metadata, sorted(metadata)) == 106
        assert module.verify(db, variants["all_hits"]["checkpoint"], metadata, None)["hits"] == 106
        assert module.verify(db, variants["top100"]["checkpoint"], metadata, 100)["hits"] == 101
        assert list(module.expected_top100(db))[-1][1] == sorted(metadata).index("t099")


@pytest.mark.parametrize("problem", ["duplicate", "id", "length", "owner", "nan", "significance", "fields", "empty"])
def test_reconstruction_rejects_bad_source(tmp_path, problem):
    line = "q\tq\t10\t10\t20\t10\t0\n"
    owner = "a.fa"
    if problem == "duplicate":
        line *= 2
    elif problem == "id":
        line = line.replace("q", "unknown", 1)
    elif problem == "length":
        line = line.replace("10", "11", 1)
    elif problem == "owner":
        owner = "wrong"
    elif problem == "nan":
        line = line.replace("20", "nan")
    elif problem == "significance":
        line = line.rstrip().rsplit("\t", 1)[0] + "\t1\n"
    elif problem == "fields":
        line += "bad\n"
    else:
        line = ""
    path = tmp_path / "hits"
    path.write_text(line)
    with sqlite3.connect(":memory:") as db:
        with pytest.raises((ValueError, sqlite3.IntegrityError)):
            module.reconstruct(db, [(path, owner)], {"q": {"length": 10, "species": "a.fa"}}, ["q"])


@pytest.mark.parametrize("problem", [None, "query", "target", "score", "missing", "extra", "dtype"])
def test_tuple_comparison_across_chunks(problem):
    rows = [(0, 0, 2.), (0, 1, 1.), (1, 0, .5)]
    queries, targets = np.array([0, 0, 1], dtype=np.int32), np.array([0, 1, 0], dtype=np.int32)
    scores = np.array([2., 1., .5], dtype=np.float64)
    if problem == "query":
        queries[-1] = 0
    elif problem == "target":
        targets[-1] = 1
    elif problem == "score":
        scores[-1] = .6
    elif problem == "missing":
        queries, targets, scores = queries[:-1], targets[:-1], scores[:-1]
    elif problem == "extra":
        rows.pop()
    elif problem == "dtype":
        queries = queries.astype(np.int64)
    if problem:
        with pytest.raises(ValueError):
            module.compare_rows(rows, queries, targets, scores, chunk_size=2)
    else:
        assert module.compare_rows(rows, queries, targets, scores, chunk_size=2) == 3
