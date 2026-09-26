import csv
from copy import deepcopy
from types import SimpleNamespace

import pytest

from benchmark_tools.search_sensitivity_hmm import input_inventory, run, write_hits


def result():
    ids = {"a.fasta": ["a0", "a1"], "b.fasta": ["b0"]}
    pairs = {}
    for q in ids:
        for t in ids:
            pairs[q, t] = SimpleNamespace(query_species=q, target_species=t,
                query_indices=[0], target_indices=[0], scores=[23.125],
                evalues=[0.0], candidate_count=2)
    return SimpleNamespace(species_ids=deepcopy(ids), pair_results=pairs), ids


def test_roundtrip_preserves_all_directed_pairs(tmp_path):
    native, ids = result()
    counts = write_hits(native, ids, tmp_path / "hits.tsv")
    with (tmp_path / "hits.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == 4
    assert sum(c["candidates"] for c in counts) == 8
    assert [r["query_id"] for r in rows] == ["a0", "a0", "b0", "b0"]
    assert all(float(r["score"]) == 23.125 and float(r["evalue"]) == 0 for r in rows)
    with pytest.raises(FileExistsError):
        write_hits(native, ids, tmp_path / "hits.tsv")


@pytest.mark.parametrize("field,value", [
    ("query_indices", [-1]), ("target_indices", [99]), ("query_indices", [0.5]),
    ("scores", []), ("scores", [float("nan")]), ("evalues", [-1]),
    ("evalues", [1e-4]), ("evalues", [float("inf")]), ("candidate_count", 0),
    ("candidate_count", 1.5), ("query_species", "wrong"),
])
def test_reject_invalid_native_hits(tmp_path, field, value):
    native, ids = result()
    setattr(native.pair_results["a.fasta", "b.fasta"], field, value)
    with pytest.raises(ValueError):
        write_hits(native, ids, tmp_path / "hits.tsv")


def test_missing_pair_and_changed_id_order(tmp_path):
    native, ids = result()
    del native.pair_results["a.fasta", "b.fasta"]
    with pytest.raises(ValueError, match="Incomplete"):
        write_hits(native, ids, tmp_path / "hits.tsv")
    native, ids = result()
    native.species_ids["a.fasta"].reverse()
    with pytest.raises(ValueError, match="gene IDs"):
        write_hits(native, ids, tmp_path / "hits.tsv")


def test_duplicate_hits(tmp_path):
    native, ids = result()
    pair = native.pair_results["a.fasta", "b.fasta"]
    for field in ("query_indices", "target_indices", "scores", "evalues"):
        setattr(pair, field, getattr(pair, field) * 2)
    with pytest.raises(ValueError, match="Repeated"):
        write_hits(native, ids, tmp_path / "hits.tsv")


def test_empty_hits_valid(tmp_path):
    native, ids = result()
    pair = native.pair_results["a.fasta", "b.fasta"]
    for field in ("query_indices", "target_indices", "scores", "evalues"):
        setattr(pair, field, [])
    assert write_hits(native, ids, tmp_path / "hits.tsv")[1]["hits"] == 0


@pytest.mark.parametrize("suffix", [".fasta", ".fa", ".faa"])
def test_fasta_inventory(tmp_path, suffix):
    (tmp_path / ("a" + suffix)).write_text(">a description\nACD\nEF\n>b\nGH\n")
    records, ids = input_inventory(tmp_path)
    assert ids == {"a" + suffix: ["a", "b"]}
    assert len(records) == 1 and len(records[0]["sha256"]) == 64


@pytest.mark.parametrize("content", ["", ">a\n", ">\nACD\n", "ACD\n", ">a\nAC\n>a\nCC\n", ">a\n>b\nCC\n", ">a\nA C\n"])
def test_bad_fasta(tmp_path, content):
    (tmp_path / "a.fasta").write_text(content)
    with pytest.raises(ValueError):
        input_inventory(tmp_path)


def test_reject_extra_files_and_symlinks(tmp_path):
    path = tmp_path / "truth.json"
    path.write_text("{}")
    with pytest.raises(ValueError):
        input_inventory(tmp_path)
    path.unlink()
    (tmp_path / "a.fasta").symlink_to("missing")
    with pytest.raises(ValueError):
        input_inventory(tmp_path)


def test_require_isolated_interpreter(tmp_path):
    with pytest.raises(RuntimeError, match="-I"):
        run(tmp_path, tmp_path / "output")
