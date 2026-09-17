import gzip
import sqlite3

import pytest

from benchmark_tools.audit_qfo_vgnc_mapping import mapped_reference, reference_data, validate_raw


def compressed(tmp_path, name, text):
    path = tmp_path / name
    with gzip.open(path, "wt") as stream:
        stream.write(text)
    return path


def test_reference_last_label_and_duplicate_rejection(tmp_path):
    path = compressed(tmp_path, "ref.gz", "1\t2\tA\n1\t3\tB\n")
    truth, labels = reference_data(path)
    assert labels == {1: "B", 2: "A", 3: "B"}
    assert len(truth) == 2
    path = compressed(tmp_path, "dup.gz", "1\t2\tA\n2\t1\tA\n")
    with pytest.raises(ValueError, match="Duplicate"):
        reference_data(path)


def test_mapping_is_bijective_and_complete(tmp_path):
    path = tmp_path / "db.sqlite"
    with sqlite3.connect(path) as db:
        db.execute("CREATE TABLE proteomes (prot_nr INT, uniprot_id TEXT, species TEXT)")
        db.executemany("INSERT INTO proteomes VALUES (?, ?, ?)", [(1, "a", "s1"), (2, "b", "s2")])
    truth = {(1, 2): "A"}
    pairs, annotations, digest, duplicates = mapped_reference(path, truth, {1: "A", 2: "A"})
    assert pairs == {("a", "b")}
    assert annotations == {"a": ("A", "s1"), "b": ("A", "s2")}
    assert len(digest) == 64
    assert duplicates == {"identical": 0, "alias": 0}
    with sqlite3.connect(path) as db:
        db.execute("INSERT INTO proteomes VALUES (1, 'a', 's1')")
    assert mapped_reference(path, truth, {1: "A", 2: "A"})[3]["identical"] == 1
    with pytest.raises(ValueError, match="Missing"):
        mapped_reference(path, truth, {1: "A", 2: "A", 3: "B"})
    with sqlite3.connect(path) as db:
        db.execute("UPDATE proteomes SET uniprot_id = 'a'")
    with pytest.raises(ValueError, match="Non-bijective"):
        mapped_reference(path, truth, {1: "A", 2: "A"})


def test_cross_label_true_positive_is_retained(tmp_path):
    path = compressed(tmp_path, "raw.gz", "a\tb\tTP\tA\tB\ts1\ts2\n")
    result = validate_raw(path, {("a", "b")}, {"a": ("A", "s1"), "b": ("B", "s2")})
    assert result["counts"] == {"TP": 1, "FP": 0, "FN": 0}


def test_last_alias_and_conflicting_species(tmp_path):
    path = tmp_path / "aliases.sqlite"
    with sqlite3.connect(path) as db:
        db.execute("CREATE TABLE proteomes (prot_nr INT, uniprot_id TEXT, species TEXT)")
        db.executemany("INSERT INTO proteomes VALUES (?, ?, ?)",
                       [(1, "a", "s1"), (1, "a_s1", "s1"), (2, "b", "s2")])
    pairs, annotations, _, extras = mapped_reference(path, {(1, 2): "A"}, {1: "A", 2: "A"})
    assert pairs == {("a_s1", "b")}
    assert "a" not in annotations
    assert extras == {"identical": 0, "alias": 1}
    with sqlite3.connect(path) as db:
        db.execute("INSERT INTO proteomes VALUES (1, 'other', 's3')")
    with pytest.raises(ValueError, match="Conflicting species"):
        mapped_reference(path, {(1, 2): "A"}, {1: "A", 2: "A"})


@pytest.mark.parametrize("rows, error", [
    ("a\tb\tTP\tX\tA\ts1\ts2\n", "annotation"),
    ("a\tb\tFP\tA\tA\ts1\ts2\n", "eligibility"),
    ("", "partition"),
    ("a\tb\tTP\tA\tA\ts1\ts2\na\tb\tTP\tA\tA\ts1\ts2\n", "duplicate"),
    ("a\tb\tTP\tA\tA\ts1\ts2\na\tb\tFN\tA\tA\ts1\ts2\n", "partition"),
])
def test_invalid_raw_rejected(tmp_path, rows, error):
    path = compressed(tmp_path, "raw.gz", rows)
    with pytest.raises(ValueError, match=error):
        validate_raw(path, {("a", "b")}, {"a": ("A", "s1"), "b": ("A", "s2")})


def test_eligible_cross_family_false_positive(tmp_path):
    annotations = {"a": ("A", "s1"), "b": ("B", "s2"),
                   "c": ("A", "s2"), "d": ("B", "s1")}
    path = compressed(tmp_path, "raw.gz", "a\tc\tFN\tA\tA\ts1\ts2\na\tb\tFP\tA\tB\ts1\ts2\n")
    assert validate_raw(path, {("a", "c")}, annotations)["counts"]["FP"] == 1
