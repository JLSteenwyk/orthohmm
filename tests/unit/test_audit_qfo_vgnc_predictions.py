import gzip
import sqlite3

import pytest

from benchmark_tools.audit_qfo_vgnc_predictions import classify, compare_raw, read_predictions


def test_truth_fp_and_unscored_predictions():
    labels = {1: "A", 2: "A", 3: "B", 4: "B", 5: "C"}
    species = {1: "s1", 2: "s2", 3: "s1", 4: "s2", 5: "s2"}
    predictions = {(1, 2), (1, 4), (1, 5)}
    result = classify({(1, 2), (3, 4)}, predictions, labels, species)
    assert result == {"TP": {(1, 2)}, "FN": {(3, 4)}, "FP": {(1, 4)}}
    assert predictions - result["TP"] - result["FP"] == {(1, 5)}


def test_native_truth_and_fp_can_overlap():
    result = classify({(1, 4)}, {(1, 4)},
                      {1: "A", 2: "A", 3: "B", 4: "B"},
                      {1: "s1", 2: "s2", 3: "s1", 4: "s2"})
    assert result["TP"] == result["FP"] == {(1, 4)}


def test_query_matches_native_orientation_subset_and_set_semantics(tmp_path):
    database = tmp_path / "predictions.sqlite"
    with sqlite3.connect(database) as db:
        db.execute("CREATE TABLE proteomes (prot_nr INT, uniprot_id TEXT, species TEXT)")
        db.execute("CREATE TABLE orthologs (prot_nr1 INT, prot_nr2 INT)")
        db.executemany("INSERT INTO proteomes VALUES (?, ?, ?)",
                       [(1, "a", "s1"), (1, "alias", "s1"), (2, "b", "s2"), (3, "c", "s3")])
        db.executemany("INSERT INTO orthologs VALUES (?, ?)",
                       [(1, 2), (1, 2), (2, 1), (1, 1), (1, 99), (3, 1)])
    metadata, predictions = read_predictions(database, {1: "A", 2: "A", 3: "A"})
    assert metadata[1] == ("alias", "s1")
    assert predictions == {(1, 2)}


@pytest.mark.parametrize("rows, expected_error", [
    ("a\tb\tTP\tA\tA\ts1\ts2\n", None),
    ("", "missing"),
    ("a\tb\tFN\tA\tA\ts1\ts2\n", "mismatch"),
    ("a\tb\tTP\tA\tA\ts1\ts2\na\tb\tTP\tA\tA\ts1\ts2\n", "duplicate"),
])
def test_exact_raw_pair_comparison(tmp_path, rows, expected_error):
    raw = tmp_path / "raw.gz"
    with gzip.open(raw, "wt") as stream:
        stream.write(rows)
    expected = {"TP": {(1, 2)}, "FN": set(), "FP": set()}
    metadata = {1: ("a", "s1"), 2: ("b", "s2")}
    if expected_error:
        with pytest.raises(ValueError, match=expected_error):
            compare_raw(raw, expected, metadata)
    else:
        compare_raw(raw, expected, metadata)


def test_equal_counts_do_not_hide_wrong_false_positive(tmp_path):
    raw = tmp_path / "wrong_fp.gz"
    with gzip.open(raw, "wt") as stream:
        stream.write("a\tc\tFP\tA\tB\ts1\ts2\n")
    expected = {"TP": set(), "FN": set(), "FP": {(1, 2)}}
    with pytest.raises(ValueError, match="1 missing, 1 extra"):
        compare_raw(raw, expected, {1: ("a", "s1"), 2: ("b", "s2"), 3: ("c", "s2")})
