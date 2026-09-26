from fractions import Fraction

import pytest

from benchmark_tools.score_search_sensitivity import aggregate, read_hits, score, select_cutoff


GENES = {"a": ("s1.fasta", 10), "b": ("s2.fasta", 20), "c": ("s1.fasta", 30), "d": ("s3.fasta", 40)}
FAMILIES = {"a": "f1", "b": "f1", "c": "f2", "d": "f3"}
DENOMINATORS = {"f1": 2, "f2": 0, "f3": 0}


def test_directed_exclusion_and_inclusive_diamond_threshold():
    value = score([("a", "a", 0), ("a", "c", 0), ("a", "b", .1),
                   ("b", "a", .2), ("a", "d", 0)], GENES, FAMILIES, DENOMINATORS, .1)
    assert value["exact_recall"] == "1/2"
    assert value["total_hits"] == 2 and value["homolog_hits"] == 1 and value["nonhomolog_hits"] == 1
    assert value["eligible_families"] == 1


def test_family_weighting_not_pair_pooling():
    genes = {k: (k, 10) for k in "abcdef"}
    families = {k: "f1" if k in "ab" else "f2" for k in genes}
    s = score([("a", "b", 0), ("b", "a", 0)], genes, families, {"f1": 2, "f2": 12}, 1)
    assert s["exact_recall"] == "1/2" and s["homology_denominator"] == 14


def row(split, hmm, diamond):
    def value(number):
        return dict(exact_recall=str(Fraction(number)), total_hits=0, nonhomolog_hits=0,
                    homolog_hits=0, homology_denominator=10)
    return dict(split=split, condition="test", hmm=value(hmm), diamond=[value(x) for x in diamond])


def test_selection_ignores_reporting_and_uses_exact_tie():
    rows = [row("calibration", "1/2", ["1/4", "3/4"]), row("reporting", 1, [0, 1])]
    assert select_cutoff(rows, [1e-20, 1e-4]) == 0
    assert aggregate(rows[:1], 0)["exact_difference"] == "-1/4"


def test_dataset_equal_weighting():
    rows = [row("calibration", 1, [0]), row("calibration", 0, [0])]
    rows[0]["hmm"]["homology_denominator"] = 10000
    assert aggregate(rows, 0)["hmm_recall"] == .5


@pytest.mark.parametrize("text", [
    "a\tb\t10\t20\t30\t40\t0\textra\n", "a\tb\t11\t20\t30\t40\t0\n",
    "unknown\tb\t10\t20\t30\t40\t0\n", "a\tb\t10\t20\tnan\t40\t0\n",
    "a\tb\t10\t20\t30\t40\t-1\n", "a\tb\t10\t20\t30\t40\t2\n",
    "a\tb\t10\t20\t30\t40\t0\na\tb\t10\t20\t30\t40\t0\n",
])
def test_invalid_diamond(tmp_path, text):
    path = tmp_path / "hits"
    path.write_text(text)
    with pytest.raises(ValueError):
        read_hits(path, GENES, "diamond", "s2.fasta")


def test_valid_diamond_and_target_check(tmp_path):
    path = tmp_path / "hits"
    path.write_text("a\tb\t10\t20\t30\t40\t1\n")
    assert read_hits(path, GENES, "diamond", "s2.fasta") == [("a", "b", 1)]
    with pytest.raises(ValueError):
        read_hits(path, GENES, "diamond", "s1.fasta")


def test_hmm_strict_cutoff_and_header(tmp_path):
    path = tmp_path / "hits"
    header = "query_species\ttarget_species\tquery_id\ttarget_id\tscore\tevalue\n"
    path.write_text(header + "s1.fasta\ts2.fasta\ta\tb\t3.5\t0\n")
    assert read_hits(path, GENES, "hmm") == [("a", "b", 0)]
    path.write_text(header + "s1.fasta\ts2.fasta\ta\tb\t3.5\t0.0001\n")
    with pytest.raises(ValueError):
        read_hits(path, GENES, "hmm")
    path.write_text("wrong\n")
    with pytest.raises(ValueError):
        read_hits(path, GENES, "hmm")


def test_zero_hit_file_is_valid(tmp_path):
    path = tmp_path / "hits"
    path.write_text("")
    assert read_hits(path, GENES, "diamond", "s2.fasta") == []


def test_no_calibration_and_no_eligible_families():
    with pytest.raises(ValueError):
        select_cutoff([row("reporting", 0, [0])], [1])
    with pytest.raises(ValueError):
        score([], GENES, FAMILIES, {"f1": 0}, 1)
