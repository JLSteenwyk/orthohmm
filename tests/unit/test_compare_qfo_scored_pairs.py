import gzip

import pytest

from benchmark_tools.compare_qfo_scored_pairs import audit, compare, read_scores


def raw(path, rows, metric="GO"):
    with gzip.open(path, "wt") as stream:
        stream.write(f"# {metric} Similarities between orthologs from fixture\n"
                     f"# Computing timestamp: test\n"
                     f"# Protein ID 1<tab>Protein ID 2<tab>{metric} Similarity\n" + rows)
    return path


@pytest.mark.parametrize("metric", ["GO", "EC"])
def test_exact_serialization_and_direction(tmp_path, metric):
    path = raw(tmp_path/"a.gz", "B\tA\t0.123456\nC\tA\t1.000000\n", metric)
    assert read_scores(path, metric) == {("A", "B"): 123456, ("A", "C"): 1000000}


@pytest.mark.parametrize("rows", ["", "A\tA\t0.100000\n", "A\tB\tnan\n",
    "A\tB\t1.000001\n", "A\tB\t0.1\n", "A\tB\t-0.100000\n",
    "A\tB\t0.100000\tX\n", "A\tB\t0.100000\nB\tA\t0.100000\n"])
def test_bad_rows(tmp_path, rows):
    with pytest.raises(ValueError):
        read_scores(raw(tmp_path/"a.gz", rows), "GO")


def test_wrong_header(tmp_path):
    with pytest.raises(ValueError, match="header"):
        read_scores(raw(tmp_path/"a.gz", "A\tB\t0.100000\n", "EC"), "GO")


def test_decomposition_keeps_original_denominators():
    left = {("A", "B"): 500000, ("A", "C"): 1000000}
    right = {("A", "B"): 400000, ("A", "D"): 200000, ("A", "E"): 300000}
    result = compare(left, right)
    assert result["original_mean_difference"] == pytest.approx(.45)
    assert sum(result["original_mean_difference_components"].values()) == pytest.approx(.45)
    assert result["shared_conditional_mean_difference"] == pytest.approx(.1)
    assert result["shared_pairs_with_different_serialized_scores"] == 1
    assert result["maximum_shared_absolute_difference_millionths"] == 100000


def test_disjoint_and_identical():
    a = {("A", "B"): 500000}
    same = compare(a, a)
    assert same["original_mean_difference"] == 0
    assert same["left_only_pairs"] == same["right_only_pairs"] == 0
    disjoint = compare(a, {("C", "D"): 500000})
    assert disjoint["shared_conditional_mean_difference"] is None
    assert disjoint["maximum_shared_absolute_difference_millionths"] is None


def test_audit_and_no_overwrite(tmp_path):
    left = raw(tmp_path/"a.gz", "A\tB\t0.100000\n")
    right = raw(tmp_path/"b.gz", "B\tA\t0.100000\n")
    output = tmp_path/"result.json"
    result = audit(left, right, "GO", output)
    assert result["uncertainty_admitted"] is False
    assert result["result"]["shared_pairs"] == 1
    with pytest.raises(FileExistsError):
        audit(left, right, "GO", output)
