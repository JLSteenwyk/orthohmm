import gzip
import math

import pytest

from benchmark_tools.audit_qfo_fas_samples import read_sample, summarize, render_table


def sample(tmp_path, text):
    path = tmp_path / "raw.gz"
    with gzip.open(path, "wt") as stream:
        stream.write(text)
    return path


def test_recomputed_statistics_and_dependent_proteins(tmp_path):
    pairs = read_sample(sample(tmp_path, "Acc1\tAcc2\tFAS\nA\tB\t0.2\nA\tC\t0.8\n"))
    result = summarize(pairs, {"metric_y": 0.5, "stderr_y": 0.3, "metric_x": 100})
    assert result["sample_fraction"] == 0.02
    assert result["sample_proteins"] == 3
    assert result["proteins_in_multiple_sample_pairs"] == 1
    assert result["maximum_sample_protein_degree"] == 2
    assert math.isclose(result["native_pair_iid_sem"], 0.3)


@pytest.mark.parametrize("row", ["B\tA\t0.5", "A\tA\t0.5", "A_X\tB\t0.5",
                                  "A\tB\tnan", "A\tB\tinf", "A\tB\t-0.1",
                                  "A\tB\t1.1", "A\tB", "A\tB\t0.5\textra"])
def test_malformed_rows_rejected(tmp_path, row):
    with pytest.raises(ValueError):
        read_sample(sample(tmp_path, "Acc1\tAcc2\tFAS\n" + row + "\n"))


def test_duplicate_pair_rejected(tmp_path):
    with pytest.raises(ValueError, match="Duplicate"):
        read_sample(sample(tmp_path, "Acc1\tAcc2\tFAS\nA\tB\t0.5\nA\tB\t0.5\n"))


@pytest.mark.parametrize("key,value", [("metric_y", 0.6), ("stderr_y", 0.2),
                                       ("metric_x", 1), ("metric_x", 2.5), ("metric_x", float("nan"))])
def test_aggregate_mutations_rejected(key, value):
    participant = {"metric_y": 0.5, "stderr_y": 0.3, "metric_x": 100}
    participant[key] = value
    with pytest.raises(ValueError):
        summarize({("A", "B"): 0.2, ("A", "C"): 0.8}, participant)


def test_header_and_small_sample_rejected(tmp_path):
    for text in ("", "wrong\n", "Acc1\tAcc2\tFAS\nA\tB\t0.2\n"):
        with pytest.raises(ValueError):
            read_sample(sample(tmp_path, text))


def test_table_uses_raw_derived_statistics():
    row = summarize({("A", "B"): 0.2, ("A", "C"): 0.8},
                    {"metric_y": 0.5, "stderr_y": 0.3, "metric_x": 100})
    assert "| test | 2 | 100 | 2.0000 | 0.500000000 |" in render_table(
        {"methods": [{"method": "test", **row}]})
