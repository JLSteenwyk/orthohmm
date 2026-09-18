import gzip
import math
import statistics

import pytest

from benchmark_tools.audit_qfo_go_ec import read_raw, validate


def write_raw(tmp_path, rows, metric="GO"):
    path = tmp_path / "raw.gz"
    with gzip.open(path, "wt") as stream:
        stream.write(f"# {metric} Similarities between orthologs from test\n# Computing timestamp: fixture\n"
                     f"# Protein ID 1<tab>Protein ID 2<tab>{metric} Similarity\n" + rows)
    return path


@pytest.mark.parametrize("metric", ["GO", "EC"])
def test_native_direction_need_not_be_accession_sorted(tmp_path, metric):
    result = read_raw(write_raw(tmp_path, "B\tA\t0.200000\nC\tA\t0.800000\n", metric), metric)
    assert result["pairs"] == 2
    assert result["mean_from_rounded_raw"] == 0.5
    assert math.isclose(result["pair_iid_sem_from_rounded_raw"], 0.3)
    assert result["proteins_in_multiple_pairs"] == 1
    assert result["maximum_protein_degree"] == 2
    validated = validate(result, {"metric_x": 2, "metric_y": .5, "stderr_y": 3.81}, 12.7)
    assert math.isclose(validated["native_95_half_width_from_rounded_raw"], 3.81)


@pytest.mark.parametrize("rows", ["A\tA\t0.100000\n", "A\tB\tnan\n", "A\tB\t1.100000\n",
                                  "A\tB\t0.1\n", "A\tB\t0.100000\nB\tA\t0.200000\n",
                                  "A\tB\t0.100000\n", "A\tB\t0.100000\textra\n"])
def test_invalid_raw_rejected(tmp_path, rows):
    with pytest.raises(ValueError):
        read_raw(write_raw(tmp_path, rows), "GO")


def test_wrong_endpoint_header_rejected(tmp_path):
    with pytest.raises(ValueError, match="header"):
        read_raw(write_raw(tmp_path, "A\tB\t0.100000\n", "EC"), "GO")


@pytest.mark.parametrize("key,value", [("metric_x", 3), ("metric_y", .50001),
                                       ("stderr_y", .3), ("stderr_y", float("nan"))])
def test_incorrect_aggregate_and_sem_interpretation_rejected(key, value):
    summary = {"pairs": 2, "mean_from_rounded_raw": .5, "pair_iid_sem_from_rounded_raw": .3}
    participant = {"metric_x": 2, "metric_y": .5, "stderr_y": 3.81}
    participant[key] = value
    with pytest.raises(ValueError):
        validate(summary, participant, 12.7)


def test_six_decimal_rounding_bound(tmp_path):
    full = [.20000049, .79999951, .40000049]
    rows = "".join(f"A\t{b}\t{v:.6f}\n" for b, v in zip("BCD", full))
    summary = read_raw(write_raw(tmp_path, rows), "GO")
    participant = {"metric_x": 3, "metric_y": statistics.mean(full),
                   "stderr_y": 4.3 * statistics.stdev(full) / math.sqrt(3)}
    validate(summary, participant, 4.3)
