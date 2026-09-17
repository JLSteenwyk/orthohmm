import gzip

import pytest

from benchmark_tools.audit_qfo_swiss_counts import HEADER, read_raw, statistics, verify_native


def test_native_prior():
    assert statistics({"TP": 0, "FP": 0, "FN": 0}) == {"PPV": 0.5, "TPR": 0.5, "F1": 0.5}
    assert statistics({"TP": 3, "FP": 1, "FN": 5})["TPR"] == 5 / 12


def test_raw_truth_and_members(tmp_path):
    path = tmp_path / "raw.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(HEADER + "\nA\tx\ty\tTP\nA\tx\tz\tFN\nA\ty\tz\tFP\n")
    counts, truth, members = read_raw(path, ["A"])
    assert counts["A"] == {"TP": 1, "FN": 1, "FP": 1, "TN": 0}
    assert truth["A", "x", "z"] is True
    assert truth["A", "y", "z"] is False
    assert members == {"A": {"x", "y", "z"}}


@pytest.mark.parametrize("rows", ["A\tx\ty\tTP\nA\ty\tx\tTP\n", "A\tx\tx\tTP\n",
                                    "B\tx\ty\tTP\n", "A\tx\ty\tBAD\n", "A\tx\ty\n", ""])
def test_bad_raw(tmp_path, rows):
    path = tmp_path / "raw.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(HEADER + "\n" + rows)
    with pytest.raises(ValueError):
        read_raw(path, ["A"])


def test_macro_statistic_not_mean_family_f1():
    counts = {"A": {"TP": 16, "FP": 0, "FN": 16}, "B": {"TP": 16, "FP": 16, "FN": 0}}
    native = {}
    for family, row in counts.items():
        for metric in ("TPR", "PPV"):
            native["SwissTrees-" + family, metric] = {"metrics": {"value": statistics(row)[metric]}}
    for metric in ("TPR", "PPV"):
        native["SwissTrees", metric] = {"metrics": {"value": 0.7}}
    families, aggregate = verify_native(counts, native)
    assert aggregate["F1"] == pytest.approx(0.7)
    assert aggregate["F1"] != pytest.approx(sum(r["F1"] for r in families.values()) / 2)
    native["SwissTrees", "TPR"]["metrics"]["value"] = 0.6
    with pytest.raises(ValueError, match="aggregate"):
        verify_native(counts, native)
