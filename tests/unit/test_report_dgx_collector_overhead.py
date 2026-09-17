import pytest

from benchmark_tools.report_dgx_collector_overhead import MODES, paired_summary


def rows():
    return [{"index": i, "mode": mode, "wall_s": 100., "cpu_s": 1900.} for i, mode in enumerate(MODES)]


def test_counterbalanced_pairs_and_budget():
    data = rows()
    data[1]["wall_s"] = 101
    data[2]["wall_s"] = 102
    data[5]["wall_s"] = 99
    result = paired_summary(data)
    assert result["median_wall_inflation"] == pytest.approx(0.01)
    assert result["range_wall_inflation"] == pytest.approx([-0.01, 0.02])
    assert result["numerical_wall_budget_met"]


def test_unfavorable_pair_not_removed():
    data = rows()
    data[2]["wall_s"] = 111
    result = paired_summary(data)
    assert result["median_wall_inflation"] == 0
    assert not result["numerical_wall_budget_met"]
    assert len(result["pairs"]) == 3


def test_incomplete_inventory_rejected():
    with pytest.raises(ValueError):
        paired_summary(rows()[:-1])
