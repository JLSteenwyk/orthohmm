from copy import deepcopy

import pytest

from benchmark_tools.readback_boundary_clear_search import ARMS,ORDERS,timing_summary


def case():
    repeats = [{"repeat":index,"order":list(order),"arms":{ARMS[0]:{"wall_s":2.,"process_cpu_s":8.},
        ARMS[1]:{"wall_s":1.8,"process_cpu_s":7.2}}} for index,order in enumerate(ORDERS)]
    return {"repeats":repeats,"ratios":{"wall_s":.9,"process_cpu_s":.9}}


def test_raw_native_search_medians_reproduce_without_using_reported_ratios():
    assert timing_summary(case()) == {"wall_s":.9,"process_cpu_s":.9}


@pytest.mark.parametrize("change", ["ratio","order","missing","bad_cost","bad_repeat","bad_arm"])
def test_native_search_readback_rejects_corrupt_or_incomplete_cost_records(change):
    row = deepcopy(case())
    if change == "ratio": row["ratios"]["wall_s"] = .91
    elif change == "order": row["repeats"][0]["order"] = list(ARMS[::-1])
    elif change == "missing": row["repeats"].pop()
    elif change == "bad_cost": row["repeats"][0]["arms"][ARMS[0]]["wall_s"] = float("nan")
    elif change == "bad_repeat": row["repeats"][0]["repeat"] = 1
    else: row["repeats"][0]["arms"].pop(ARMS[0])
    with pytest.raises(ValueError):
        timing_summary(row)
