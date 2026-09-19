import copy

import pytest

from benchmark_tools.run_lineage_lifecycle_control import assess


def args():
    return ["/user.slice/manager", "/system.slice/job_1/step_batch",
        {"cgroup": "0::/user.slice/manager/app.slice/test.service\n",
         "cpu_seconds": .751, "affinity": [3]}, "test.service", True,
        {"aggregate_deltas": [{"cpu_usec": 760000}]},
        {"root_minus_target_cpu_usec": 780000}, 3]


def test_expected_response():
    result = assess(*args())
    assert result["descendant_cpu_retention_response"]
    assert result["outside_target_response"]
    assert not result["scientific_timings_admitted"]


@pytest.mark.parametrize("change", ["membership", "unit", "cpu_low", "cpu_high", "affinity", "not_removed"])
def test_invalid_control(change):
    values = args()
    if change == "membership":
        values[2]["cgroup"] = "0::/another/test.service\n"
    elif change == "unit":
        values[3] = "other.service"
    elif change == "cpu_low":
        values[2]["cpu_seconds"] = .1
    elif change == "cpu_high":
        values[2]["cpu_seconds"] = 2
    elif change == "affinity":
        values[2]["affinity"] = [3, 4]
    else:
        values[4] = False
    with pytest.raises(ValueError):
        assess(*values)


def test_negative_or_insufficient_response_not_dropped():
    values = copy.deepcopy(args())
    values[5]["aggregate_deltas"][-1]["cpu_usec"] = 400000
    values[6]["root_minus_target_cpu_usec"] = -10000
    result = assess(*values)
    assert not result["descendant_cpu_retention_response"]
    assert not result["outside_target_response"]
    assert result["signed_outside_target_cpu_s"] == -.01
