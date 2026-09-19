import copy
import json
from pathlib import Path

import pytest

from benchmark_tools.run_lineage_lifecycle_control import assess
from benchmark_tools.probe_cgroup_lineage import compare


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


def test_retained_dgx_raw_snapshots_reproduce_all_three_controls():
    root = Path(__file__).resolve().parents[2] / "benchmark_tools/results/lineage_lifecycle_21989"
    identity = json.loads((root / "identity.json").read_text())
    report = json.loads((root / "report.json").read_text())
    assert identity["job"] == 21989 and identity["affinity"] == [0, 1]
    assert [r["index"] for r in report["trials"]] == [0, 1, 2]
    assert report["all_controls_met"] and report["sources_unchanged"]
    for index, row in enumerate(report["trials"]):
        directory = root / f"trial_{index}"
        before = json.loads((directory / "before.json").read_text())
        after = json.loads((directory / "after.json").read_text())
        removal = json.loads((directory / "removal.json").read_text())
        service = json.loads((directory / "service.json").read_text())
        assert service["returncode"] == 0
        load = json.loads(service["stdout"])
        assert load == row["load"]
        assert not removal["observations"][-1]["exists"]
        assert removal["scope"] == load["cgroup"].strip()[3:]
        removed_at = removal["observations"][-1]["monotonic_ns"]
        assert before["target"]["rows"][-1]["finished_ns"] < removed_at
        assert removed_at < after["manager"]["rows"][0]["started_ns"]
        assert compare(before["manager"], after["manager"]) == row["manager"]
        assert compare(before["target"], after["target"]) == row["outside"]
        assert assess(identity["manager"], identity["target"], load,
            f"orthohmm-lineage-21989-{index}.service", True, row["manager"], row["outside"], 0) == row["result"]
        # Independently recompute the two response scalars from retained integer counters.
        def delta(kind, position):
            def counter(point):
                raw = point[kind]["rows"][position]["raw"]
                return int(dict(line.split() for line in raw.splitlines())["usage_usec"])
            return counter(after) - counter(before)
        assert delta("manager", -1) / 1e6 == row["result"]["manager_cpu_s"]
        outside = (delta("target", 0) - delta("target", -1)) / 1e6
        assert outside == row["result"]["signed_outside_target_cpu_s"]
        assert json.loads((directory / "result.json").read_text()) == row
