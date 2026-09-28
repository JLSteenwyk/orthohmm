from copy import deepcopy

import pytest

from benchmark_tools.probe_threadripper_report_size import expand, repeat


def templates():
    row = {"value": [1]}
    a = dict(schema="threadripper_scaling_v3", point_records=[row],
             screening=dict(intervals=[row], original_screening=dict(
                 hierarchy_intervals=[row], original_threshold_screen=dict(intervals=[row]))))
    b = dict(context=dict(intervals=[row], observations=2))
    return a, b


def test_expands_every_duration_dependent_list_without_aliasing():
    a, b = templates()
    original = deepcopy((a, b))
    result, context = expand(a, b, 5)
    assert (a, b) == original
    assert len(result["point_records"]) == 5
    s = result["screening"]
    o = s["original_screening"]
    for rows in (s["intervals"], o["hierarchy_intervals"], o["original_threshold_screen"]["intervals"], context["context"]["intervals"]):
        assert len(rows) == 4
        rows[0]["value"].append(2)
        assert rows[1]["value"] == [1]
    assert s["narrow_flagged_intervals"] == list(range(4))
    assert o["original_threshold_screen"]["flagged_intervals"] == list(range(4))
    for report in (result, context):
        assert report["synthetic"]
        assert report["status"] == "synthetic_report_size_probe"
        assert not report["scientific_timings_admitted"]


@pytest.mark.parametrize("count", [0, 1, -1, True, 2.0, 85803])
def test_invalid_count(count):
    with pytest.raises(ValueError):
        expand(*templates(), count)


def test_empty_template_rejected():
    with pytest.raises(ValueError):
        repeat([], 5)


def test_historical_schema_rejected():
    a, b = templates()
    a["schema"] = "other"
    with pytest.raises(ValueError):
        expand(a, b, 5)
