from copy import deepcopy

import pytest

from benchmark_tools import describe_dual_cpu_aggregation as module
from tests.unit.test_measure_native_hierarchy_step import evidence
from tests.unit.test_probe_dual_cpu_brackets import paired


@pytest.mark.parametrize("n,width,expected", [(12, 5, [(0, 5), (5, 10), (10, 12)]),
    (10, 5, [(0, 5), (5, 10)]), (1, 60, [(0, 1)])])
def test_complete_endpoint_partition_including_remainder(n, width, expected):
    assert module.blocks(n, width) == expected
    assert [i for a, b in expected for i in range(a, b)] == list(range(n))


@pytest.mark.parametrize("n,width", [(0, 1), (1, 0), (True, 5), (3, 2.5)])
def test_invalid_partition(n, width):
    with pytest.raises(ValueError): module.blocks(n, width)


def fixture(evidence):
    points = paired(evidence)
    single = module.compare(*points, 21816)
    flags = [0] if not single["narrow"]["screen_passed"] else []
    screen = dict(narrow_intervals=[single["narrow"]], narrow_flagged_intervals=flags,
        original_screening=dict(original_threshold_screen=dict(intervals=[single["outer"]])),
        observation_window=module.compare(*points, 21816, enforce_gap=False))
    run = dict(index=0, method="test", job_id=21816, screening=screen, narrow_flagged_intervals=flags)
    return run, dict(job_id=21816, points=points, screening=deepcopy(screen))


def test_real_counter_fixture_replayed_without_erasing_originals(evidence):
    run, measured = fixture(evidence)
    original = deepcopy(measured)
    result = module.aggregate(run, measured)
    assert measured == original
    assert len(result["original_intervals"]) == 1
    for scale in result["scales"]:
        assert len(scale["blocks"]) == 1
        block = scale["blocks"][0]
        assert block["remainder"] is True and block["observation_steps"] == 1
        assert block["result"] == result["whole_observation_window"]
    assert result["scientific_timings_admitted"] is False


def test_coarser_passing_blocks_keep_burst_flags_and_final_partial_block(monkeypatch):
    calls = []
    def compare(a, b, job, enforce_gap=True):
        calls.append((a, b, enforce_gap))
        flagged = (a, b) == (3, 4)
        screen = dict(screen_passed=not flagged, wall_s=b-a)
        return dict(narrow=screen, outer=screen)
    monkeypatch.setattr(module, "compare", compare)
    points = list(range(13))
    singles = [compare(a, b, 1) for a, b in zip(points, points[1:])]
    screen = dict(narrow_intervals=[r["narrow"] for r in singles], narrow_flagged_intervals=[3],
        original_screening=dict(original_threshold_screen=dict(intervals=[r["outer"] for r in singles])),
        observation_window=compare(0, 12, 1, enforce_gap=False))
    run = dict(index=0, method="test", job_id=1, screening=screen, narrow_flagged_intervals=[3])
    result = module.aggregate(run, dict(job_id=1, points=points, screening=screen))
    scale = result["scales"][0]
    assert result["original_narrow_flagged_intervals"] == [3]
    assert scale["narrow_flagged_blocks"] == 0
    assert scale["flagged_original_intervals_in_narrow_passing_blocks"] == 1
    assert scale["blocks"][-1]["observation_steps"] == 2
    assert scale["blocks"][-1]["remainder"] is True
    assert (0, 5, False) in calls and (10, 12, False) in calls


@pytest.mark.parametrize("fault", ["stored", "counter", "flags", "whole", "job"])
def test_corrupt_evidence_rejected(evidence, fault):
    run, measured = fixture(evidence)
    if fault == "stored":
        measured["screening"]["narrow_intervals"][0]["wall_s"] += 1
    elif fault == "counter":
        measured["points"][1]["hierarchy_host_after"]["cpu_ticks"]["user"] += 1
    elif fault == "flags":
        run["narrow_flagged_intervals"] = [999]
    elif fault == "whole":
        for screen in (run["screening"], measured["screening"]):
            screen["observation_window"] = {}
    else:
        measured["job_id"] += 1
    with pytest.raises(ValueError):
        module.aggregate(run, measured)
