import pytest

from benchmark_tools import describe_cpu_window_scales as module


@pytest.mark.parametrize("count,width,expected", [(2, 5, [(0, 1)]),
    (11, 5, [(0, 5), (5, 10)]), (12, 5, [(0, 5), (5, 10), (10, 11)])])
def test_disjoint_endpoint_grid_retains_tail(count, width, expected):
    assert module.windows(count, width) == expected
    covered = [i for start, end in expected for i in range(start, end)]
    assert covered == list(range(count-1))


@pytest.mark.parametrize("count,width", [(1, 5), (10, 0), (True, 1), (10, 1.5)])
def test_invalid_grid(count, width):
    with pytest.raises(ValueError):
        module.windows(count, width)


def test_uses_endpoint_counters_and_preserves_flags(monkeypatch):
    calls = []

    def compare(left, right, job, enforce_gap):
        calls.append((left, right))
        assert not enforce_gap
        return dict(narrow=dict(screen_passed=True, signed_unassigned_average_cores=.1))

    monkeypatch.setattr(module, "compare", compare)
    measured = dict(points=list(range(12)), job_id=1, screening=dict(narrow_flagged_intervals=[6],
        narrow_intervals=[dict(screen_passed=i != 6) for i in range(11)],
        observation_window=compare(0, 11, 1, False),
        original_screening=dict(original_threshold_screen=dict(whole_command_screen={"original": True}))))
    calls.clear()
    result = module.describe(measured)
    assert result["original_flag_indices"] == [6]
    assert calls == [(0, 5), (5, 10), (10, 11), (0, 10), (10, 11), (0, 11), (0, 11)]
    assert result["scales"][0]["windows"][1]["original_flag_indices"] == [6]
    assert result["scales"][0]["windows"][-1]["partial"] is True
    measured["screening"]["observation_window"] = {}
    with pytest.raises(ValueError, match="does not reproduce"):
        module.describe(measured)
