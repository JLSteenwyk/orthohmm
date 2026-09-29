import pytest

from benchmark_tools.summarize_fixture_memory_scopes import peak, summarize


def observation(value, scope="/job"):
    return dict(raw={"memory.peak": str(value)}, errors=[], scope=scope, started_ns=1, finished_ns=2)


def replay():
    return dict(status="threadripper_scaling_measurement_replayed",
                native_completion=dict(status="anchor_only_at_boundaries"),
                memory=observation(10, "/job/step"),
                job_memory=dict(before=observation(5), after=observation(15)),
                report_finalization=dict(job_memory=observation(20)))


def test_scopes_remain_separate_and_final_unknown():
    row = summarize(replay())
    assert row["final_whole_job_peak_bytes"] is None
    assert row["scientific_timings_admitted"] is False
    assert row["measurements"]["native_step"]["bytes"] == 10
    assert row["measurements"]["job_through_reporting"]["bytes"] == 20


@pytest.mark.parametrize("value", ["-1", "NaN", "max", "1.5", ""])
def test_bad_peak(value):
    with pytest.raises(ValueError):
        peak(observation(value))


def test_incomplete_observation():
    row = observation(10)
    row["errors"] = ["failed read"]
    with pytest.raises(ValueError):
        peak(row)


def test_wrong_scope():
    row = replay()
    row["memory"]["scope"] = "/other/step"
    with pytest.raises(ValueError, match="scopes"):
        summarize(row)


def test_decreasing_peak():
    row = replay()
    row["report_finalization"]["job_memory"]["raw"]["memory.peak"] = "14"
    with pytest.raises(ValueError, match="peaks"):
        summarize(row)
