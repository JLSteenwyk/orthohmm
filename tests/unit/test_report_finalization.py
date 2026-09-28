import json

import pytest

from benchmark_tools.report_finalization import observe, validate


def fixture():
    after = dict(scope="/job", finished_ns=9)
    receipt = dict(schema="threadripper_report_finalization_v1", status="reporting_completed",
        job_id=42, scope="/job", started_ns=10, finished_ns=20, reporting_wall_s=1e-8,
        job_memory=dict(started_ns=21), scientific_timings_admitted=False)
    return receipt, after


def test_validate_phase():
    receipt, after = fixture()
    assert validate(receipt, 42, after) == receipt["job_memory"]


@pytest.mark.parametrize("key,value", [
    ("schema", "unknown"), ("status", "reporting_failed"), ("job_id", True),
    ("job_id", 43), ("scope", "/other"), ("scientific_timings_admitted", True),
    ("started_ns", 8), ("finished_ns", 21), ("started_ns", True),
    ("reporting_wall_s", True), ("reporting_wall_s", 1), ("reporting_wall_s", float("nan")),
])
def test_reject_invalid_phase(key, value):
    receipt, after = fixture()
    receipt[key] = value
    with pytest.raises(ValueError):
        validate(receipt, 42, after)


@pytest.mark.parametrize("body_failure,memory_failure", [(False, False), (True, False), (False, True), (True, True)])
def test_observation_and_failure_preservation(tmp_path, body_failure, memory_failure):
    times = iter([10, 20])
    events = []
    def read(scope):
        events.append("memory")
        assert scope == "/job"
        if memory_failure:
            raise OSError("memory unavailable")
        return dict(started_ns=21)
    def run():
        with observe(tmp_path, 42, "/job", read, clock=lambda: next(times)):
            events.append("body")
            if body_failure:
                raise RuntimeError("serialization failed")
    if body_failure:
        with pytest.raises(RuntimeError, match="serialization"):
            run()
    elif memory_failure:
        with pytest.raises(OSError, match="memory"):
            run()
    else:
        run()
    receipt = json.loads((tmp_path / "report_finalization.json").read_text())
    assert receipt["reporting_wall_s"] == 1e-8
    assert events == ["body", "memory"]
    assert receipt["status"] == ("reporting_memory_observation_failed" if memory_failure else
                                  "reporting_failed" if body_failure else "reporting_completed")
    if body_failure:
        assert receipt["error_type"] == "RuntimeError"
    if memory_failure:
        assert receipt["memory_error_type"] == "OSError"
