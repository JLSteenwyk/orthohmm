import json

import pytest

from benchmark_tools import replay_full_node_control as module
from tests.unit.test_validate_full_node_control import fixture


def save(path, value):
    path.write_text(json.dumps(value))


def archive(tmp_path, monkeypatch, mode="steady"):
    args = fixture(mode)
    work = tmp_path / "workload"
    work.mkdir()
    save(work / "workload_go.json", {"go": True})
    save(work / "workload_ready.json", args["ready"])
    save(work / "workload_done.json", args["done"])
    for before, after in zip(args["ready"]["workers"], args["done"]["workers"]):
        save(work / f'ready_{after["cpu"]}.json', before)
        save(work / f'done_{after["cpu"]}.json', after)
    controller = dict(status="completed")
    if mode == "contended":
        save(work / "competitor_ready.json", args["competitor_ready"])
        save(work / "competitor_done.json", args["competitor"])
        controller["competitor_exit_code"] = 0
    save(tmp_path / "controller.json", controller)
    points = [dict(native_membership="native", host=[
        dict(started_monotonic_ns=t-1, raw=dict(cgroup_membership="batch")),
        dict(finished_monotonic_ns=t+1)]) for t in (0, 5_000_000_000, 15_000_000_000, 25_000_000_000)]
    screens = [dict(screen_passed=True, reasons=[]) for _ in range(3)]
    report = dict(points=points, native=args["native"], screening=dict(narrow_intervals=screens))

    def measurement(directory, job, command, expected_timeout_s):
        assert job == 17 and expected_timeout_s == 60
        assert command == ["/pinned/python", "-B", "/recipe/benchmark_tools/full_node_control_workload.py",
                           "--directory", "/remote/trial/workload", "--worker",
                           "churn" if mode == "churn" else "steady"]
        return dict(measured=report)

    monkeypatch.setattr(module, "replay_measurement", measurement)
    validation = module.validate_witnesses(**args)
    expected = dict(mode=mode, job_id=17, status="workload_validated", validation=validation,
        common_intervals=[1], common_narrow_flagged=[],
        positive_control_detected=False if mode == "contended" else None,
        scientific_timings_admitted=False)
    save(tmp_path / "trial.json", expected)
    return expected


def replay(path, mode="steady"):
    return module.replay(path, "/remote/trial", "/recipe", "/pinned/python", mode, 17)


@pytest.mark.parametrize("mode", ["steady", "churn", "contended"])
def test_replays_workload_and_common_window(tmp_path, monkeypatch, mode):
    expected = archive(tmp_path, monkeypatch, mode)
    result = replay(tmp_path, mode)
    assert result["trial"] == expected
    assert not result["scientific_timings_admitted"]
    if mode == "contended":
        assert result["trial"]["positive_control_detected"] is False


@pytest.mark.parametrize("fault", ["raw_worker", "go", "failed_worker", "extra_worker", "summary", "competitor", "symlink"])
def test_archive_inconsistency_rejected(tmp_path, monkeypatch, fault):
    archive(tmp_path, monkeypatch)
    work = tmp_path / "workload"
    if fault == "raw_worker":
        save(work / "done_0.json", {})
    elif fault == "go":
        save(work / "workload_go.json", {"go": 1})
    elif fault == "failed_worker":
        save(work / "failed_0.json", {})
    elif fault == "extra_worker":
        save(work / "ready_99.json", {})
    elif fault == "summary":
        save(tmp_path / "trial.json", {})
    elif fault == "competitor":
        save(work / "competitor_done.json", {})
    else:
        (tmp_path / "link").symlink_to(work / "done_0.json")
    with pytest.raises(ValueError):
        replay(tmp_path)


def test_evidence_mutation_during_replay_rejected(tmp_path, monkeypatch):
    archive(tmp_path, monkeypatch)
    original = module.validate_witnesses

    def change(*args, **kwargs):
        result = original(*args, **kwargs)
        save(tmp_path / "controller.json", {"status": "changed"})
        return result

    monkeypatch.setattr(module, "validate_witnesses", change)
    with pytest.raises(ValueError):
        replay(tmp_path)
