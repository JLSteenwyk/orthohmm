import pytest

from benchmark_tools.summarize_fixture_memory_scopes import native_cpu, peak, summarize


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


def cpu_fixture():
    def sample(start, usage):
        return dict(started_monotonic_ns=start, finished_monotonic_ns=start + 1,
                    errors=[], raw=dict(boot_id="boot", online_cpus="0-31",
                    cgroup_membership="0::/job_1/step_0/user/task_0\n"),
                    optional={"cgroup_cpu.stat": f"usage_usec {usage}\nuser_usec {usage - 2}\nsystem_usec 2\n"})
    done = dict(snapshots=[sample(1, 10), sample(8, 30)], started_ns=3, finished_ns=7)
    completion = dict(status="anchor_only_at_boundaries", errors=[], command_finished_ns=7,
                      before=dict(scope="/sys/fs/cgroup/job_1/step_0/user"),
                      after=dict(scope="/sys/fs/cgroup/job_1/step_0/user"))
    return done, completion


def test_native_cpu_retains_scope_and_raw_differences():
    done, completion = cpu_fixture()
    result = native_cpu(done, 1, completion)
    assert result["delta_usec"] == dict(usage_usec=20, user_usec=20, system_usec=0)
    assert result["cpu_seconds"] == .000020
    assert result["final_whole_job_cpu_seconds"] is None
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("defect", ["time", "scope", "boot", "missing", "duplicate", "decrease", "error", "completion"])
def test_reject_invalid_cpu_evidence(defect):
    done, completion = cpu_fixture()
    after = done["snapshots"][1]
    if defect == "time":
        done["finished_ns"] = 9
    elif defect == "scope":
        completion["before"]["scope"] = "/sys/fs/cgroup/job_2/step_0/user"
    elif defect == "boot":
        after["raw"]["boot_id"] = "other"
    elif defect == "missing":
        after["optional"]["cgroup_cpu.stat"] = "usage_usec 30\n"
    elif defect == "duplicate":
        after["optional"]["cgroup_cpu.stat"] += "usage_usec 30\n"
    elif defect == "decrease":
        after["optional"]["cgroup_cpu.stat"] = "usage_usec 5\nuser_usec 3\nsystem_usec 2\n"
    elif defect == "error":
        after["errors"] = ["failed"]
    else:
        completion["command_finished_ns"] = 6
    with pytest.raises(ValueError):
        native_cpu(done, 1, completion)
