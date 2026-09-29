import copy

import pytest

from benchmark_tools import derive_threadripper_resources as module


def fixture():
    def sample(start, usage):
        return dict(started_monotonic_ns=start, finished_monotonic_ns=start + 1,
            errors=[], raw=dict(boot_id="boot", online_cpus="0-31",
            cgroup_membership="0::/job_1/step_0/user/task_0\n"),
            optional={"cgroup_cpu.stat": f"usage_usec {usage}\nuser_usec {usage - 2}\nsystem_usec 2\n"})
    done = dict(snapshots=[sample(1, 10), sample(8, 30)], started_ns=3, finished_ns=7)
    completion = dict(status="anchor_only_at_boundaries", errors=[], command_finished_ns=7,
        before=dict(scope="/sys/fs/cgroup/job_1/step_0/user"),
        after=dict(scope="/sys/fs/cgroup/job_1/step_0/user"))
    def memory(value, scope):
        return dict(raw={"memory.peak": str(value)}, errors=[], scope=scope, started_ns=10, finished_ns=11)
    report = dict(status="threadripper_scaling_measurement_replayed", native_wall_s=4e-9,
        measured=dict(job_id=1, native=copy.deepcopy(done)), native_completion=completion,
        memory=memory(10, "/job_1/step_0"),
        job_memory=dict(before=memory(5, "/job_1"), after=memory(15, "/job_1")),
        report_finalization=dict(job_memory=memory(20, "/job_1")),
        native_outcome="exited_zero", native_exit_code=0)
    return report, done


def test_primary_values_and_scopes_keep_job_peak_separate():
    report, done = fixture()
    result = module.endpoints(report, done, 1)
    assert result["primary"] == dict(wall_seconds=4e-9, cpu_seconds=.000020, peak_memory_bytes=10)
    assert result["memory"]["measurements"]["job_through_reporting"]["bytes"] == 20
    assert result["memory"]["final_whole_job_peak_bytes"] is None
    assert result["primary_scopes"] == module.SCOPES
    assert result["controlled_timing_admitted"] is False


@pytest.mark.parametrize("defect", ["job", "done", "wall", "peak", "cpu"])
def test_bad_resource_bindings_rejected(defect):
    report, done = fixture()
    if defect == "job":
        report["measured"]["job_id"] = 2
    elif defect == "done":
        report["measured"]["native"]["finished_ns"] = 6
    elif defect == "wall":
        report["native_wall_s"] = 1
    elif defect == "peak":
        report["memory"]["raw"]["memory.peak"] = "-1"
    else:
        done["snapshots"][1]["optional"]["cgroup_cpu.stat"] = "usage_usec 1\nuser_usec 0\nsystem_usec 1\n"
        report["measured"]["native"] = copy.deepcopy(done)
    with pytest.raises(ValueError):
        module.endpoints(report, done, 1)


def test_failed_native_outcome_is_preserved():
    report, done = fixture()
    report.update(native_outcome="exited_nonzero", native_exit_code=1)
    result = module.endpoints(report, done, 1)
    assert result["native_exit_code"] == 1
    assert result["native_outcome"] == "exited_nonzero"
    assert result["controlled_timing_admitted"] is False


def test_existing_output_refused_before_reading_inputs(tmp_path):
    with pytest.raises(FileExistsError):
        module.derive(tmp_path, 1, tmp_path / "absent", "unused", tmp_path)


def test_protocol_digest_checked_before_parsing(tmp_path):
    path = tmp_path / "protocol.json"
    path.write_text("invalid JSON")
    with pytest.raises(ValueError, match="digest"):
        module.protocol(path, "wrong")
