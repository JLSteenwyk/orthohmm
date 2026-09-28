from copy import deepcopy
import json
from io import StringIO

import pytest

from benchmark_tools.review_threadripper_process_stream import evaluate
from benchmark_tools.command_host_monitor import HostMonitor


def fixture():
    service = dict(pid=10, created=1., cgroup="/service", name="service")
    observer = dict(pid=20, created=2., cgroup="/job/step", name="observer")
    policy = dict(schema="threadripper_process_policy_v2", boot_id="test-boot",
        review_reference="synthetic only", ordinary_processes=[dict(service,
            classification="ordinary_background", reason="synthetic fixture")])
    rows = []
    for index, t in enumerate([10., 12., 14.]):
        processes = [dict(p, user_s=0., system_s=0., observed_monotonic_s=t + .1,
            kernel_identity=dict(pid=p["pid"], tgid=p["pid"], kthread=0,
                started_monotonic_s=t + .15, finished_monotonic_s=t + .2))
            for p in [service, observer]]
        rows.append(dict(index=index, observer_pid=20, interval={"forged": "ignored"},
            snapshot=dict(schema="threadripper_typed_process_snapshot_v1",
                boot_id="test-boot", started_monotonic_s=t, finished_monotonic_s=t + .3,
                errors=[], processes=processes)))
    return policy, rows


def run(policy, rows, **kwargs):
    options = dict(boot_id="test-boot", job_scope="/job", observer_pid=20,
                   launch=11., end=13., maximum_foreign_average_cores=.1,
                   maximum_sample_period_s=3.)
    options.update(kwargs)
    return evaluate((json.dumps(r) for r in rows), policy, **options)


def test_all_intervals_recomputed_not_trusted_and_not_admission():
    policy, rows = fixture()
    saved = deepcopy(rows)
    result = run(policy, rows)
    assert result["sampled_process_policy_satisfied"]
    assert result["intervals"] == result["policy_matched_intervals"] == 2
    assert not result["scientific_timings_admitted"]
    assert not result["controlled_workload_verified"]
    assert rows == saved


def test_real_monitor_serialization_composes_with_stream_reviewer():
    policy, rows = fixture()
    samples = iter(r["snapshot"] for r in rows)
    handle = StringIO()
    monitor = HostMonitor(handle, "/job", sample_fn=lambda: next(samples))
    monitor.observer_pid = 20
    for _ in rows:
        monitor.observe()
    assert monitor.summary(11., 13.)["command_bracketed_by_samples"]
    handle.seek(0)
    result = evaluate(handle, policy, boot_id="test-boot", job_scope="/job",
        observer_pid=20, launch=11., end=13., maximum_foreign_average_cores=.1,
        maximum_sample_period_s=3.)
    assert result["sampled_process_policy_satisfied"]


@pytest.mark.parametrize("change", ["idle_unknown", "cpu", "name", "boot", "type",
    "error", "index", "observer", "counter", "gap", "missing", "duplicate",
    "observation_error", "untyped", "malformed", "overlap"])
def test_middle_stream_corruption_never_disappears_in_endpoint_check(change):
    policy, rows = fixture()
    middle = rows[1]["snapshot"]
    process = middle["processes"][0]
    if change == "idle_unknown":
        new = deepcopy(process)
        new["pid"] = 30
        new["kernel_identity"].update(pid=30, tgid=30)
        middle["processes"].append(new)
    elif change == "cpu":
        process["user_s"] = 1.
        rows[2]["snapshot"]["processes"][0]["user_s"] = 1.
    elif change == "name": process["name"] = "other"
    elif change == "boot": middle["boot_id"] = "other"
    elif change == "type": del process["kernel_identity"]
    elif change == "error": middle["errors"].append({"type": "AccessDenied"})
    elif change == "index": rows[1]["index"] = 3
    elif change == "observer": rows[1]["observer_pid"] = 30
    elif change == "counter": rows[0]["snapshot"]["processes"][0]["user_s"] = 1.
    elif change == "gap":
        return_result = run(policy, rows, maximum_sample_period_s=1.)
        assert not return_result["sampled_process_policy_satisfied"]
        return
    elif change == "missing": rows.pop(1)
    elif change == "duplicate": rows.insert(1, deepcopy(rows[1]))
    elif change == "observation_error": rows[1] = dict(observation_error="OSError")
    elif change == "untyped": del middle["schema"]
    elif change == "malformed": rows[1] = []
    else: middle["started_monotonic_s"] = 10.
    assert not run(policy, rows)["sampled_process_policy_satisfied"]


@pytest.mark.parametrize("rows", [[], [dict(observation_error="OSError")]])
def test_absent_evidence_is_not_quiet(rows):
    assert not run(fixture()[0], rows)["sampled_process_policy_satisfied"]


@pytest.mark.parametrize("options", [dict(launch=10.), dict(end=15.)])
def test_native_bracketing_required(options):
    assert not run(*fixture(), **options)["sampled_process_policy_satisfied"]


@pytest.mark.parametrize("options", [dict(maximum_foreign_average_cores=True),
    dict(maximum_foreign_average_cores=float("nan")), dict(maximum_sample_period_s=0.),
    dict(launch=14., end=13.)])
def test_invalid_prospective_limits_rejected(options):
    with pytest.raises(ValueError):
        run(*fixture(), **options)


def test_constant_memory_single_pass_iterator():
    policy, rows = fixture()
    class Once:
        def __init__(self): self.used = False
        def __iter__(self):
            assert not self.used
            self.used = True
            for i in range(1000):
                row = deepcopy(rows[0])
                row["index"] = i
                sample = row["snapshot"]
                shift = 2. * i
                for key in ["started_monotonic_s", "finished_monotonic_s"]:
                    sample[key] += shift
                for p in sample["processes"]:
                    p["observed_monotonic_s"] += shift
                    for key in ["started_monotonic_s", "finished_monotonic_s"]:
                        p["kernel_identity"][key] += shift
                yield json.dumps(row)
    result = evaluate(Once(), policy, boot_id="test-boot", job_scope="/job",
        observer_pid=20, launch=11., end=2007., maximum_foreign_average_cores=.1,
        maximum_sample_period_s=3.)
    assert result["sampled_process_policy_satisfied"]
    assert result["intervals"] == 999
