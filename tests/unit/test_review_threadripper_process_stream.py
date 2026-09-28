from copy import deepcopy
import json
from io import StringIO

import pytest

from benchmark_tools.review_threadripper_process_stream import evaluate
from benchmark_tools.command_host_monitor import HostMonitor
from benchmark_tools import review_threadripper_process_stream as stream


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


def test_stream_accepts_native_exec_but_retains_foreign_migration_failure():
    policy, rows = fixture()
    for row in rows:
        child = deepcopy(row["snapshot"]["processes"][1])
        child.update(pid=30, name="python" if row["index"] == 0 else "FastTree")
        child["kernel_identity"].update(pid=30, tgid=30)
        row["snapshot"]["processes"].append(child)
    assert run(policy, rows)["sampled_process_policy_satisfied"]
    rows[1]["snapshot"]["processes"][-1]["cgroup"] = "/outside"
    assert not run(policy, rows)["sampled_process_policy_satisfied"]


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


def bound_fixture(tmp_path):
    def put(name, data):
        p = tmp_path / name
        p.write_text(json.dumps(data))
        return stream.record(p)
    policy, rows = fixture()
    for row in rows:
        row["snapshot"]["processes"][1]["cgroup"] = "/slurm/job_42/step_0/user/task_0"
    process_ref = put("process_policy.json", policy)
    support = put("support.json", {"synthetic": True})
    configuration = put("service_configuration.json", {"synthetic_service": "fixed"})
    environment = put("policy.json", dict(schema="threadripper_environment_policy_v1",
        decision="reviewed", host="bizon", plan_sha256="synthetic-plan",
        process_policy=process_ref, evidence=[support], configuration_files=[configuration],
        maximum_foreign_average_cores=.1,
        maximum_pressure_sample_period_s=3., maximum_pressure_percent=dict(cpu=10., memory=10., io=10.),
        maximum_sample_period_s=3.))
    preflight = put("preflight.json", dict(schema="threadripper_environment_preflight_v1",
        decision="passed", job_id=42, index=0, environment_policy=environment,
        plan_sha256="synthetic-plan", boot_id="test-boot", evidence=[support]))
    put("ready.json", dict(cgroup="0::/slurm/job_42/step_0/user/task_0\n"))
    put("done.json", dict(started_ns=11_000_000_000, finished_ns=13_000_000_000))
    (tmp_path / "host_processes.jsonl").write_text("".join(json.dumps(r) + "\n" for r in rows))
    for i, t in enumerate([10, 12, 14]):
        host = []
        for offset in [0, .2]:
            host.append(dict(started_monotonic_ns=int((t + offset) * 1e9),
                finished_monotonic_ns=int((t + offset + .1) * 1e9), errors=[],
                raw=dict(boot_id="test-boot", online_cpus="0-31",
                         cgroup_membership="0::/slurm/job_42/step_0\n"),
                optional={"host_" + r + "_pressure":
                    "some avg10=0 avg60=0 avg300=0 total=0\nfull avg10=0 avg60=0 avg300=0 total=0\n"
                    for r in ["cpu", "memory", "io"]}))
        put(f"point_{i:06d}.json", dict(host=host))
    return environment, preflight


def test_bound_audit_keeps_evidence_and_refuses_overwrite(tmp_path):
    refs = bound_fixture(tmp_path)
    ref, result = stream.audit(tmp_path, *refs, job_id=42, index=0)
    stream.check(ref)
    assert result["sampled_process_policy_satisfied"]
    assert result["sampled_environment_policy_satisfied"]
    assert result["job_id"] == 42 and not result["scientific_timings_admitted"]
    assert result["configuration_endpoint_hashes_verified"]
    assert stream.record(tmp_path / "service_configuration.json") in result["evidence"]
    for item in result["evidence"]:
        stream.check(item)
    with pytest.raises(FileExistsError):
        stream.audit(tmp_path, *refs, job_id=42, index=0)


@pytest.mark.parametrize("inventory", [None, [], "not-a-list"])
def test_bound_audit_requires_declared_configuration_inventory(tmp_path, inventory):
    policy_ref, preflight_ref = bound_fixture(tmp_path)
    policy_path = tmp_path / "policy.json"
    policy = json.loads(policy_path.read_text())
    if inventory is None:
        del policy["configuration_files"]
    else:
        policy["configuration_files"] = inventory
    policy_path.write_text(json.dumps(policy))
    policy_ref = stream.record(policy_path)
    preflight_path = tmp_path / "preflight.json"
    preflight = json.loads(preflight_path.read_text())
    preflight["environment_policy"] = policy_ref
    preflight_path.write_text(json.dumps(preflight))
    with pytest.raises(ValueError, match="configuration"):
        stream.audit(tmp_path, policy_ref, stream.record(preflight_path), job_id=42, index=0)
    assert not (tmp_path / "process_stream_review.json").exists()


@pytest.mark.parametrize("during_review", [False, True])
def test_configuration_change_rejected_without_incidental_preflight_reference(tmp_path, monkeypatch, during_review):
    refs = bound_fixture(tmp_path)
    configuration = tmp_path / "service_configuration.json"
    real = stream.evaluate
    def change(*args, **kwargs):
        result = real(*args, **kwargs)
        configuration.write_text("changed")
        return result
    if during_review:
        monkeypatch.setattr(stream, "evaluate", change)
    else:
        configuration.write_text("changed")
    with pytest.raises(ValueError, match="identity changed"):
        stream.audit(tmp_path, *refs, job_id=42, index=0)
    assert not (tmp_path / "process_stream_review.json").exists()


@pytest.mark.parametrize("change", ["policy", "preflight", "process", "support", "wrong_job", "wrong_index"])
def test_bound_audit_rejects_changed_or_wrong_attempt(tmp_path, change):
    refs = bound_fixture(tmp_path)
    kwargs = dict(job_id=42, index=0)
    if change == "wrong_job": kwargs["job_id"] = 43
    elif change == "wrong_index": kwargs["index"] = 1
    else:
        path = tmp_path / {"policy": "policy.json", "preflight": "preflight.json",
                           "process": "process_policy.json", "support": "support.json"}[change]
        path.write_text(path.read_text() + " ")
    with pytest.raises((ValueError, RuntimeError)):
        stream.audit(tmp_path, *refs, **kwargs)
    assert not (tmp_path / "process_stream_review.json").exists()


def test_bound_audit_retains_negative_verdict(tmp_path):
    refs = bound_fixture(tmp_path)
    path = tmp_path / "host_processes.jsonl"
    rows = path.read_text().splitlines()
    path.write_text(rows[0] + "\n" + rows[-1] + "\n")
    ref, result = stream.audit(tmp_path, *refs, job_id=42, index=0)
    assert not result["sampled_process_policy_satisfied"]
    assert (tmp_path / "process_stream_review.json").exists()


def test_changed_stream_during_review_rejected(tmp_path, monkeypatch):
    refs = bound_fixture(tmp_path)
    real = stream.evaluate
    def changed(*args, **kwargs):
        result = real(*args, **kwargs)
        with (tmp_path / "host_processes.jsonl").open("a") as f:
            f.write("{}\n")
        return result
    monkeypatch.setattr(stream, "evaluate", changed)
    with pytest.raises((ValueError, RuntimeError)):
        stream.audit(tmp_path, *refs, job_id=42, index=0)
    assert not (tmp_path / "process_stream_review.json").exists()


def test_pressure_failure_fails_combined_review_with_clean_processes(tmp_path):
    refs = bound_fixture(tmp_path)
    path = tmp_path / "point_000001.json"
    point = json.loads(path.read_text())
    point["host"][0]["errors"].append(dict(field="host_io_pressure", type="OSError"))
    path.write_text(json.dumps(point))
    _, result = stream.audit(tmp_path, *refs, job_id=42, index=0)
    assert result["sampled_process_policy_satisfied"]
    assert not result["sampled_environment_policy_satisfied"]
    assert not result["pressure_review"]["sampled_pressure_policy_satisfied"]
