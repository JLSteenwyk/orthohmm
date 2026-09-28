import json
from copy import deepcopy
from types import SimpleNamespace

import pytest

from benchmark_tools import observe_threadripper_process_identity as capture
from benchmark_tools.review_threadripper_process_policy import review


@pytest.mark.parametrize("text", ["", "Pid: 10\nTgid: 10\n", "Pid: 11\nTgid: 10\nKthread: 1",
    "Pid: 10\nTgid: 11\nKthread: 1", "Pid: 10\nTgid: 10\nKthread: 2",
    "Pid: 10\nTgid: 10\nKthread: 1\nKthread: 1"])
def test_bad_status(text):
    with pytest.raises(ValueError):
        capture.kernel_fields(text, 10)


@pytest.mark.parametrize("flag", [0, 1])
def test_status(flag):
    result = capture.kernel_fields(f"Name: kworker\nPid: 10\nTgid: 10\nKthread: {flag}\n", 10)
    assert result == dict(pid=10, tgid=10, kthread=flag)


def fixture():
    policy = dict(schema="threadripper_process_policy_v2", boot_id="boot", review_reference="synthetic",
        ordinary_processes=[dict(pid=10, created=1., cgroup="/", name="worker-old",
                                 classification="reviewed_kernel_thread", reason="synthetic kernel")])
    def sample(t, name):
        rows = [dict(pid=10, created=1., cgroup="/", name=name),
                dict(pid=20, created=2., cgroup="/job/step", name="observer")]
        for row in rows:
            row.update(user_s=0., system_s=0., observed_monotonic_s=t + .1,
                kernel_identity=dict(pid=row["pid"], tgid=row["pid"], kthread=int(row["pid"] == 10),
                                     started_monotonic_s=t + .2, finished_monotonic_s=t + .3))
        return dict(schema="threadripper_typed_process_snapshot_v1", boot_id="boot", processes=rows,
                    errors=[], started_monotonic_s=t, finished_monotonic_s=t + .5)
    return policy, sample(10., "worker-old"), sample(12., "worker-new")


def run(*args):
    return review(*args, boot_id="boot", job_scope="/job", observer_pid=20)


def test_verified_reviewed_kernel_rename_keeps_cpu_visible():
    policy, before, after = fixture()
    after["processes"][0]["user_s"] = 3.
    result = run(policy, before, after)
    assert result["process_policy_matched"]
    assert result["reviewed_kernel_name_changes"] == [dict(pid=10, before="worker-old", after="worker-new")]
    assert result["cpu_diagnostic"]["sum_observed_foreign_average_cores"] == 1.5
    assert not result["controlled_workload_verified"] and not result["scientific_timings_admitted"]


@pytest.mark.parametrize("change", ["exec", "replacement", "step", "birth", "exit"])
def test_native_job_churn_is_not_foreign_work(change):
    policy, before, after = fixture()
    for sample in [before, after]:
        row = deepcopy(sample["processes"][1])
        row.update(pid=30, name="native-child", user_s=5.)
        row["kernel_identity"].update(pid=30, tgid=30)
        sample["processes"].append(row)
    row = after["processes"][-1]
    if change == "exec": row["name"] = "FastTree"
    elif change == "replacement": row.update(created=9., name="mcl", user_s=0.)
    elif change == "step": row["cgroup"] = "/job/other-step"
    elif change == "birth": before["processes"].pop()
    else: after["processes"].pop()
    result = run(policy, before, after)
    assert result["process_policy_matched"]
    assert result["cpu_diagnostic"]["sum_observed_foreign_average_cores"] == 0.
    if change in {"exec", "replacement", "step"}:
        assert result["observed_in_job_identity_changes"][0]["pid"] == 30
    assert not result["scientific_timings_admitted"]


@pytest.mark.parametrize("change", ["leaves_job", "counter", "error", "observer", "backwards_creation"])
def test_native_scope_rule_does_not_waive_uncertain_or_external_changes(change):
    policy, before, after = fixture()
    for sample in [before, after]:
        row = deepcopy(sample["processes"][1])
        row.update(pid=30, name="native", user_s=5.)
        row["kernel_identity"].update(pid=30, tgid=30)
        sample["processes"].append(row)
    if change == "leaves_job": after["processes"][-1]["cgroup"] = "/other"
    elif change == "counter": after["processes"][-1]["user_s"] = 0.
    elif change == "error": after["errors"].append(dict(pid=30, type="NoSuchProcess"))
    elif change == "backwards_creation": after["processes"][-1]["created"] = 1.
    else:
        after["processes"][1]["name"] = "different-observer"
        with pytest.raises(ValueError, match="Observer identity changed"):
            run(policy, before, after)
        return
    assert not run(policy, before, after)["process_policy_matched"]


@pytest.mark.parametrize("change", ["spoofed_name", "reused_pid", "migration", "counter", "unapproved",
                                   "ordinary_kernel", "death"])
def test_kernel_rule_does_not_hide_other_changes(change):
    policy, before, after = fixture()
    if change == "spoofed_name": after["processes"][0]["kernel_identity"]["kthread"] = 0
    elif change == "reused_pid": after["processes"][0]["created"] = 9.
    elif change == "migration": after["processes"][0]["cgroup"] = "/job/step"
    elif change == "counter": before["processes"][0]["user_s"] = 2.
    elif change == "unapproved": policy["ordinary_processes"] = []
    elif change == "ordinary_kernel": policy["ordinary_processes"][0]["classification"] = "ordinary_background"
    else: after["processes"].pop(0)
    assert not run(policy, before, after)["process_policy_matched"]


@pytest.mark.parametrize("change", ["missing", "wrong_boot", "wrong_pid", "bool_type", "time", "error"])
def test_missing_or_invalid_type_evidence_rejected(change):
    data = fixture()
    row = data[2]["processes"][0]
    if change == "missing": del row["kernel_identity"]
    elif change == "wrong_boot": data[2]["boot_id"] = "other"
    elif change == "wrong_pid": row["kernel_identity"]["pid"] = 30
    elif change == "bool_type": row["kernel_identity"]["kthread"] = True
    elif change == "time": row["kernel_identity"]["finished_monotonic_s"] = 99.
    else: row["kernel_identity_error"] = "AccessDenied"
    with pytest.raises(ValueError):
        run(*data)


@pytest.mark.parametrize("change", [None, "pid_reuse", "group", "type", "name", "disappeared"])
def test_live_detail_binding(monkeypatch, change):
    calls = []
    def process(pid):
        calls.append(pid)
        final = len(calls) > 1
        return SimpleNamespace(create_time=lambda: 9. if final and change == "pid_reuse" else 1.,
            is_running=lambda: change != "disappeared", name=lambda: "other" if change == "name" else "user")
    monkeypatch.setattr(capture.psutil, "Process", process)
    groups = iter(["/other", "/changed" if change == "group" else "/other"])
    monkeypatch.setattr(capture, "membership", lambda pid: next(groups))
    flags = iter([0, 1 if change == "type" else 0])
    monkeypatch.setattr(capture.Path, "read_text", lambda self: f"Pid: 10\nTgid: 10\nKthread: {next(flags)}")
    row = dict(pid=10, created=1., cgroup="/other", name="user")
    if change:
        with pytest.raises(ValueError): capture.details(row)
    else:
        assert capture.details(row)["kthread"] == 0
        assert calls == [10, 10]


def test_capture_failure_preserved_and_no_overwrite(tmp_path, monkeypatch):
    monkeypatch.setattr(capture, "enriched_snapshot", lambda: (_ for _ in ()).throw(RuntimeError()))
    output = tmp_path / "capture.json"
    with pytest.raises(RuntimeError): capture.capture(output)
    assert json.loads(output.read_text())["status"] == "capture_failed"
    with pytest.raises(FileExistsError): capture.capture(output)
