from copy import deepcopy

import pytest

from benchmark_tools.review_threadripper_process_policy import review


def fixture():
    service = dict(pid=10, created=1., cgroup="/system.slice/service", name="service")
    policy = dict(schema="threadripper_process_policy_v1", boot_id="boot-fixture",
                  review_reference="synthetic test only", ordinary_processes=[
                      dict(service, classification="ordinary_background", reason="synthetic service")])
    def snapshot(t):
        rows = [service, dict(pid=20, created=2., cgroup="/job/step", name="observer")]
        return dict(started_monotonic_s=t, finished_monotonic_s=t + .2, errors=[],
                    processes=[dict(row, user_s=0., system_s=0., observed_monotonic_s=t + .1)
                               for row in rows])
    return policy, snapshot(10.), snapshot(12.)


def run(policy, before, after, **kwargs):
    return review(policy, before, after, **dict(
        dict(boot_id="boot-fixture", job_scope="/job", observer_pid=20), **kwargs))


def test_matching_inventory_is_not_quiet_host_or_admission():
    data = fixture()
    saved = deepcopy(data)
    result = run(*data)
    assert result["process_policy_matched"]
    assert not result["controlled_workload_verified"]
    assert not result["scientific_timings_admitted"]
    assert data == saved


@pytest.mark.parametrize("cpu", [0., .01, 10.])
def test_unreviewed_process_blocks_even_when_idle(cpu):
    policy, before, after = fixture()
    policy["ordinary_processes"] = []
    after["processes"][0]["user_s"] = cpu
    result = run(policy, before, after)
    assert not result["process_policy_matched"]
    assert len(result["unreviewed_processes"]) == 2


@pytest.mark.parametrize("field,value", [("created", 3.), ("name", "python"),
                                        ("cgroup", "/job/step"), ("cgroup", "/other")])
def test_reuse_exec_name_and_migration_not_silently_accepted(field, value):
    policy, before, after = fixture()
    after["processes"][0][field] = value
    result = run(policy, before, after)
    assert not result["process_policy_matched"]
    assert result["changed_processes"]


def test_similar_cgroup_prefix_is_not_job_membership():
    policy, before, after = fixture()
    for sample in (before, after):
        sample["processes"][0]["cgroup"] = "/job_other"
    assert not run(policy, before, after)["process_policy_matched"]


@pytest.mark.parametrize("change", ["birth", "death", "permission", "decreased"])
def test_missing_or_changed_observation_blocks(change):
    policy, before, after = fixture()
    if change == "birth":
        after["processes"].append(dict(after["processes"][0], pid=30))
    elif change == "death":
        after["processes"].pop(0)
    elif change == "permission":
        after["errors"].append(dict(pid=30, type="AccessDenied"))
    else:
        before["processes"][0]["user_s"] = 1.
    assert not run(policy, before, after)["process_policy_matched"]


@pytest.mark.parametrize("change", ["wrong_boot", "root_scope", "missing_observer", "outside_observer",
    "reused_observer", "duplicate_pid", "nan_cpu", "bool_pid", "out_of_bounds", "overlap",
    "duplicate_policy", "missing_reason", "wrong_classification", "noncanonical", "empty",
    "wrong_policy_type", "double_slash"])
def test_invalid_evidence_rejected(change):
    policy, before, after = fixture()
    kwargs = {}
    if change == "wrong_boot": policy["boot_id"] = "other"
    elif change == "root_scope": kwargs["job_scope"] = "/"
    elif change == "missing_observer": after["processes"].pop()
    elif change == "outside_observer": after["processes"][1]["cgroup"] = "/other"
    elif change == "reused_observer": after["processes"][1]["created"] = 3.
    elif change == "duplicate_pid": after["processes"].append(dict(after["processes"][0]))
    elif change == "nan_cpu": after["processes"][0]["user_s"] = float("nan")
    elif change == "bool_pid": after["processes"][0]["pid"] = True
    elif change == "out_of_bounds": after["processes"][0]["observed_monotonic_s"] = 40.
    elif change == "overlap": before["finished_monotonic_s"] = 12.
    elif change == "duplicate_policy": policy["ordinary_processes"] *= 2
    elif change == "missing_reason": policy["ordinary_processes"][0]["reason"] = " "
    elif change == "wrong_classification": policy["ordinary_processes"][0]["classification"] = "scientific"
    elif change == "noncanonical": after["processes"][0]["cgroup"] = "/system.slice/../job"
    elif change == "wrong_policy_type": policy["ordinary_processes"] = {}
    elif change == "double_slash": kwargs["job_scope"] = "//job"
    else: after["processes"] = []
    with pytest.raises(ValueError):
        run(policy, before, after, **kwargs)


def test_cpu_diagnostic_does_not_override_policy_or_certify_eligibility():
    policy, before, after = fixture()
    after["processes"][0]["user_s"] = 10.
    result = run(policy, before, after)
    assert result["process_policy_matched"]
    assert result["cpu_diagnostic"]["status"] == "competing_cpu_observed"
    assert not result["controlled_workload_verified"]
    assert not result["scientific_timings_admitted"]
