from copy import deepcopy

import pytest

from benchmark_tools.replay_threadripper_scaling import replay_affinity


def point():
    return dict(native_membership="0::/slurm/job_2/step_0/user/task_0\n",
                root_context={"host_after": {"finished_ns": 10}}, thread_affinity=dict(
                    scope="/sys/fs/cgroup/slurm/job_2/step_0/user", allowed_cpus=list(range(32)),
                    full_run_affinity_verified=False, started_ns=11, finished_ns=12,
                    initial_tids=[4], final_tids=[4], threads=[dict(tid=4, start_ticks=1,
                    cgroup="0::/slurm/job_2/step_0/user/task_0\n", affinity=[0], outside_cpus=[])],
                    errors=[], violating_tids=[], status="observed_within_affinity"))


def test_replay_all_outcomes():
    p = point()
    assert replay_affinity(p, 2) == "observed_within_affinity"
    a = p["thread_affinity"]
    a["threads"][0].update(affinity=[96], outside_cpus=[96])
    a.update(violating_tids=[4], status="violation")
    assert replay_affinity(p, 2) == "violation"
    a.update(threads=[], violating_tids=[], status="incomplete", errors=[dict(error="exited")])
    assert replay_affinity(p, 2) == "incomplete"


@pytest.mark.parametrize("change", ["scope", "policy", "admission", "window", "gap", "outside", "status", "membership", "duplicate"])
def test_rejects_tampering(change):
    p = deepcopy(point())
    a = p["thread_affinity"]
    if change == "scope":
        a["scope"] += "1"
    elif change == "policy":
        a["allowed_cpus"] = list(range(64))
    elif change == "admission":
        a["full_run_affinity_verified"] = True
    elif change == "window":
        a["started_ns"] = 9
    elif change == "gap":
        a["threads"] = []
    elif change == "outside":
        a["threads"][0]["affinity"] = [96]
    elif change == "status":
        a["status"] = "violation"
    elif change == "membership":
        a["threads"][0]["cgroup"] = "0::/slurm/job_2/step_01/user/task_0\n"
    else:
        a["threads"].append(deepcopy(a["threads"][0]))
    with pytest.raises(ValueError):
        replay_affinity(p, 2)
