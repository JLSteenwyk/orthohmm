import copy

import pytest

from benchmark_tools.validate_full_node_control import common_intervals, validate_witnesses


def fixture(mode="steady"):
    cpus = list(range(20))
    workers, ready_workers = [], []
    for cpu in cpus:
        identity = dict(pid=100+cpu, parent_pid=90, cpu=cpu, allowed=cpus,
                        affinity=[cpu], membership="native")
        ready_workers.append(identity)
        workers.append(dict(identity, final_affinity=[cpu], final_membership="native",
            mode="churn" if mode == "churn" else "steady", duration_s=20, creation_cap=200000,
            started_ns=2_000_000_000, finished_ns=22_000_000_000, self_cpu_s=19,
            waited_child_user_s=0, waited_child_system_s=0, creations=10 if mode == "churn" else 0,
            creation_cap_reached=False))
    ready = dict(pid=90, cpus=cpus, mode=workers[0]["mode"], membership="native", workers=ready_workers)
    done = dict(ready, workers=workers, duration_s=20, creation_cap=200000,
                statuses={str(row["pid"]): 0 for row in workers},
                started_ns=1_500_000_000, finished_ns=23_000_000_000)
    args = dict(mode=mode, ready=ready, done=done, native=dict(exit_code=0, timed_out=False,
                started_ns=1_000_000_000, finished_ns=24_000_000_000),
                native_membership="native", observer_membership="batch")
    if mode == "contended":
        identity = dict(pid=200, parent_pid=80, affinity=[0], membership="batch")
        args.update(competitor_ready=identity, competitor=dict(identity, final_affinity=[0],
            final_membership="batch", started_ns=2_100_000_000, finished_ns=22_100_000_000,
            self_cpu_s=9, waited_child_user_s=0, waited_child_system_s=0, creations=0,
            creation_cap_reached=False))
    return copy.deepcopy(args)


@pytest.mark.parametrize("mode", ["steady", "churn", "contended"])
def test_valid_witnesses(mode):
    result = validate_witnesses(**fixture(mode))
    assert result["status"] == "workload_witnesses_validated"
    assert result["scientific_timings_admitted"] is False


def test_explicit_competitor_scope_does_not_weaken_historical_default():
    args = fixture("contended")
    scope = "0::/user.slice/user-1000.slice/user@1000.service/app.slice/owned.service\n"
    args["competitor_ready"]["membership"] = scope
    args["competitor"].update(membership=scope, final_membership=scope)
    with pytest.raises(ValueError):
        validate_witnesses(**args)
    assert validate_witnesses(**args, expected_competitor_membership=scope)["status"] == "workload_witnesses_validated"
    args["competitor"]["final_membership"] = "batch"
    with pytest.raises(ValueError):
        validate_witnesses(**args, expected_competitor_membership=scope)


@pytest.mark.parametrize("field,value", [("self_cpu_s", float("nan")), ("self_cpu_s", 0),
    ("parent_pid", 999), ("membership", "batch"), ("final_affinity", [1]),
    ("finished_ns", 21_000_000_000), ("creation_cap_reached", True),
    ("creations", True), ("mode", "churn"), ("duration_s", 10), ("allowed", [0])])
def test_invalid_worker_rejected(field, value):
    args = fixture()
    args["done"]["workers"][0][field] = value
    with pytest.raises(ValueError):
        validate_witnesses(**args)


@pytest.mark.parametrize("field,value", [("self_cpu_s", 4.9), ("self_cpu_s", 21.1),
    ("membership", "native"), ("affinity", [1]), ("started_ns", 3_000_000_000)])
def test_invalid_competitor_rejected(field, value):
    args = fixture("contended")
    args["competitor"][field] = value
    with pytest.raises(ValueError):
        validate_witnesses(**args)


def test_missing_and_failed_status_rejected():
    for statuses in ({}, {str(100+i): 1 if i == 0 else 0 for i in range(20)}):
        args = fixture()
        args["done"]["statuses"] = statuses
        with pytest.raises(ValueError, match="exit status"):
            validate_witnesses(**args)


def test_common_interval_requires_entire_host_window():
    points = [dict(host=[dict(started_monotonic_ns=t-2), dict(finished_monotonic_ns=t+2)])
              for t in (0, 10, 20, 30)]
    assert common_intervals(points, dict(common_started_ns=0, common_finished_ns=30)) == [1]


@pytest.mark.parametrize("mode", ["steady", "contended"])
def test_empty_readiness_is_not_evidence(mode):
    args = fixture(mode)
    if mode == "steady":
        args["ready"]["workers"][0] = {}
    else:
        args["competitor_ready"] = {}
    with pytest.raises(ValueError, match="readiness"):
        validate_witnesses(**args)
