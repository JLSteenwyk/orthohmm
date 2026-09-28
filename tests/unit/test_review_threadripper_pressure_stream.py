from copy import deepcopy

import pytest

from benchmark_tools.review_threadripper_pressure_stream import evaluate


def fixture():
    points = []
    for t in [10, 12, 14]:
        samples = []
        for offset in [0, .2]:
            samples.append(dict(started_monotonic_ns=int((t + offset) * 1e9),
                finished_monotonic_ns=int((t + offset + .1) * 1e9), errors=[],
                raw=dict(boot_id="boot\n", online_cpus="0-31\n",
                         cgroup_membership="0::/slurm/job_42/step_0\n"),
                optional={"host_" + r + "_pressure":
                    "some avg10=0.00 avg60=0.00 avg300=0.00 total=0\n"
                    "full avg10=0.00 avg60=0.00 avg300=0.00 total=0\n"
                    for r in ["cpu", "memory", "io"]}))
        points.append(dict(host=samples))
    return points


def run(points, **kwargs):
    options = dict(boot_id="boot", job_scope="/slurm/job_42", launch_ns=11_000_000_000,
                   end_ns=13_000_000_000, limits=dict(cpu=10., memory=10., io=10.),
                   maximum_period_s=3.)
    options.update(kwargs)
    return evaluate(iter(points), **options)


def test_complete_stream_is_not_timing_admission():
    points = fixture()
    saved = deepcopy(points)
    result = run(points)
    assert result["sampled_pressure_policy_satisfied"]
    assert result["intervals"] == 2
    assert not result["scientific_timings_admitted"]
    assert not result["controlled_workload_verified"]
    assert saved == points


@pytest.mark.parametrize("resource", ["cpu", "memory", "io"])
def test_burst_not_hidden_by_averaging(resource):
    points = fixture()
    for p in points[1:]:
        for s in p["host"]:
            key = "host_" + resource + "_pressure"
            s["optional"][key] = s["optional"][key].replace("total=0", "total=300000", 1)
    result = run(points)
    assert not result["sampled_pressure_policy_satisfied"]
    assert result["maximum_observed_some_percent"][resource] == pytest.approx(15.)
    assert result["failures"][resource + "_pressure_bound_exceeded"] == 1


@pytest.mark.parametrize("change", ["read_error", "boot", "group", "cpu_set", "missing",
    "counter", "overlap", "brackets", "gap", "late", "early", "malformed"])
def test_incomplete_or_changed_evidence_rejected(change):
    points = fixture()
    s = points[1]["host"][0]
    if change == "read_error": s["errors"].append(dict(field="unrelated", type="OSError"))
    elif change == "boot": s["raw"]["boot_id"] = "other"
    elif change == "group": s["raw"]["cgroup_membership"] = "0::/other\n"
    elif change == "cpu_set": s["raw"]["online_cpus"] = "0\n"
    elif change == "missing": del s["optional"]["host_io_pressure"]
    elif change == "counter":
        points[0]["host"][0]["optional"]["host_io_pressure"] = s["optional"]["host_io_pressure"].replace("total=0", "total=1")
    elif change == "overlap": s["started_monotonic_ns"] = 10_000_000_000
    elif change == "brackets": points[1]["host"].pop()
    elif change == "gap": points.pop(1)
    elif change == "late": return assert_failed(run(points, launch_ns=10_000_000_000))
    elif change == "early": return assert_failed(run(points, end_ns=15_000_000_000))
    else: s["optional"]["host_cpu_pressure"] = "garbage"
    assert_failed(run(points))


def assert_failed(result):
    assert not result["sampled_pressure_policy_satisfied"]


@pytest.mark.parametrize("points", [[], [{}]])
def test_absent_evidence_never_passes(points):
    assert_failed(run(points))


@pytest.mark.parametrize("options", [dict(maximum_period_s=0), dict(limits=dict(cpu=1)),
    dict(limits=dict(cpu=float("nan"), io=1, memory=1)), dict(launch_ns=True),
    dict(end_ns=1), dict(limits=dict(cpu=101, io=1, memory=1))])
def test_invalid_bounds_rejected(options):
    with pytest.raises(ValueError):
        run(fixture(), **options)
