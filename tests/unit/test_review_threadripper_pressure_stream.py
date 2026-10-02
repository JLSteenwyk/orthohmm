from copy import deepcopy

import pytest

from benchmark_tools.review_threadripper_pressure_stream import evaluate, native_pressure_role


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
@pytest.mark.parametrize("role", ["eligibility", "diagnostic_only"])
def test_burst_not_hidden_by_averaging(resource, role):
    points = fixture()
    for p in points[1:]:
        for s in p["host"]:
            key = "host_" + resource + "_pressure"
            s["optional"][key] = s["optional"][key].replace("total=0", "total=300000", 1)
    result = run(points, pressure_role=role)
    assert result["maximum_observed_some_percent"][resource] == pytest.approx(15.)
    if role == "eligibility":
        assert not result["sampled_pressure_policy_satisfied"]
        assert result["failures"][resource + "_pressure_bound_exceeded"] == 1
        assert "sampled_pressure_evidence_satisfied" not in result
    else:
        assert result["sampled_pressure_evidence_satisfied"]
        assert result["diagnostic_threshold_exceedances"][resource + "_pressure_bound_exceeded"] == 1
        assert result["failures"] == {}
        assert not result["diagnostic_thresholds_satisfied"]
        assert not result["pressure_thresholds_used_for_eligibility"]
        assert "sampled_pressure_policy_satisfied" not in result


@pytest.mark.parametrize("change", ["read_error", "boot", "group", "cpu_set", "missing",
    "counter", "overlap", "brackets", "gap", "late", "early", "malformed"])
@pytest.mark.parametrize("role", ["eligibility", "diagnostic_only"])
def test_incomplete_or_changed_evidence_rejected(change, role):
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
    elif change == "late": return assert_failed(run(points, launch_ns=10_000_000_000, pressure_role=role))
    elif change == "early": return assert_failed(run(points, end_ns=15_000_000_000, pressure_role=role))
    else: s["optional"]["host_cpu_pressure"] = "garbage"
    assert_failed(run(points, pressure_role=role))


def assert_failed(result):
    key = ("sampled_pressure_evidence_satisfied" if result["schema"].endswith("_v2")
           else "sampled_pressure_policy_satisfied")
    assert result[key] is False


@pytest.mark.parametrize("points", [[], [{}]])
def test_absent_evidence_never_passes(points):
    assert_failed(run(points))


@pytest.mark.parametrize("options", [dict(maximum_period_s=0), dict(limits=dict(cpu=1)),
    dict(limits=dict(cpu=float("nan"), io=1, memory=1)), dict(launch_ns=True),
    dict(end_ns=1), dict(limits=dict(cpu=101, io=1, memory=1))])
def test_invalid_bounds_rejected(options):
    with pytest.raises(ValueError):
        run(fixture(), **options)


@pytest.mark.parametrize("policy,expected", [
    (dict(schema="threadripper_environment_policy_v1"), "eligibility"),
    (dict(schema="threadripper_environment_policy_v2", native_pressure_role="diagnostic_only"), "diagnostic_only"),
    (dict(schema="threadripper_environment_policy_v2"), None),
    (dict(schema="threadripper_environment_policy_v1", native_pressure_role="diagnostic_only"), None),
    (dict(schema="threadripper_environment_policy_v2", native_pressure_role="eligibility"), None),
    (dict(schema="unknown"), None)])
def test_versioned_policy_role_is_explicit(policy, expected):
    if expected is None:
        with pytest.raises(ValueError):
            native_pressure_role(policy)
    else:
        assert native_pressure_role(policy) == expected


@pytest.mark.parametrize("role", [None, True, [], "ignore_evidence", ""])
def test_unknown_pressure_roles_rejected(role):
    with pytest.raises(ValueError):
        run(fixture(), pressure_role=role)


def long_boundary_points():
    points = [fixture()[0], fixture()[-1]]
    for sample in points[-1]['host']:
        sample['started_monotonic_ns'] += 96 * 10**9
        sample['finished_monotonic_ns'] += 96 * 10**9
    return points


def test_boundary_pressure_is_not_periodic_coverage():
    points = long_boundary_points()
    original = deepcopy(points)
    periodic = run(points, end_ns=109 * 10**9, pressure_role='diagnostic_only')
    assert periodic['failures'] == {'pressure_sample_period_exceeded': 1}
    assert not periodic['sampled_pressure_evidence_satisfied']
    result = run(points, end_ns=109 * 10**9, pressure_role='diagnostic_only', boundary_only=True)
    assert result['schema'] == 'threadripper_boundary_pressure_review_v1'
    assert result['sampled_pressure_evidence_satisfied'] and result['native_interval_bracketed']
    assert result['points'] == 2 and result['intervals'] == 1
    assert result['maximum_observed_period_s'] == 100.
    assert not result['periodic_pressure_cadence_checked']
    assert result['bounds']['maximum_point_duration_s'] == 3.
    assert 'maximum_period_s' not in result['bounds']
    assert 'maximum_observed_some_percent' not in result
    assert not result['scientific_timings_admitted']
    assert points == original


@pytest.mark.parametrize('resource', ['cpu', 'memory', 'io'])
def test_boundary_stalls_are_retained_interval_average_diagnostics(resource):
    points = long_boundary_points()
    for sample in points[-1]['host']:
        key = 'host_' + resource + '_pressure'
        sample['optional'][key] = sample['optional'][key].replace('total=0', 'total=15000000', 1)
    result = run(points, end_ns=109 * 10**9, pressure_role='diagnostic_only', boundary_only=True)
    assert result['interval_average_some_percent'][resource] == pytest.approx(15.)
    assert result['sampled_pressure_evidence_satisfied']
    assert result['diagnostic_threshold_exceedances'][resource + '_pressure_bound_exceeded'] == 1
    assert not result['diagnostic_thresholds_satisfied']


@pytest.mark.parametrize('defect', ['missing', 'extra', 'read_error', 'boot', 'group', 'counter', 'duration', 'late', 'early'])
def test_boundary_mode_keeps_evidence_integrity_and_exact_boundaries(defect):
    points = long_boundary_points()
    options = dict(end_ns=109 * 10**9, pressure_role='diagnostic_only', boundary_only=True)
    if defect == 'missing': points.pop()
    elif defect == 'extra': points.insert(1, fixture()[1])
    elif defect == 'read_error': points[-1]['host'][0]['errors'].append(dict(field='pressure', type='OSError'))
    elif defect == 'boot': points[-1]['host'][0]['raw']['boot_id'] = 'other'
    elif defect == 'group': points[-1]['host'][0]['raw']['cgroup_membership'] = '0::/other\n'
    elif defect == 'counter':
        key = 'host_io_pressure'
        points[0]['host'][0]['optional'][key] = points[0]['host'][0]['optional'][key].replace('total=0', 'total=1', 1)
    elif defect == 'duration':
        points[-1]['host'][-1]['started_monotonic_ns'] += 4 * 10**9
        points[-1]['host'][-1]['finished_monotonic_ns'] += 4 * 10**9
    elif defect == 'late': options['launch_ns'] = 10 * 10**9
    else: options['end_ns'] = 111 * 10**9
    assert not run(points, **options)['sampled_pressure_evidence_satisfied']


@pytest.mark.parametrize('options', [dict(boundary_only=True), dict(boundary_only=1),
                                    dict(boundary_only=None), dict(boundary_only='boundary')])
def test_boundary_pressure_cannot_relax_eligibility_policy(options):
    with pytest.raises(ValueError, match='Boundary pressure'):
        run(long_boundary_points(), **options)
