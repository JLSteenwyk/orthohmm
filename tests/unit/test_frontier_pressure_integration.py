from copy import deepcopy
import json

import pytest

from benchmark_tools import measure_native_frontier_step as module
from benchmark_tools import measure_frontier_boundary_step as boundary
from tests.unit.test_measure_native_frontier_step import extend
from tests.unit.test_measure_native_hierarchy_step import evidence, setup_measure
from tests.unit.test_probe_native_pressure import raw


def add_pressure(points):
    for index, point in enumerate(points):
        start = point['frontier']['root'][1]['finished_ns'] + 10
        hosts = [deepcopy(point['host'][1]) for _ in range(2)]
        for host, offset in zip(hosts, (0, 100)):
            host['started_monotonic_ns'] = start + offset
            host['finished_monotonic_ns'] = start + offset + 1
        point['native_pressure'] = dict(
            host=hosts, native_membership=point['native_membership'],
            scope=module.interval_point(point, 21816)['native_cpu_scope'],
            scope_identity=[32, 100], ticks=point['ticks'],
            pressure={resource: dict(raw=raw(index * 100, index * 40),
                totals=dict(some=index * 100, full=index * 40),
                started_ns=start + 20 + i * 10, finished_ns=start + 21 + i * 10)
                for i, resource in enumerate(('cpu', 'memory', 'io'))})


def test_whole_command_and_intervals_preserve_cpu_screen(evidence):
    points, done = evidence
    extend(points, done)
    before = module.evaluate(points, done, 21816)
    add_pressure(points)
    after = module.evaluate(points, done, 21816)
    assert all(after[key] == value for key, value in before.items())
    assert after['native_pressure_whole_command']['native_stall_usec']['cpu'] == {'some': 100, 'full': 40}
    assert len(after['native_pressure_intervals']) == 1
    assert not after['scientific_timings_admitted']


@pytest.mark.parametrize('fault', ['missing', 'scope', 'identity', 'before', 'after', 'boot', 'reset'])
def test_inconsistent_pressure_rejected(evidence, fault):
    points, done = evidence
    extend(points, done)
    add_pressure(points)
    pressure = points[1]['native_pressure']
    if fault == 'missing':
        del points[1]['native_pressure']
    elif fault == 'scope':
        pressure['native_membership'] = pressure['native_membership'].replace('step_0', 'step_9')
    elif fault == 'identity':
        pressure['scope_identity'][1] += 1
    elif fault == 'before':
        pressure['host'][0]['started_monotonic_ns'] = 0
    elif fault == 'after':
        pressure['host'][1]['finished_monotonic_ns'] = points[1]['host'][1]['finished_monotonic_ns'] + 1
    elif fault == 'boot':
        for host in pressure['host']:
            host['raw']['boot_id'] = 'different'
    else:
        points[0]['native_pressure']['pressure']['cpu'].update(raw=raw(200, 100), totals=dict(some=200, full=100))
    with pytest.raises(ValueError):
        module.evaluate(points, done, 21816)


@pytest.mark.parametrize('exit_code', [0, 7, 124])
@pytest.mark.parametrize('collector', [module, boundary])
def test_opt_in_worker_lifecycle(tmp_path, monkeypatch, evidence, exit_code, collector):
    import tests.unit.test_measure_native_hierarchy_step as tests
    points, done = evidence
    extend(points, done)
    add_pressure(points)
    monkeypatch.setattr(tests, 'module', collector)
    if collector is boundary:
        monkeypatch.setattr(collector, 'read_hierarchy', None, raising=False)
    directory = setup_measure(tmp_path, monkeypatch, evidence, exit_code)
    pending = iter(points)
    def read(*args, native_pressure=False):
        assert native_pressure
        return next(pending)
    monkeypatch.setattr(collector, 'read_frontier_point', read)
    result = collector.measure(['/usr/bin/true'], directory, 21816, 20, 96 * 1024**3, 60, 1., native_pressure=True)
    assert result['native']['exit_code'] == exit_code
    assert 'native_pressure_whole_command' in result['screening']
    assert (directory / 'release.json').exists()
    assert not result['scientific_timings_admitted']


def test_reader_retains_pressure_inside_outer_bracket(monkeypatch, evidence):
    points, done = evidence
    initial = deepcopy(points[0])
    extend(points, done)
    add_pressure(points)
    expected = points[0]
    monkeypatch.setattr(module, 'read_hierarchy', lambda *args: deepcopy(initial))
    monkeypatch.setattr(module, 'frontier_snapshot', lambda *args: deepcopy(expected['frontier']))
    monkeypatch.setattr(module, 'read_pressure', lambda *args: deepcopy(expected['native_pressure']))
    monkeypatch.setattr(module, 'host_snapshot', lambda: deepcopy(expected['host'][1]))
    read_text = module.Path.read_text
    monkeypatch.setattr(module.Path, 'read_text', lambda path, *a, **kw:
                        initial['native_membership'] if str(path) == '/proc/123/cgroup'
                        else read_text(path, *a, **kw))
    assert module.read_frontier_point(123, initial['native_membership'], 21816, native_pressure=True) == expected


def test_pressure_read_failure_retained(tmp_path, monkeypatch, evidence):
    points, done = evidence
    initial = deepcopy(points[0])
    extend(points, done)
    monkeypatch.setattr(module, 'read_hierarchy', lambda *a: deepcopy(initial))
    monkeypatch.setattr(module, 'frontier_snapshot', lambda *a: deepcopy(points[0]['frontier']))
    def fail(*args):
        raise OSError('pressure unavailable')
    monkeypatch.setattr(module, 'read_pressure', fail)
    output = tmp_path / 'failed.json'
    with pytest.raises(OSError, match='unavailable'):
        module.read_frontier_point(123, initial['native_membership'], 21816, output, native_pressure=True)
    saved = json.loads(output.read_text())
    assert saved['status'] == 'invalid_pressure_observation'
    assert not saved['scientific_timings_admitted']
    assert 'frontier' in saved['preceding_point']


def test_boundary_pressure_matches_periodic_without_interval_claim(evidence):
    points, done = evidence
    extend(points, done)
    add_pressure(points)
    result = boundary.evaluate(points, done, 21816)
    assert result['native_pressure_whole_command'] == module.evaluate(points, done, 21816)['native_pressure_whole_command']
    assert result['flagged_intervals'] is None
    assert not result['interval_screening_available']
    del points[1]['native_pressure']
    with pytest.raises(ValueError, match='missing'):
        boundary.evaluate(points, done, 21816)
