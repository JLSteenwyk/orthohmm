import json
from pathlib import Path

import pytest

from benchmark_tools.results import review_shared_prenative_failure_20261004 as reviewer
from benchmark_tools import threadripper_environment_worker as worker


@pytest.fixture
def case(tmp_path):
    run = dict(index=17, method='orthofinder_3_1_5_full', proteomes=12, repeat=1,
        measurement_directory=str(tmp_path / 'measurement'))
    request = dict(index=17, job_id=42, execution_scope=reviewer.executor.SHARED_SCOPE)
    documents = dict(
        preflight=dict(index=17, job_id=42, decision='failed', whole_run_observer_ready=False,
            execution_scope=reviewer.executor.SHARED_SCOPE, foreign_cpu_used_for_eligibility=False,
            uncontended_timing=False),
        detail=dict(status='failed', error_type='ValueError', native_release_authorized=False,
            error='Frozen input/source identity changed: ' + str(tmp_path / 'measurement/host_processes.jsonl'),
            snapshots=[dict(boot_id='boot', started_monotonic_s=25., finished_monotonic_s=28.),
                dict(boot_id='boot', started_monotonic_s=31., finished_monotonic_s=35.)],
            available_memory_bytes=[256 * 1024**3, 255 * 1024**3],
            process_review=dict(cpu_diagnostic=dict(sum_observed_foreign_average_cores=42.))),
        release=dict(status='environment_release_failed', error_type='ValueError',
            error='Execution evidence identity or decision differs'),
        wrapper=dict(status='verified_wrapper_failed', error_type='ValueError',
            error='Execution evidence identity or decision differs', scientific_results_admitted=False),
        result=dict(index=17, job_id=42, status='executor_failed',
            error='Environmental worker was not successfully joined at release', automatic_retry=False,
            next_submission_authorized=False, scientific_timings_admitted=False),
        lifecycle=dict(terminal=True, exit_code=1, status='failed_or_cancelled'),
        go=dict(abort=True), aborted=dict(status='observer_did_not_release_native'),
        policy=dict(minimum_available_memory_bytes=128 * 1024**3))
    for side in ('before', 'after'):
        documents['wrapper'][side] = dict(runtime=dict(status='runtime_and_lookup_checked',
            scientific_execution_authorized=False,
            lookup={name:dict(status='native_python_lookup_matches') for name in ('orthohmm', 'orthofinder')}))
    rows = [dict(index=i, observer_pid=123, interval=None if i == 0 else {},
        snapshot=dict(boot_id='boot', started_monotonic_s=float(i * 30))) for i in range(2)]
    return run, request, documents, rows


def test_abort_has_no_native_measurement_or_progression(case):
    result = reviewer.assess(*case, [])
    assert result['native_outcome'] == 'not_started'
    assert result['resources'] is None
    assert result['status'] == 'pre_native_infrastructure_failure_reviewed'
    for field in ('comparative_timing_eligible', 'automatic_retry',
                  'next_submission_authorized', 'scientific_timings_admitted', 'uncontended_timing'):
        assert result[field] is False


@pytest.mark.parametrize('change', ['identity', 'job', 'gate', 'artifact', 'error', 'wrapper',
    'lifecycle', 'samples', 'observer', 'boot', 'overlap', 'capacity', 'runtime', 'authorization'])
def test_refuses_borrowed_or_contradictory_failure(case, change):
    run, request, docs, rows = case
    present = []
    if change == 'identity': run['repeat'] = 2
    elif change == 'job': docs['preflight']['job_id'] = 43
    elif change == 'gate': docs['go'] = {'go': True}
    elif change == 'artifact': present = ['native.log']
    elif change == 'error': docs['detail']['error'] = 'Insufficient safe available memory'
    elif change == 'wrapper': docs['wrapper']['status'] = 'command_exited_zero'
    elif change == 'lifecycle': docs['lifecycle']['terminal'] = False
    elif change == 'samples': rows.pop()
    elif change == 'observer': rows[1]['observer_pid'] = 124
    elif change == 'boot': rows[1]['snapshot']['boot_id'] = 'another'
    elif change == 'overlap': docs['detail']['snapshots'][1]['finished_monotonic_s'] = 29.
    elif change == 'capacity': docs['detail']['available_memory_bytes'][0] = 64 * 1024**3
    elif change == 'runtime': docs['wrapper']['after']['runtime']['lookup']['orthofinder']['status'] = 'mismatch'
    else: docs['result']['next_submission_authorized'] = True
    with pytest.raises(ValueError):
        reviewer.assess(run, request, docs, rows, present)


def test_shared_initial_copy_survives_append_while_old_stream_pin_fails(tmp_path, monkeypatch):
    ready = dict(pid=456, cgroup='scope')
    row = dict(pid=123)
    sample = dict(index=0, interval=None, observer_pid=123,
        snapshot=dict(processes=[row, dict(pid=456)], errors=[]))
    (tmp_path / 'ready.json').write_text(json.dumps(ready))
    stream = tmp_path / 'host_processes.jsonl'
    stream.write_text(json.dumps(sample) + '\n')
    monkeypatch.setattr(worker, 'job_scope', lambda *args: 'scope')
    monkeypatch.setattr(worker, 'membership', lambda *args: 'scope')
    monkeypatch.setattr(worker, 'scoped_path', lambda *args: 'scope')
    monkeypatch.setattr(worker, 'live_identity', lambda *args: None)
    pins = worker.collector_ready(tmp_path, 'scope', 42, shared_host=True)
    old_stream_pin = worker.record(stream)
    worker.check(pins[1])
    with stream.open('a') as handle:
        handle.write(json.dumps(dict(sample, index=1, interval={})) + '\n')
    with pytest.raises(ValueError, match='Frozen input/source identity changed'):
        worker.check(old_stream_pin)
    worker.check(pins[1])
    assert worker.collector_ready(tmp_path, 'scope', 42, shared_host=True) == pins
    assert Path(pins[1]['path']).name == 'preflight_initial_process_sample.json'


def test_actual_receipt_preserves_failure_and_null_resources():
    root = Path(__file__).resolve().parents[2]
    report = json.loads((root / 'benchmark_tools/results/threadripper_shared_prenative_failure_22413.json').read_text())
    assert report['job_id'] == 22413 and report['index'] == 17
    assert report['scheduler_state'] == 'FAILED' and report['scheduler_exit_code'] == '1:0'
    assert report['native_outcome'] == 'not_started' and report['resources'] is None
    assert report['first_process_interval_seconds'] == 30.00103936693631
    assert report['available_memory_bytes'] == [332716515328, 332481789952]
    assert report['next_submission_authorized'] is False
    assert report['terminal_accounting']['returncode'] == 0
    assert '22413|FAILED|1:0|2026-10-04T01:19:39|2026-10-04T01:22:42' in report['terminal_accounting']['stdout']
    reviewer.executor.check(report['source'])
