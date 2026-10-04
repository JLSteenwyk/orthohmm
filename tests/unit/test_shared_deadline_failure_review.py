import copy
import json
from pathlib import Path

import pytest

from benchmark_tools.results import review_shared_deadline_failure_20261004 as audit


@pytest.fixture
def retained():
    root = Path(__file__).resolve().parents[2]
    work = root / 'benchmarks/results/threadripper_scaling_v1'
    measurement = work / 'run_20/measurement'
    session = work / 'sessions/run_20'
    if not (session / 'result.json').exists():
        pytest.skip('Retained raw failure evidence is local-only')
    read = lambda path: json.loads(path.read_text())
    run = read(root / 'benchmark_tools/results/threadripper_private_commands_20260928.json')['runs'][20]
    request = read(root / 'benchmarks/work/threadripper_shared_execution_20261003/request_run_20.json')
    paths = dict(preflight=session / 'environment_preflight.json', detail=session / 'environment_worker_evidence.json',
        lifecycle=session / 'environment_worker_lifecycle.json', result=session / 'result.json',
        wrapper=measurement.parent / 'verification.json', release=measurement / 'environment_release.json',
        go=measurement / 'go.json', aborted=measurement / 'aborted_before_native.json',
        marker=measurement / 'environment_review_requested.json', ready=measurement / 'ready.json',
        initial=measurement / 'preflight_initial_process_sample.json',
        policy=root / 'benchmark_tools/results/threadripper_shared_environment_policy_prenative_20261004.json')
    documents = {key: read(path) for key, path in paths.items()}
    rows = [json.loads(line) for line in (measurement / 'host_processes.jsonl').read_text().splitlines()]
    return run, request, documents, rows, []


def test_actual_timeout_has_no_native_resources(retained):
    value = audit.assess(*retained)
    assert value['resources'] is None and value['native_outcome'] == 'not_started'
    assert value['request_to_review_start_seconds'] == 11.799922912
    assert value['request_to_response_seconds'] == 21.377395496
    assert value['response_deadline_seconds'] == 20
    assert not value['comparative_timing_eligible'] and not value['automatic_retry']
    assert not value['next_submission_authorized']


@pytest.mark.parametrize('key,field,value', [
    ('preflight', 'decision', 'passed'), ('preflight', 'job_id', 22415),
    ('preflight', 'publication_deadline_expired', False),
    ('preflight', 'observation_finished_unix_ns', 1791098320000000000),
    ('marker', 'index', 19), ('marker', 'wait_seconds', 30),
    ('marker', 'requested_unix_ns', True), ('detail', 'status', 'failed'),
    ('detail', 'error_type', 'TimeoutError'), ('detail', 'error', 'unrelated failure'),
    ('release', 'error_type', 'ValueError'), ('wrapper', 'status', 'verified_wrapper_passed'),
    ('lifecycle', 'terminal', False), ('result', 'automatic_retry', True),
    ('go', 'abort', False), ('ready', 'pid', 999),
    ('detail', 'available_memory_bytes', [1, 1]),
])
def test_contradictory_failure_is_refused(retained, key, field, value):
    args = copy.deepcopy(retained)
    args[2][key][field] = value
    with pytest.raises(ValueError):
        audit.assess(*args)


def test_native_artifacts_refuse_pre_native_classification(retained):
    with pytest.raises(ValueError, match='Native artifacts'):
        audit.assess(*retained[:4], ['native.log'])


def test_changed_initial_sample_is_refused(retained):
    args = copy.deepcopy(retained)
    args[2]['initial']['index'] = 1
    with pytest.raises(ValueError, match='initial sample'):
        audit.assess(*args)


def test_empty_preparation_directory_is_not_native_output(tmp_path):
    output = tmp_path / 'output'
    output.mkdir()
    run = dict(measurement_directory=str(tmp_path),
        configuration=dict(output=str(output), metrics=str(tmp_path / 'metrics.json')))
    assert audit.native_artifacts(run) == []
    (output / '.unexpected').touch()
    assert audit.native_artifacts(run) == [str(output)]


@pytest.mark.parametrize('name', ['native.log', 'done.json', 'lineage_report.json', 'point_0001.json', 'metrics.json'])
def test_every_native_marker_is_detected(tmp_path, name):
    run = dict(measurement_directory=str(tmp_path),
        configuration=dict(output=str(tmp_path / 'output'), metrics=str(tmp_path / 'metrics.json')))
    (tmp_path / name).touch()
    assert audit.native_artifacts(run) == [str(tmp_path / name)]


def test_output_symlink_is_not_an_empty_preparation_directory(tmp_path):
    empty = tmp_path / 'empty'
    empty.mkdir()
    output = tmp_path / 'output'
    output.symlink_to(empty)
    run = dict(measurement_directory=str(tmp_path),
        configuration=dict(output=str(output), metrics=str(tmp_path / 'metrics.json')))
    assert audit.native_artifacts(run) == [str(output)]
