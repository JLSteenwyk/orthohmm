import json
from pathlib import Path

import pytest

from benchmark_tools.prepare_scaling_inputs import planned_runs
from benchmark_tools.results import continue_shared_prepared_panel_20261004 as launcher
from benchmark_tools.results import review_shared_prepared_panel_20261004 as reviewer
from benchmark_tools.results import snapshot_shared_prepared_panel_20261004 as snapshot


@pytest.mark.parametrize('index', [True, 0, 1, 17, 18, 19, 20, 28, -1])
def test_continuation_cannot_retry_or_skip_unreviewed_prefix(index):
    with pytest.raises(ValueError): launcher.history(index, {})


def test_history_uses_retained_resolutions_not_fake_success(monkeypatch):
    captured = []
    monkeypatch.setattr(launcher.executor, 'record', lambda path: dict(path=str(path)))
    def read(pin):
        index = int(Path(pin['path']).parent.name.removeprefix('review_run_'))
        return dict(status='shared_attempt_independently_reviewed', index=index,
            shared_host_resources_reviewed=True, original_environment_protocol_passed=True,
            session=dict(path=f'session_{index}'))
    monkeypatch.setattr(launcher.executor, 'read', read)
    def bind(plan, refs, **kwargs):
        captured.append((plan, refs, kwargs))
        return dict(progress=dict(status='next_identity_requires_preflight', index=21))
    monkeypatch.setattr(launcher.executor, 'bind', bind)
    refs, _ = launcher.history(21, dict(sha256='plan'))
    assert len(refs) == 21
    assert refs[0] == dict(path=str(launcher.RESOLVED))
    assert refs[17] == dict(path=str(launcher.ABORT_RESOLVED[17]))
    assert refs[20] == dict(path=str(launcher.ABORT_RESOLVED[20]))
    assert refs[18:20] == [dict(path=f'session_{i}') for i in (18, 19)]
    assert refs[1:17] == [dict(path=f'session_{i}') for i in range(1, 17)]
    assert captured[0][2]['allocation_mode'] == 'shared'


@pytest.mark.parametrize('state,index', [('unresolved_prior_attempt', 21),
    ('next_identity_requires_preflight', 17), ('all_attempts_reviewed', None)])
def test_wrong_next_history_position_never_launches(monkeypatch, state, index):
    monkeypatch.setattr(launcher.executor, 'record', lambda path: dict(path=str(path)))
    monkeypatch.setattr(launcher.executor, 'read', lambda pin: dict(
        status='shared_attempt_independently_reviewed', index=int(Path(pin['path']).parent.name[-2:]),
        shared_host_resources_reviewed=True, original_environment_protocol_passed=True, session=pin))
    monkeypatch.setattr(launcher.executor, 'bind', lambda *a, **k: dict(progress=dict(status=state,index=index)))
    with pytest.raises(ValueError): launcher.history(21, {})


def test_reviewer_keeps_existing_raw_checks_but_uses_new_bound_launcher():
    current = Path(reviewer.__file__).read_text()
    previous = Path(reviewer.__file__).with_name('review_shared_prenative_panel_20261004.py').read_text()
    expected = previous.replace('continue_shared_prenative_panel_20261004', 'continue_shared_prepared_panel_20261004')
    expected = expected.replace('not 0 <= index', 'not 21 <= index')
    expected = expected.replace("        ready = executor.read(executor.record(measurement / 'ready.json'))",
        "        handoff = launch.observer_handoff(session, measurement, request_ref, policy_ref, index, job, terminal=True)\n"
        "        handoff_ref = launch.save(directory / 'prepared_handoff.json', handoff)\n"
        "        ready = executor.read(executor.record(measurement / 'ready.json'))")
    expected = expected.replace('environment=[preflight_ref, retained_ref, env_ref, policy_ref]',
        "environment=[preflight_ref, retained_ref, env_ref, policy_ref, handoff_ref, *handoff['evidence']]")
    expected = expected.replace('terminal_accounting=accounting_ref,', 'terminal_accounting=accounting_ref, prepared_handoff=handoff_ref,')
    assert current == expected
    assert reviewer.launch is launcher and snapshot.launcher is launcher


def launch_fixture(tmp_path, monkeypatch, *, existing=False, corrupt_held=False):
    work = tmp_path / 'work'
    work.mkdir()
    monkeypatch.setattr(launcher, 'WORK', work)
    prep_path = work / 'preparation.json'
    monkeypatch.setattr(launcher, 'PREPARATION', prep_path)
    refs = {}
    for key in ('lookup', 'binding', 'plan', 'recipe', 'policy', 'readiness', 'resources'):
        path = work / (key + '.json')
        path.write_text('{}\n')
        refs[key] = launcher.executor.record(path)
    prep_path.write_text(json.dumps(refs))
    run = dict(planned_runs()[21], measurement_directory=str(tmp_path / 'measurements/run_21/measurement'))
    plan = dict(runs=[{}]*21 + [run] + [{}]*5)
    monkeypatch.setattr(launcher, 'load', lambda: ({}, {}, plan, refs['plan']))
    prefix = [dict(path='resolved_0'), dict(path='resolved_17')]
    monkeypatch.setattr(launcher, 'history', lambda index, plan_ref: (prefix, dict(progress=dict(index=index))))
    monkeypatch.setattr(launcher.executor, 'select', lambda request, root, job:
        (run, {}, dict(progress=dict(index=21)), []))
    calls = []
    def command(label, argv):
        calls.append((label, argv))
        if label == 'submission_21': return '999\n'
        if label == 'held_controller_21':
            return (f'JobId=999 JobState=PENDING Reason=JobHeldUser NumCPUs={32 if corrupt_held else 64} '
                f'NumTasks=1 CPUs/Task=64 MinMemoryNode=128G TimeLimit={launcher.executor.TIME_LIMIT} '
                f'Requeue=0 Command={launcher.SCRIPT} WorkDir={launcher.ROOT}')
        if label == 'pinned_controller_21':
            return 'Comment=' + launcher.executor.record(work / 'request_run_21.json')['sha256']
        return ''
    monkeypatch.setattr(launcher.first, 'command', command)
    if existing: (work / 'request_run_21.json').write_text('{}')
    return work, refs, prefix, calls


def test_launch_binds_new_lookup_preserved_prefix_and_no_retry(tmp_path, monkeypatch):
    work, refs, prefix, calls = launch_fixture(tmp_path, monkeypatch)
    result = launcher.launch(21)
    request = json.loads((work / 'request_run_21.json').read_text())
    assert request['runtime_lookup'] == refs['lookup']
    assert request['lookup_sha256'] == launcher.LOOKUP_SHA
    assert request['history'] == prefix and request['index'] == 21 and request['job_id'] == 999
    assert result['method'] == 'orthohmm_high_sensitivity' and result['repeat'] == 2
    assert result['automatic_retry'] is False and result['scientific_timings_admitted'] is False
    assert [label for label, _ in calls] == ['submission_21', 'held_controller_21',
        'request_comment_21', 'pinned_controller_21', 'release_21']
    assert '--hold' in calls[0][1]


def test_existing_attempt_is_not_submitted_again(tmp_path, monkeypatch):
    _, _, _, calls = launch_fixture(tmp_path, monkeypatch, existing=True)
    with pytest.raises(FileExistsError): launcher.launch(21)
    assert calls == []


def test_wrong_held_limits_remain_held_without_request_or_release(tmp_path, monkeypatch):
    work, _, _, calls = launch_fixture(tmp_path, monkeypatch, corrupt_held=True)
    with pytest.raises(ValueError): launcher.launch(21)
    assert [label for label, _ in calls] == ['submission_21', 'held_controller_21']
    assert not (work / 'request_run_21.json').exists()


def test_preparation_refuses_overwrite_before_load(tmp_path, monkeypatch):
    path = tmp_path / 'existing.json'
    path.write_text('{}')
    monkeypatch.setattr(launcher, 'PREPARATION', path)
    monkeypatch.setattr(launcher, 'load', lambda: pytest.fail('Do not rebind an existing preparation'))
    with pytest.raises(FileExistsError): launcher.prepare()


def handoff_fixture(tmp_path):
    session, measurement = tmp_path / 'session', tmp_path / 'measurement'
    session.mkdir()
    measurement.mkdir()
    request = launcher.save(tmp_path / 'request.json', dict(job_id=999, index=21))
    policy = launcher.save(tmp_path / 'policy.json', {})
    started = dict(status='started', request=request, policy=policy, slurm_job_id='999',
        automatic_retry=False, scientific_timings_admitted=False, pid=101, started_unix_ns=1_000_000_000)
    launcher.save(session / 'environment_worker_started.json', started)
    prepared = launcher.save(session / 'environment_worker_prepared.json', dict(
        schema='threadripper_environment_worker_prepared_v1',
        status='prepared_waiting_for_release_request', request=request, policy=policy,
        job_id=999, index=21, pid=101, job_scope='/slurm/job_999', boot_id='boot',
        prepared_unix_ns=2_000_000_000, native_release_authorized=False, scientific_timings_admitted=False))
    launcher.save(measurement / 'environment_review_requested.json', dict(job_id=999, index=21,
        request=request, review_path=str(session / 'environment_preflight.json'), wait_seconds=20,
        native_released=False, requested_unix_ns=3_000_000_000))
    preflight = launcher.save(session / 'environment_preflight.json', dict(job_id=999, index=21,
        decision='passed', environment_policy=policy, whole_run_observer_ready=True,
        execution_scope=launcher.executor.SHARED_SCOPE, uncontended_timing=False, boot_id='boot',
        observation_started_unix_ns=4_000_000_000, observation_finished_unix_ns=7_000_000_000,
        evidence=[prepared]))
    launcher.save(measurement / 'environment_release.json', dict(status='environment_review_bound',
        review=preflight, observational_validity_independently_established=False,
        scientific_timings_admitted=False, checked_unix_ns=8_000_000_000))
    launcher.save(measurement / 'go.json', dict(go=True))
    launcher.save(measurement / 'ready.json', dict(cgroup='0::/slurm/job_999/step_0/user/task\n'))
    log = session / 'environment_worker.log'
    log.write_text('owned worker joined\n')
    lifecycle = dict(started, status='completed', prepared=prepared, exit_code=0, terminal=True,
        parent_error_type=None, finished_unix_ns=9_000_000_000, log=launcher.executor.record(log))
    launcher.save(session / 'environment_worker_lifecycle.json', lifecycle)
    return session, measurement, request, policy


def rewrite(path, **fields):
    data = json.loads(path.read_text())
    data.update(fields)
    path.write_text(json.dumps(data))


@pytest.mark.parametrize('terminal', [False, True])
def test_actual_receipt_order_replays_without_live_observer(tmp_path, terminal):
    session, measurement, request, policy = handoff_fixture(tmp_path)
    result = launcher.observer_handoff(session, measurement, request, policy, 21, 999, terminal=terminal)
    assert result['prepared_before_release_request'] is True
    assert result['worker_terminal_reviewed'] is terminal
    assert result['preparation_to_request_seconds'] == 1
    assert result['request_to_bound_review_seconds'] == 5
    assert len(result['evidence']) == (8 if terminal else 7)


@pytest.mark.parametrize('file,fields', [
    ('environment_worker_prepared.json', dict(pid=102)),
    ('environment_worker_prepared.json', dict(job_id=998)),
    ('environment_worker_prepared.json', dict(index=20)),
    ('environment_worker_prepared.json', dict(boot_id='other')),
    ('environment_worker_prepared.json', dict(job_scope='/slurm/job_998')),
    ('environment_worker_prepared.json', dict(request={})),
    ('environment_worker_prepared.json', dict(policy={})),
    ('environment_worker_prepared.json', dict(native_release_authorized=True)),
    ('environment_worker_prepared.json', dict(prepared_unix_ns=5_000_000_000)),
    ('environment_worker_prepared.json', dict(prepared_unix_ns=True)),
    ('environment_worker_started.json', dict(pid=True)),
    ('environment_worker_started.json', dict(slurm_job_id='998')),
    ('environment_worker_lifecycle.json', dict(prepared={})),
    ('environment_worker_lifecycle.json', dict(status='failed_or_cancelled')),
    ('environment_worker_lifecycle.json', dict(terminal=False)),
    ('environment_worker_lifecycle.json', dict(exit_code=1)),
    ('environment_worker_lifecycle.json', dict(cleanup_requested=True)),
    ('environment_worker_lifecycle.json', dict(kill_requested=True)),
    ('environment_worker_lifecycle.json', dict(finished_unix_ns=6_000_000_000)),
])
def test_borrowed_stale_or_aborted_preparation_is_not_admitted(tmp_path, file, fields):
    session, measurement, request, policy = handoff_fixture(tmp_path)
    rewrite(session / file, **fields)
    with pytest.raises(ValueError):
        launcher.observer_handoff(session, measurement, request, policy, 21, 999, terminal=True)


@pytest.mark.parametrize('file,fields', [
    ('environment_review_requested.json', dict(wait_seconds=21)),
    ('environment_review_requested.json', dict(native_released=True)),
    ('environment_review_requested.json', dict(request={})),
    ('environment_review_requested.json', dict(requested_unix_ns=8_000_000_000)),
    ('environment_release.json', dict(status='environment_release_failed')),
    ('environment_release.json', dict(review={})),
    ('environment_release.json', dict(checked_unix_ns=6_000_000_000)),
    ('go.json', dict(go=False)),
    ('go.json', dict(go=1)),
    ('go.json', dict(abort=True)),
])
def test_failed_mismatched_or_reordered_release_is_not_admitted(tmp_path, file, fields):
    session, measurement, request, policy = handoff_fixture(tmp_path)
    rewrite(measurement / file, **fields)
    with pytest.raises(ValueError):
        launcher.observer_handoff(session, measurement, request, policy, 21, 999)


def test_live_snapshot_does_not_require_post_measurement_lifecycle(tmp_path):
    session, measurement, request, policy = handoff_fixture(tmp_path)
    (session / 'environment_worker_lifecycle.json').unlink()
    assert launcher.observer_handoff(session, measurement, request, policy, 21, 999)['worker_terminal_reviewed'] is False
    with pytest.raises(FileNotFoundError):
        launcher.observer_handoff(session, measurement, request, policy, 21, 999, terminal=True)


def test_preflight_must_pin_prepared_receipt(tmp_path):
    session, measurement, request, policy = handoff_fixture(tmp_path)
    rewrite(session / 'environment_preflight.json', evidence=[])
    preflight = launcher.executor.record(session / 'environment_preflight.json')
    rewrite(measurement / 'environment_release.json', review=preflight)
    with pytest.raises(ValueError, match='bind the preparation'):
        launcher.observer_handoff(session, measurement, request, policy, 21, 999)


@pytest.mark.parametrize('observation_end,checked', [
    (23_000_000_000, 24_000_000_000), (7_000_000_000, 125_000_000_000)])
def test_response_and_freshness_bounds_are_not_relaxed(tmp_path, observation_end, checked):
    session, measurement, request, policy = handoff_fixture(tmp_path)
    rewrite(session / 'environment_preflight.json', observation_finished_unix_ns=observation_end)
    preflight = launcher.executor.record(session / 'environment_preflight.json')
    rewrite(measurement / 'environment_release.json', review=preflight, checked_unix_ns=checked)
    with pytest.raises(ValueError, match='unchanged'):
        launcher.observer_handoff(session, measurement, request, policy, 21, 999)


def test_checked_time_after_response_deadline_is_allowed_within_freshness(tmp_path):
    session, measurement, request, policy = handoff_fixture(tmp_path)
    rewrite(measurement / 'environment_release.json', checked_unix_ns=24_000_000_000)
    assert launcher.observer_handoff(session, measurement, request, policy, 21, 999)['native_gate_go'] is True
