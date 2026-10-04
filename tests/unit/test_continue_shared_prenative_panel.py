import json
from pathlib import Path

import pytest

from benchmark_tools.prepare_scaling_inputs import planned_runs
from benchmark_tools.results import continue_shared_prenative_panel_20261004 as launcher
from benchmark_tools.results import review_shared_prenative_panel_20261004 as reviewer
from benchmark_tools.results import snapshot_shared_prenative_panel_20261004 as snapshot


@pytest.mark.parametrize('index', [True, 0, 1, 17, 28, -1])
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
        return dict(progress=dict(status='next_identity_requires_preflight', index=18))
    monkeypatch.setattr(launcher.executor, 'bind', bind)
    refs, _ = launcher.history(18, dict(sha256='plan'))
    assert len(refs) == 18
    assert refs[0] == dict(path=str(launcher.RESOLVED))
    assert refs[17] == dict(path=str(launcher.PRENATIVE_RESOLVED))
    assert refs[1:17] == [dict(path=f'session_{i}') for i in range(1, 17)]
    assert captured[0][2]['allocation_mode'] == 'shared'


@pytest.mark.parametrize('state,index', [('unresolved_prior_attempt', 18),
    ('next_identity_requires_preflight', 17), ('all_attempts_reviewed', None)])
def test_wrong_next_history_position_never_launches(monkeypatch, state, index):
    monkeypatch.setattr(launcher.executor, 'record', lambda path: dict(path=str(path)))
    monkeypatch.setattr(launcher.executor, 'read', lambda pin: dict(
        status='shared_attempt_independently_reviewed', index=int(Path(pin['path']).parent.name[-2:]),
        shared_host_resources_reviewed=True, original_environment_protocol_passed=True, session=pin))
    monkeypatch.setattr(launcher.executor, 'bind', lambda *a, **k: dict(progress=dict(status=state,index=index)))
    with pytest.raises(ValueError): launcher.history(18, {})


def test_reviewer_keeps_existing_raw_checks_but_uses_new_bound_launcher():
    previous = Path(reviewer.__file__).with_name('review_shared_threadripper_repaired_20261003.py')
    assert Path(reviewer.__file__).read_text() == previous.read_text().replace(
        'continue_shared_threadripper_panel_20261003', 'continue_shared_prenative_panel_20261004')
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
    run = dict(planned_runs()[18], measurement_directory=str(tmp_path / 'measurements/run_18/measurement'))
    plan = dict(runs=[{}]*18 + [run] + [{}]*8)
    monkeypatch.setattr(launcher, 'load', lambda: ({}, {}, plan, refs['plan']))
    prefix = [dict(path='resolved_0'), dict(path='resolved_17')]
    monkeypatch.setattr(launcher, 'history', lambda index, plan_ref: (prefix, dict(progress=dict(index=index))))
    monkeypatch.setattr(launcher.executor, 'select', lambda request, root, job:
        (run, {}, dict(progress=dict(index=18)), []))
    calls = []
    def command(label, argv):
        calls.append((label, argv))
        if label == 'submission_18': return '999\n'
        if label == 'held_controller_18':
            return (f'JobId=999 JobState=PENDING Reason=JobHeldUser NumCPUs={32 if corrupt_held else 64} '
                f'NumTasks=1 CPUs/Task=64 MinMemoryNode=128G TimeLimit={launcher.executor.TIME_LIMIT} '
                f'Requeue=0 Command={launcher.SCRIPT} WorkDir={launcher.ROOT}')
        if label == 'pinned_controller_18':
            return 'Comment=' + launcher.executor.record(work / 'request_run_18.json')['sha256']
        return ''
    monkeypatch.setattr(launcher.first, 'command', command)
    if existing: (work / 'request_run_18.json').write_text('{}')
    return work, refs, prefix, calls


def test_launch_binds_new_lookup_preserved_prefix_and_no_retry(tmp_path, monkeypatch):
    work, refs, prefix, calls = launch_fixture(tmp_path, monkeypatch)
    result = launcher.launch(18)
    request = json.loads((work / 'request_run_18.json').read_text())
    assert request['runtime_lookup'] == refs['lookup']
    assert request['lookup_sha256'] == launcher.LOOKUP_SHA
    assert request['history'] == prefix and request['index'] == 18 and request['job_id'] == 999
    assert result['method'] == 'orthofinder_3_1_5_full' and result['repeat'] == 2
    assert result['automatic_retry'] is False and result['scientific_timings_admitted'] is False
    assert [label for label, _ in calls] == ['submission_18', 'held_controller_18',
        'request_comment_18', 'pinned_controller_18', 'release_18']
    assert '--hold' in calls[0][1]


def test_existing_attempt_is_not_submitted_again(tmp_path, monkeypatch):
    _, _, _, calls = launch_fixture(tmp_path, monkeypatch, existing=True)
    with pytest.raises(FileExistsError): launcher.launch(18)
    assert calls == []


def test_wrong_held_limits_remain_held_without_request_or_release(tmp_path, monkeypatch):
    work, _, _, calls = launch_fixture(tmp_path, monkeypatch, corrupt_held=True)
    with pytest.raises(ValueError): launcher.launch(18)
    assert [label for label, _ in calls] == ['submission_18', 'held_controller_18']
    assert not (work / 'request_run_18.json').exists()


def test_preparation_refuses_overwrite_before_load(tmp_path, monkeypatch):
    path = tmp_path / 'existing.json'
    path.write_text('{}')
    monkeypatch.setattr(launcher, 'PREPARATION', path)
    monkeypatch.setattr(launcher, 'load', lambda: pytest.fail('Do not rebind an existing preparation'))
    with pytest.raises(FileExistsError): launcher.prepare()
