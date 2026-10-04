import hashlib
import json
from pathlib import Path

import pytest

from benchmark_tools.results import continue_shared_goal_reaffirmed_20261004 as current


def test_only_exact_authorized_addition_is_accepted():
    old = b'Goal\n  >\nExisting full scope\n'
    new = b'Goal\n  >\n' + current.ADDITION.encode('ascii') + b'Existing full scope\n'
    current.reaffirmed(old, new)
    for changed in (old, new + b'Extra change\n', new.replace(b'full scope', b'narrow scope')):
        with pytest.raises(ValueError):
            current.reaffirmed(old, changed)


@pytest.mark.parametrize('refs', [[], [{'path': 'other'}], [{'path': 'old'}, {'path': 'old'}]])
def test_pin_replacement_requires_exactly_one_document(refs):
    with pytest.raises(ValueError):
        current.replace_pin(refs, {'path': 'old'}, {'path': 'new'})


def test_pin_replacement_preserves_order_and_other_bindings():
    refs = [{'path': 'runtime'}, {'path': 'old'}, {'path': 'calibration'}]
    replaced = current.replace_pin(refs, {'path': 'old'}, {'path': 'new'})
    assert replaced == [{'path': 'runtime'}, {'path': 'new'}, {'path': 'calibration'}]
    assert refs[1] == {'path': 'old'}


def test_preparation_refuses_existing_output_before_any_runtime_work(tmp_path, monkeypatch):
    path = tmp_path / 'existing.json'
    path.write_text('{}')
    monkeypatch.setattr(current, 'FILES', {'preparation': path})
    monkeypatch.setattr(current.previous, 'load', lambda: pytest.fail('Must not prepare again'))
    with pytest.raises(FileExistsError):
        current.prepare()


def test_preparation_changes_only_goal_and_policy_references(tmp_path, monkeypatch):
    old = b'Goal\n  >\nExisting full scope\n'
    goal = tmp_path / 'goal.txt'
    goal.write_bytes(b'Goal\n  >\n' + current.ADDITION.encode('ascii') + b'Existing full scope\n')
    monkeypatch.setattr(current, 'GOAL', goal)
    monkeypatch.setattr(current, 'ROOT', tmp_path)
    monkeypatch.setattr(current, 'FILES', {key: tmp_path / (key + '_new.json')
        for key in ('policy', 'readiness', 'preparation')})
    monkeypatch.setattr(current.previous, 'load', lambda: None)
    old_goal = dict(path=str(goal), bytes=len(old), sha256=hashlib.sha256(old).hexdigest())
    other = tmp_path / 'other.json'
    other.write_text('{}')
    other_ref = current.executor.record(other)
    refs = {}
    for key in ('lookup', 'binding', 'plan', 'recipe', 'resources'):
        refs[key] = current.previous.save(tmp_path / (key + '.json'), dict(unchanged=key))
    policy = dict(decision='reviewed', minimum_available_memory_bytes=137438953472,
        maximum_sample_period_s=35, evidence=[old_goal, other_ref])
    refs['policy'] = current.previous.save(tmp_path / 'policy.json', policy)
    ready = dict(decision='passed', isolation_required=False, environment_policy=refs['policy'],
        evidence=[old_goal, refs['policy'], other_ref])
    refs['readiness'] = current.previous.save(tmp_path / 'ready.json', ready)
    parent = dict(refs, source_commit='old', next_index=21, prefix=[other_ref])
    parent_ref = current.previous.save(tmp_path / 'parent.json', parent)
    monkeypatch.setattr(current.previous, 'PREPARATION', Path(parent_ref['path']))
    source = Path(current.__file__)
    def git(argv, **kwargs):
        if argv[1] == 'rev-parse': return 'new\n'
        if argv[-1].startswith('old:'): return old
        return source.read_bytes() if source.name in argv[-1] else goal.read_bytes()
    monkeypatch.setattr(current.subprocess, 'check_output', git)
    # The temporary root is intentionally different from the module location.
    original_relative = Path.relative_to
    monkeypatch.setattr(Path, 'relative_to', lambda self, root:
        Path(self.name) if self == source else original_relative(self, root))
    result = current.executor.read(current.prepare())
    assert {k: result[k] for k in refs if k not in ('policy', 'readiness')} == {
        k: refs[k] for k in refs if k not in ('policy', 'readiness')}
    policy_new = current.executor.read(result['policy'])
    ready_new = current.executor.read(result['readiness'])
    assert policy_new == dict(policy, evidence=[current.executor.record(goal), other_ref])
    assert ready_new == dict(ready, environment_policy=result['policy'],
        evidence=[current.executor.record(goal), result['policy'], other_ref])
    assert result['prefix'] == parent['prefix'] and result['next_index'] == 21
    assert current.executor.read(parent_ref) == parent
    assert current.executor.read(refs['policy']) == policy
    assert result['documentary_reaffirmation']['only_documentary_evidence_changed'] is True


def launch_fixture(tmp_path, monkeypatch, *, existing=False, bad_limits=False, bad_selection=False):
    work = tmp_path / 'work'
    work.mkdir()
    monkeypatch.setattr(current, 'WORK', work)
    prep = work / 'preparation.json'
    monkeypatch.setattr(current, 'FILES', {'preparation': prep})
    refs = {}
    for key in ('lookup', 'binding', 'plan', 'recipe', 'policy', 'readiness', 'resources'):
        path = work / (key + '.json')
        path.write_text(json.dumps(dict(decision='reviewed' if key == 'policy' else 'passed',
            review_reference='fixture', evidence=[])))
        refs[key] = current.executor.record(path)
    prep.write_text(json.dumps(refs))
    run = dict(index=24, method='orthohmm_satellite_v2', proteomes=12, repeat=2,
        measurement_directory=str(tmp_path / 'measurements/run_24/measurement'))
    plan = dict(runs=[{}] * 24 + [run] + [{}] * 2)
    monkeypatch.setattr(current.previous, 'load', lambda: ({}, {}, plan, refs['plan']))
    prefix = [{'path': 'resolved0'}, {'path': 'resolved17'}, {'path': 'resolved20'}]
    monkeypatch.setattr(current.previous, 'history', lambda *a: (prefix, {}))
    reviews = []
    monkeypatch.setattr(current.executor, 'review', lambda ref, expected: reviews.append((ref, expected)))
    monkeypatch.setattr(current.executor, 'select', lambda *a:
        ({} if bad_selection else run, {}, dict(progress=dict(index=24)), []))
    calls = []
    def command(label, argv):
        calls.append((label, argv))
        if label == 'submission_24': return '999\n'
        if label == 'held_controller_24':
            return (f'JobId=999 JobState=PENDING Reason=JobHeldUser NumCPUs={32 if bad_limits else 64} '
                f'NumTasks=1 CPUs/Task=64 MinMemoryNode=128G TimeLimit={current.executor.TIME_LIMIT} '
                f'Requeue=0 Command={current.previous.SCRIPT} WorkDir={current.ROOT}')
        if label == 'pinned_controller_24':
            return 'Comment=' + current.executor.record(work / 'request_run_24.json')['sha256']
        return ''
    monkeypatch.setattr(current.previous.first, 'command', command)
    if existing: (work / 'request_run_24.json').write_text('{}')
    return work, refs, prefix, reviews, calls


def test_launch_preserves_history_limits_and_live_gate(tmp_path, monkeypatch):
    work, refs, prefix, reviews, calls = launch_fixture(tmp_path, monkeypatch)
    result = current.launch(24)
    request = json.loads((work / 'request_run_24.json').read_text())
    assert request['readiness_review'] == refs['readiness']
    assert request['history'] == prefix and request['index'] == 24
    assert result['preparation'] == current.executor.record(work / 'preparation.json')
    assert result['native_phase'] == 'requires_live_observation'
    assert result['automatic_retry'] is False and result['uncontended_timing'] is False
    assert len(reviews) == 2
    assert [label for label, _ in calls] == ['submission_24', 'held_controller_24',
        'request_comment_24', 'pinned_controller_24', 'release_24']
    assert '--hold' in calls[0][1]


@pytest.mark.parametrize('mode', ['existing', 'bad_limits', 'bad_selection'])
def test_invalid_attempt_or_placement_never_releases(tmp_path, monkeypatch, mode):
    _, _, _, _, calls = launch_fixture(tmp_path, monkeypatch, **{mode: True})
    with pytest.raises(FileExistsError if mode == 'existing' else ValueError):
        current.launch(24)
    assert not any(label.startswith('release') for label, _ in calls)
    if mode == 'existing': assert not calls


def test_earlier_science_and_review_sources_are_not_modified():
    assert current.previous.__name__.endswith('continue_shared_prepared_panel_20261004')
    assert current.FILES['preparation'] != current.previous.PREPARATION
    source = Path(current.__file__).read_text()
    assert 'previous.prepare()' not in source
