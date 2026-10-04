import json
from pathlib import Path

import pytest

from benchmark_tools.bind_threadripper_panel_history import bind, REVIEWS
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.prepare_scaling_inputs import planned_runs
from tests.unit.test_verify_threadripper_controller import RAW


def put(path, data):
    path.write_text(json.dumps(data))
    return record(path)


def setup(tmp_path):
    plan = put(tmp_path / "plan.json", dict(runs=planned_runs()))
    controller = put(tmp_path / "controller.json", dict(
        command=["scontrol", "show", "job", "42", "--oneliner"], returncode=0,
        stdout=RAW.replace("RUNNING", "COMPLETED")))
    support = put(tmp_path / "support.json", dict(diagnostic=True))
    reviews = {k: put(tmp_path / (k+".json"), dict(schema="threadripper_panel_review_v1",
        index=0, job_id=42, plan_sha256=plan["sha256"], category=k,
        decision="passed", evidence=[support])) for k in REVIEWS}
    session = dict(schema="threadripper_panel_session_v1", index=0, job_id=42,
        plan_sha256=plan["sha256"], controller=controller, reviews=reviews,
        phase="terminal", native_outcome="exited_zero")
    session["native_audit"] = put(tmp_path / "native_audit.json", dict(job_id=42,
        native_outcome="exited_zero", status="native_success_outputs_verified"))
    return plan, session


def run(tmp_path, plan, session):
    ref = put(tmp_path / "session.json", session)
    return bind(plan, [ref], command="/recipe/run.sh", cwd="/recipe")


def test_bound_history_advances_but_does_not_authorize(tmp_path):
    plan, session = setup(tmp_path)
    result = run(tmp_path, plan, session)
    assert result["progress"]["index"] == 1
    assert result["scientific_execution_authorized"] is False


def test_missing_review_waits(tmp_path):
    plan, session = setup(tmp_path)
    session["reviews"] = None
    assert run(tmp_path, plan, session)["progress"]["status"] == "terminal_attempt_requires_review"


@pytest.mark.parametrize("actual,expected,accepted", [
    ("1-02:00:00", "1-02:00:00", True),
    ("1-00:00:00", "1-02:00:00", False),
    ("1-02:00:00", "1-00:00:00", False)])
def test_frozen_scheduler_envelope_is_propagated(tmp_path, actual, expected, accepted):
    plan, session = setup(tmp_path)
    path = tmp_path / "controller.json"
    controller = json.loads(path.read_text())
    controller["stdout"] = controller["stdout"].replace("1-00:00:00", actual)
    session["controller"] = put(path, controller)
    ref = put(tmp_path / "session.json", session)
    def invoke():
        return bind(plan, [ref], command="/recipe/run.sh", cwd="/recipe", time_limit=expected)
    if accepted:
        result = invoke()
        assert result["progress"]["index"] == 1
        assert result["controller_policy"]["time_limit"] == expected
        assert not result["scientific_execution_authorized"]
    else:
        with pytest.raises(ValueError, match="allocation"):
            invoke()


@pytest.mark.parametrize("field,value", [("job_id", 43), ("index", 1),
    ("category", "wrong"), ("plan_sha256", "wrong"), ("evidence", [])])
def test_rehashed_wrong_review_rejected(tmp_path, field, value):
    plan, session = setup(tmp_path)
    path = tmp_path / "runtime.json"
    data = json.loads(path.read_text())
    data[field] = value
    session["reviews"]["runtime"] = put(path, data)
    with pytest.raises(ValueError):
        run(tmp_path, plan, session)


@pytest.mark.parametrize("name", ["controller.json", "support.json", "runtime.json", "plan.json"])
def test_changed_evidence_rejected(tmp_path, name):
    plan, session = setup(tmp_path)
    (tmp_path / name).write_text("{}")
    with pytest.raises(ValueError):
        run(tmp_path, plan, session)


@pytest.mark.parametrize("field,value", [("returncode", 1), ("returncode", False),
    ("command", ["scontrol", "show", "job", "43", "--oneliner"]),
    ("stdout", RAW.replace("ExitCode=0:0", "ExitCode=1:0").replace("RUNNING", "COMPLETED"))])
def test_controller_contradiction_rejected(tmp_path, field, value):
    plan, session = setup(tmp_path)
    path = tmp_path / "controller.json"
    data = json.loads(path.read_text())
    data[field] = value
    session["controller"] = put(path, data)
    with pytest.raises(ValueError):
        run(tmp_path, plan, session)


@pytest.mark.parametrize("field,value", [("job_id", 43), ("native_outcome", "timed_out"),
                                       ("status", "unverified")])
def test_wrong_native_audit_rejected(tmp_path, field, value):
    plan, session = setup(tmp_path)
    path = tmp_path / "native_audit.json"
    data = json.loads(path.read_text())
    data[field] = value
    session["native_audit"] = put(path, data)
    with pytest.raises(ValueError, match="Native audit"):
        run(tmp_path, plan, session)


def resolved_setup(tmp_path):
    plan, session = setup(tmp_path)
    controller_path = tmp_path / 'controller.json'
    controller = json.loads(controller_path.read_text())
    controller['stdout'] = controller['stdout'].replace('COMPLETED', 'FAILED').replace('ExitCode=0:0', 'ExitCode=1:0').replace('NumCPUs=192', 'NumCPUs=64')
    session['controller'] = put(controller_path, controller)
    replay_ref = put(tmp_path / 'environment_replay.json', dict(processes=dict(
        failures={'sample_period_bound_exceeded': 1}, sampled_process_policy_satisfied=False,
        execution_scope='shared_host_matched_resources', job_scope='/slurm/job_42',
        command_bracketed_by_samples=True, records=3, intervals=2, policy_matched_intervals=2),
        pressure=dict(sampled_pressure_evidence_satisfied=True, failures={})))
    review_path = tmp_path / 'environment.json'
    environment = json.loads(review_path.read_text())
    environment.update(decision='failed', evidence=[replay_ref])
    session['reviews']['environment'] = put(review_path, environment)
    result_ref = put(tmp_path / 'result.json', dict(status='executor_failed', job_id=42, index=0,
        error='Whole-run sampled environment policy was not satisfied',
        wrapper=dict(status='command_exited_zero', measurement=dict(native=dict(exit_code=0, timed_out=False)))))
    original = put(tmp_path / 'original_session.json', session)
    resolution = dict(schema='threadripper_monitoring_failure_resolution_v1', index=0, job_id=42,
        plan_sha256=plan['sha256'], execution_scope='shared_host_matched_resources',
        kind='post_native_process_cadence_failure', decision='retain_excluded_attempt_and_advance',
        comparative_timing_eligible=False, automatic_retry=False, scientific_timings_admitted=False,
        original_session=original, native_audit=session['native_audit'],
        environment_review=session['reviews']['environment'], environment_replay=replay_ref,
        executor_result=result_ref, evidence=[record(tmp_path / 'support.json')],
        review_reference='Synthetic resolution, not real execution approval')
    session['resolution'] = put(tmp_path / 'resolution.json', resolution)
    return plan, session, resolution


def bind_resolved(tmp_path, plan, session, allocation_mode='shared'):
    return bind(plan, [put(tmp_path / 'session.json', session)], command='/recipe/run.sh',
        cwd='/recipe', allocation_mode=allocation_mode)


def test_bound_resolution_advances_only_shared_panel_without_rewriting_failure(tmp_path):
    plan, session, _ = resolved_setup(tmp_path)
    result = bind_resolved(tmp_path, plan, session)
    assert result['progress']['index'] == 1
    assert result['attempts'][0]['review']['environment'] == 'failed'
    assert result['attempts'][0]['scheduler_state'] == 'FAILED'
    assert result['progress']['reviewed_attempts'][0]['excluded_from_comparative_timing']
    with pytest.raises(ValueError):
        bind_resolved(tmp_path, plan, session, allocation_mode='exclusive')


@pytest.mark.parametrize('field,value', [('plan_sha256', 'wrong'), ('job_id', 43),
    ('index', False), ('kind', 'other'), ('comparative_timing_eligible', True),
    ('review_reference', ''), ('evidence', [])])
def test_rehashed_resolution_identity_or_admission_drift_rejected(tmp_path, field, value):
    plan, session, resolution = resolved_setup(tmp_path)
    resolution[field] = value
    session['resolution'] = put(tmp_path / 'resolution.json', resolution)
    with pytest.raises(ValueError):
        bind_resolved(tmp_path, plan, session)


@pytest.mark.parametrize('fault', ['identity', 'pressure', 'scope', 'bracketing', 'records', 'matched', 'counter', 'native', 'executor'])
def test_resolution_cannot_waive_missing_or_different_evidence(tmp_path, fault):
    plan, session, resolution = resolved_setup(tmp_path)
    replay = json.loads((tmp_path / 'environment_replay.json').read_text())
    if fault == 'identity': replay['processes']['failures'] = {'process_policy_mismatch': 1}
    elif fault == 'pressure': replay['pressure']['sampled_pressure_evidence_satisfied'] = False
    elif fault == 'scope': replay['processes']['job_scope'] = '/slurm/job_43'
    elif fault == 'bracketing': replay['processes']['command_bracketed_by_samples'] = False
    elif fault == 'records': replay['processes']['records'] += 1
    elif fault == 'matched': replay['processes']['policy_matched_intervals'] -= 1
    elif fault == 'counter': replay['processes']['failures']['sample_period_bound_exceeded'] = True
    replay_ref = put(tmp_path / 'environment_replay.json', replay)
    review = json.loads((tmp_path / 'environment.json').read_text())
    review['evidence'] = [replay_ref]
    session['reviews']['environment'] = put(tmp_path / 'environment.json', review)
    resolution.update(environment_replay=replay_ref, environment_review=session['reviews']['environment'])
    result = json.loads((tmp_path / 'result.json').read_text())
    if fault == 'native': result['wrapper']['measurement']['native']['exit_code'] = 1
    elif fault == 'executor': result['error'] = 'Other infrastructure failure'
    resolution['executor_result'] = put(tmp_path / 'result.json', result)
    resolution['original_session'] = put(tmp_path / 'original_session.json', {k: v for k, v in session.items() if k != 'resolution'})
    session['resolution'] = put(tmp_path / 'resolution.json', resolution)
    with pytest.raises(ValueError):
        bind_resolved(tmp_path, plan, session)


def pre_native_setup(tmp_path):
    plan, session = setup(tmp_path)
    controller = json.loads((tmp_path / 'controller.json').read_text())
    controller['stdout'] = controller['stdout'].replace('COMPLETED', 'FAILED').replace('ExitCode=0:0', 'ExitCode=1:0').replace('NumCPUs=192', 'NumCPUs=64')
    session['controller'] = put(tmp_path / 'controller.json', controller)
    support = record(tmp_path / 'support.json')
    pins = [put(tmp_path / 'go.json', {'abort': True}),
        put(tmp_path / 'aborted_before_native.json', {'status': 'observer_did_not_release_native'}),
        put(tmp_path / 'environment_release.json', {'status': 'environment_release_failed'}),
        put(tmp_path / 'environment_preflight.json', dict(decision='failed', index=0, job_id=42))]
    audit = dict(status='pre_native_infrastructure_failure_reviewed', index=0, job_id=42,
        native_outcome='not_started', execution_scope='shared_host_matched_resources',
        resources=None, comparative_timing_eligible=False, automatic_retry=False,
        scientific_timings_admitted=False, failure_kind='mutable_process_stream_identity_during_preflight',
        scheduler_state='FAILED', scheduler_exit_code='1:0', source=support, evidence=pins)
    session.pop('native_audit')
    session.update(native_outcome='not_started', pre_native_audit=put(tmp_path / 'pre_native_audit.json', audit))
    for category, decision in [('environment', 'failed'), ('resources', 'unresolved')]:
        review = json.loads((tmp_path / (category + '.json')).read_text())
        review['decision'] = decision
        session['reviews'][category] = put(tmp_path / (category + '.json'), review)
    repair = put(tmp_path / 'repair.json', dict(status='immutable_initial_observation_regression_passed',
        returncode=0, evidence=[support]))
    resolution = dict(schema='threadripper_pre_native_failure_resolution_v1', index=0, job_id=42,
        plan_sha256=plan['sha256'], execution_scope='shared_host_matched_resources',
        kind='pre_native_process_stream_identity_failure', decision='retain_excluded_attempt_and_advance',
        comparative_timing_eligible=False, automatic_retry=False, scientific_timings_admitted=False,
        original_session=put(tmp_path / 'original_session.json', session),
        pre_native_audit=session['pre_native_audit'], environment_review=session['reviews']['environment'],
        repair_validation=repair, evidence=[support], review_reference='Synthetic pre-native failure resolution')
    session['resolution'] = put(tmp_path / 'resolution.json', resolution)
    return plan, session, resolution, audit


def test_bound_pre_native_abort_keeps_failure_and_missing_resources(tmp_path):
    plan, session, _, _ = pre_native_setup(tmp_path)
    result = bind_resolved(tmp_path, plan, session)
    assert result['progress']['index'] == 1
    retained = result['progress']['reviewed_attempts'][0]
    assert retained['native_outcome'] == 'not_started'
    assert retained['scheduler_state'] == 'FAILED'
    assert retained['native_resources_available'] is False
    assert retained['original_review']['resources'] == 'unresolved'
    assert result['scientific_execution_authorized'] is False
    with pytest.raises(ValueError):
        bind_resolved(tmp_path, plan, session, allocation_mode='exclusive')


@pytest.mark.parametrize('fault', ['resources', 'kind', 'native_outcome', 'job', 'index', 'retry',
    'admission', 'source_drift', 'missing_evidence', 'success_gate', 'abort_marker',
    'passed_preflight', 'other_job_preflight', 'released', 'repair_failed',
    'repair_bool_exit', 'repair_empty', 'changed_original', 'borrowed_audit', 'no_resolution'])
def test_pre_native_resolution_cannot_disguise_or_borrow_other_evidence(tmp_path, fault):
    plan, session, resolution, audit = pre_native_setup(tmp_path)
    if fault == 'resources': audit['resources'] = {'wall_seconds': 0}
    elif fault == 'kind': audit['failure_kind'] = 'other_failure'
    elif fault == 'native_outcome': audit['native_outcome'] = 'exited_zero'
    elif fault == 'job': audit['job_id'] = 43
    elif fault == 'index': audit['index'] = False
    elif fault == 'retry': audit['automatic_retry'] = True
    elif fault == 'admission': audit['scientific_timings_admitted'] = True
    elif fault == 'source_drift': (tmp_path / 'support.json').write_text('{}')
    elif fault == 'missing_evidence': audit['evidence'] = []
    elif fault in ('success_gate', 'abort_marker', 'passed_preflight', 'other_job_preflight', 'released'):
        name = {'success_gate':'go.json', 'abort_marker':'aborted_before_native.json',
            'passed_preflight':'environment_preflight.json', 'other_job_preflight':'environment_preflight.json',
            'released':'environment_release.json'}[fault]
        data = json.loads((tmp_path / name).read_text())
        if fault == 'success_gate': data = {'go': True}
        elif fault == 'other_job_preflight': data['job_id'] = 43
        elif fault == 'passed_preflight': data['decision'] = 'passed'
        else: data['status'] = 'success'
        new_ref = put(tmp_path / name, data)
        audit['evidence'] = [new_ref if Path(pin['path']).name == name else pin for pin in audit['evidence']]
    elif fault.startswith('repair_'):
        repair = json.loads((tmp_path / 'repair.json').read_text())
        if fault == 'repair_failed': repair['returncode'] = 1
        elif fault == 'repair_bool_exit': repair['returncode'] = False
        else: repair['evidence'] = []
        resolution['repair_validation'] = put(tmp_path / 'repair.json', repair)
    session['pre_native_audit'] = put(tmp_path / 'pre_native_audit.json', audit)
    resolution['pre_native_audit'] = session['pre_native_audit']
    original = {k:v for k,v in session.items() if k != 'resolution'}
    if fault == 'changed_original': original['native_outcome'] = 'exited_zero'
    resolution['original_session'] = put(tmp_path / 'original_session.json', original)
    if fault == 'borrowed_audit': resolution['pre_native_audit'] = record(tmp_path / 'support.json')
    session['resolution'] = put(tmp_path / 'resolution.json', resolution)
    if fault == 'no_resolution':
        session.pop('resolution')
        assert bind_resolved(tmp_path, plan, session)['progress']['status'] == 'infrastructure_failure_requires_resolution'
    else:
        with pytest.raises((ValueError, KeyError)):
            bind_resolved(tmp_path, plan, session)
