import json
import os
from pathlib import Path
from types import SimpleNamespace

import psutil
import pytest

from benchmark_tools import threadripper_environment_worker as worker


def put(path, value):
    path.write_text(json.dumps(value))
    return worker.record(path)


@pytest.fixture
def setup(tmp_path, monkeypatch):
    boot = Path('/proc/sys/kernel/random/boot_id').read_text().strip()
    support = put(tmp_path / 'support.json', {'synthetic_test_only': True})
    process_policy = dict(schema='threadripper_process_policy_v2', boot_id=boot,
        review_reference='synthetic test only', ordinary_processes=[])
    process_ref = put(tmp_path / 'process.json', process_policy)
    policy = dict(schema='threadripper_environment_policy_v1', decision='reviewed', host='bizon',
        plan_sha256=worker.PLAN_SHA, evidence=[support], review_reference='synthetic test only',
        process_policy=process_ref, configuration_files=[support], maximum_foreign_average_cores=.1,
        maximum_pressure_percent=dict(cpu=0., io=0., memory=0.))
    policy_ref = put(tmp_path / 'policy.json', policy)
    readiness = put(tmp_path / 'readiness.json', dict(environment_policy=policy_ref))
    request = dict(job_id=42, index=0, recipe=support, readiness_review=readiness,
                   environment_preflight_path=str(tmp_path / 'environment_preflight.json'))
    request_ref = put(tmp_path / 'request.json', request)
    directory = tmp_path / 'measurement'
    directory.mkdir()
    marker = dict(job_id=42, index=0, request=request_ref,
        review_path=request['environment_preflight_path'], wait_seconds=20, native_released=False,
        requested_unix_ns=10**12)
    put(directory / 'environment_review_requested.json', marker)
    monkeypatch.setenv('SLURM_JOB_ID', '42')
    monkeypatch.setattr(worker.os, 'uname', lambda: SimpleNamespace(nodename='bizon'))
    monkeypatch.setattr(worker, 'job_scope', lambda pid, job: '/job')
    monkeypatch.setattr(worker, 'select', lambda *args: (dict(measurement_directory=str(directory)), None, None, [support]))
    def sample(t):
        pid = os.getpid()
        return dict(schema='threadripper_typed_process_snapshot_v1', boot_id=boot,
            started_monotonic_s=t, finished_monotonic_s=t+.5, errors=[], processes=[
                dict(pid=pid, created=1., cgroup='/job/step', name='observer',
                     observed_monotonic_s=t+.1, user_s=0., system_s=0.,
                     kernel_identity=dict(pid=pid, tgid=pid, kthread=0,
                         started_monotonic_s=t+.2, finished_monotonic_s=t+.3))])
    samples = [sample(10.), sample(12.)]
    def host(t):
        psi='some avg10=0.00 avg60=0.00 avg300=0.00 total=0\nfull avg10=0.00 avg60=0.00 avg300=0.00 total=0\n'
        return dict(started_monotonic_ns=t, finished_monotonic_ns=t+1, errors=[],
            raw=dict(boot_id=boot, online_cpus='0-191', cgroup_membership='0::/job/step'),
            optional={f'host_{resource}_pressure':psi for resource in ('cpu','memory','io')})
    hosts=[host(10**10),host(13*10**9)]
    def execute(**extra):
        sample_iter, host_iter = iter(samples), iter(hosts)
        options=dict(root=tmp_path, sample=lambda:next(sample_iter), host_sample=lambda:next(host_iter),
            image_check=lambda *args:[], collector_check=lambda *args:[support], clock=lambda:10**12+10**9,
            sleep=lambda seconds:None, waiter=lambda *args,**kwargs:None)
        return worker.respond(request_ref, policy_ref, **{**options,**extra})
    return SimpleNamespace(root=tmp_path, run=execute, samples=samples, hosts=hosts, marker=marker,
        directory=directory, request=request, policy=policy, process_policy=process_policy,
        request_ref=request_ref, policy_ref=policy_ref, support=support)


def test_success_is_only_a_preflight_and_atomic_files(setup):
    result=setup.run()
    assert result['decision']=='passed' and result['whole_run_observer_ready']
    assert not result['scientific_timings_admitted']
    assert not (setup.directory/'go.json').exists()
    for ref in result['evidence']: worker.check(ref)
    assert not list(setup.root.glob('*.pending'))
    with pytest.raises(FileExistsError): setup.run()


@pytest.mark.parametrize('change',['unknown','sampling_error','missing_type','pressure','host_error',
    'image_error','collector_error','source_drift','expired','interrupt'])
def test_failed_reviews_do_not_release_or_discard_evidence(setup,change):
    extra={}
    if change=='unknown':
        for sample in setup.samples:
            row=dict(sample['processes'][0],pid=10,cgroup='/other',name='python',
                     kernel_identity={**sample['processes'][0]['kernel_identity'],'pid':10,'tgid':10})
            sample['processes'].append(row)
    elif change=='sampling_error': setup.samples[1]['errors'].append(dict(pid=10,type='AccessDenied'))
    elif change=='missing_type': del setup.samples[1]['processes'][0]['kernel_identity']
    elif change=='pressure':
        setup.hosts[1]['optional']['host_io_pressure']=setup.hosts[1]['optional']['host_io_pressure'].replace('total=0','total=10000')
    elif change=='host_error': setup.hosts[1]['errors'].append(dict(field='host_io_pressure',type='OSError'))
    elif change in {'image_error','collector_error'}:
        def fail(*args): raise ValueError('synthetic check failure')
        extra['image_check' if change=='image_error' else 'collector_check']=fail
    elif change=='source_drift':
        def drift(*args):
            Path(setup.support['path']).write_text('changed')
            return []
        extra['image_check']=drift
    elif change=='expired':
        times=iter([10**12+10**9,10**12+21*10**9,10**12+22*10**9,10**12+23*10**9])
        extra['clock']=lambda:next(times)
    else:
        def interrupt(): raise KeyboardInterrupt()
        extra['sample']=interrupt
    if change=='interrupt':
        with pytest.raises(KeyboardInterrupt): setup.run(**extra)
        result=json.loads(Path(setup.request['environment_preflight_path']).read_text())
    else: result=setup.run(**extra)
    assert result['decision']=='failed'
    assert not (setup.directory/'go.json').exists()
    assert (setup.root/'environment_worker_evidence.json').is_file()


@pytest.mark.parametrize('key,value',[('job_id',43),('native_released',True),('requested_unix_ns',True),
                                      ('requested_unix_ns',1)])
def test_wrong_or_stale_marker_does_not_publish(setup,key,value):
    setup.marker[key]=value
    put(setup.directory/'environment_review_requested.json',setup.marker)
    with pytest.raises(ValueError): setup.run()
    assert not Path(setup.request['environment_preflight_path']).exists()


def test_policy_must_be_bound_by_readiness(setup):
    path=Path(setup.request['readiness_review']['path'])
    setup.request['readiness_review']=put(path,dict(environment_policy={'path':'/other'}))
    updated=put(Path(setup.request_ref['path']),setup.request)
    setup.request_ref.update(updated)
    put(setup.directory/'environment_review_requested.json',setup.marker)
    with pytest.raises(ValueError,match='Readiness does not bind'): setup.run()


def test_explicit_background_cpu_bound(setup,monkeypatch):
    original=worker.process_review
    def high_cpu(*args,**kwargs):
        result=original(*args,**kwargs)
        result['cpu_diagnostic']['sum_observed_foreign_average_cores']=1.
        return result
    monkeypatch.setattr(worker,'process_review',high_cpu)
    assert setup.run()['decision']=='failed'
    detail=json.loads((setup.root/'environment_worker_evidence.json').read_text())
    assert 'background CPU bound exceeded' in detail['error']


def test_deadline_crossed_after_evidence_serialization(setup):
    times=iter([10**12+10**9,10**12+2*10**9,10**12+3*10**9,10**12+21*10**9])
    result=setup.run(clock=lambda:next(times))
    assert result['decision']=='failed' and result['publication_deadline_expired']


def test_real_loaded_interpreter_identity():
    pid=os.getpid(); process=psutil.Process(pid)
    row=dict(pid=pid,created=process.create_time(),cgroup=worker.membership(pid),name=process.name())
    image=worker.record(Path(f'/proc/{pid}/exe').resolve())
    policy=dict(ordinary_processes=[dict(row,classification='ordinary_background',image=image)])
    result=worker.images(policy,dict(processes=[row]))
    assert len(result)==1 and result[0]['loaded_sha256']==image['sha256']


def test_wrong_reviewed_image_rejected(tmp_path):
    pid=os.getpid(); process=psutil.Process(pid)
    row=dict(pid=pid,created=process.create_time(),cgroup=worker.membership(pid),name=process.name())
    image=put(tmp_path/'not-the-interpreter.json',dict(synthetic_test_only=True))
    with pytest.raises(ValueError,match='differs from reviewed image'):
        worker.images(dict(ordinary_processes=[dict(row,classification='ordinary_background',image=image)]),
                      dict(processes=[row]))


def test_existing_guard_consumes_worker_response_before_budget(setup):
    from benchmark_tools.run_threadripper_scaling import EnvironmentalReleaseGuard
    (setup.directory/'environment_review_requested.json').unlink()
    calls=[]
    times=iter([10**12,10**12+2*10**9,10**12+3*10**9])
    def wait(path, seconds):
        assert seconds==20
        assert setup.run()['decision']=='passed'
    guard=EnvironmentalReleaseGuard(setup.request,setup.request_ref,lambda directory:calls.append(directory),
        clock=lambda:next(times),waiter=wait)
    guard(setup.directory)
    assert calls==[setup.directory]
    assert json.loads((setup.directory/'environment_release.json').read_text())['status']=='environment_review_bound'
    assert not (setup.directory/'go.json').exists()


@pytest.mark.parametrize('problem',[None,'released','missing_native','errors','wrong_index','extra_line'])
def test_collector_readiness_identity_and_parked_state(tmp_path,monkeypatch,problem):
    scope='/root/job_42'
    put(tmp_path/'ready.json',dict(pid=50,cgroup=f'0::{scope}/step_0/user/task_0\n'))
    row=dict(index=0,interval=None,observer_pid=60,snapshot=dict(errors=[],processes=[dict(pid=50),dict(pid=60)]))
    if problem=='released': put(tmp_path/'go.json',dict(go=True))
    elif problem=='missing_native': row['snapshot']['processes'].pop(0)
    elif problem=='errors': row['snapshot']['errors'].append(dict(pid=70,type='AccessDenied'))
    elif problem=='wrong_index': row['index']=1
    stream=tmp_path/'host_processes.jsonl'
    stream.write_text(json.dumps(row)+'\n'+(json.dumps(row)+'\n' if problem=='extra_line' else ''))
    monkeypatch.setattr(worker,'job_scope',lambda *args:scope)
    monkeypatch.setattr(worker,'membership',lambda pid:scope+'/step_0/user/task_0')
    checked=[]
    monkeypatch.setattr(worker,'live_identity',lambda r:checked.append(r['pid']))
    if problem:
        with pytest.raises(ValueError): worker.collector_ready(tmp_path,scope,42)
    else:
        refs=worker.collector_ready(tmp_path,scope,42)
        assert checked==[60,50]
        for ref in refs: worker.check(ref)
