"""Retain the index-17 pre-native abort; never submit, retry or admit resources."""

import argparse
import csv
import io
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from benchmark_tools import run_threadripper_scaling as executor
from benchmark_tools.verify_threadripper_controller import validate


def assess(run, request, documents, rows, present):
    executor.expect(run, dict(index=17, method='orthofinder_3_1_5_full', proteomes=12, repeat=1))
    executor.expect(request, dict(index=17, execution_scope=executor.SHARED_SCOPE))
    job = request['job_id']
    measurement = Path(run['measurement_directory'])
    preflight, detail = documents['preflight'], documents['detail']
    executor.expect(preflight, dict(index=17, job_id=job, decision='failed',
        whole_run_observer_ready=False, execution_scope=executor.SHARED_SCOPE,
        foreign_cpu_used_for_eligibility=False, uncontended_timing=False))
    if detail.get('error') != 'Frozen input/source identity changed: ' + str(measurement / 'host_processes.jsonl'):
        raise ValueError('Failure is not the retained mutable process-stream identity failure')
    executor.expect(detail, dict(status='failed', error_type='ValueError', native_release_authorized=False))
    executor.expect(documents['release'], dict(status='environment_release_failed', error_type='ValueError',
        error='Execution evidence identity or decision differs'))
    executor.expect(documents['wrapper'], dict(status='verified_wrapper_failed', error_type='ValueError',
        error='Execution evidence identity or decision differs', scientific_results_admitted=False))
    executor.expect(documents['result'], dict(index=17, job_id=job, status='executor_failed',
        error='Environmental worker was not successfully joined at release', automatic_retry=False,
        next_submission_authorized=False, scientific_timings_admitted=False))
    executor.expect(documents['lifecycle'], dict(terminal=True, exit_code=1, status='failed_or_cancelled'))
    if documents['go'] != {'abort': True} or documents['aborted'] != {'status': 'observer_did_not_release_native'}:
        raise ValueError('Native gate did not retain the expected abort')
    if present:
        raise ValueError('Native artifacts exist; cannot classify this as a pre-native abort')
    if len(rows) != 2 or [row['index'] for row in rows] != [0, 1] or rows[0]['interval'] is not None:
        raise ValueError('Retained process stream differs from the two-sample failure')
    observer = rows[0]['observer_pid']
    if type(observer) is not int or observer != rows[1]['observer_pid']:
        raise ValueError('Process observer changed identity')
    samples = [row['snapshot'] for row in rows]
    snapshots = detail['snapshots']
    if len(snapshots) != 2 or len({s['boot_id'] for s in [*samples, *snapshots]}) != 1:
        raise ValueError('Preflight and collector observations do not share a boot')
    first, second = [s['started_monotonic_s'] for s in samples]
    if not first < snapshots[0]['started_monotonic_s'] < second < snapshots[1]['finished_monotonic_s']:
        raise ValueError('The periodic append did not overlap the retained preflight')
    memory = detail['available_memory_bytes']
    minimum = documents['policy']['minimum_available_memory_bytes']
    if len(memory) != 2 or any(type(value) is not int or value < minimum for value in memory):
        raise ValueError('Failure has a separate capacity problem')
    for side in ('before', 'after'):
        checked = documents['wrapper'][side]['runtime']
        executor.expect(checked, dict(status='runtime_and_lookup_checked', scientific_execution_authorized=False))
        for name in ('orthohmm', 'orthofinder'):
            executor.expect(checked['lookup'][name], dict(status='native_python_lookup_matches'))
    return dict(status='pre_native_infrastructure_failure_reviewed', index=17, job_id=job,
        method=run['method'], proteomes=12, repeat=1, native_outcome='not_started',
        resources=None, comparative_timing_eligible=False, automatic_retry=False,
        next_submission_authorized=False, scientific_timings_admitted=False,
        execution_scope=executor.SHARED_SCOPE, uncontended_timing=False,
        failure_kind='mutable_process_stream_identity_during_preflight',
        available_memory_bytes=memory, minimum_available_memory_bytes=minimum,
        preflight_foreign_average_cores=detail['process_review']['cpu_diagnostic']['sum_observed_foreign_average_cores'],
        retained_process_samples=len(rows), first_process_interval_seconds=second-first,
        limitations=['A pre-native abort is neither failed OrthoFinder inference nor a native resource measurement.',
            'Recorded runtime checks are retained, not repeated or promoted to a full successful-run review.',
            'The failed preflight and scheduler outcome remain failed; no native endpoints are imputed.',
            'A targeted monitoring repair and explicit history resolution are required before progression.',
            'Contention remains accepted with unknown, potentially method-dependent distortion.'])


def collect(terminal_controller=None):
    plan_path = ROOT / 'benchmark_tools/results/threadripper_private_commands_20260928.json'
    plan_ref = executor.record(plan_path)
    run = executor.read_frozen(plan_path, executor.PRIVATE_PLAN_SHA)['runs'][17]
    work = ROOT / 'benchmarks/work/threadripper_shared_execution_20261003'
    request_ref = executor.record(work / 'request_run_17.json')
    request = executor.read(request_ref)
    job = request['job_id']
    argv = ['scontrol', 'show', 'job', str(job), '--oneliner']
    accounting = None
    if terminal_controller is None:
        started = time.time_ns()
        query = subprocess.run(argv, text=True, capture_output=True, timeout=5, check=True)
        controller = dict(command=argv, returncode=query.returncode, stdout=query.stdout,
            stderr=query.stderr, started_unix_ns=started, finished_unix_ns=time.time_ns())
    else:
        controller_ref = executor.record(terminal_controller)
        controller = executor.read(controller_ref)
        executor.expect(controller, dict(command=argv, returncode=0))
        accounting_argv = ['sacct', '-j', str(job), '--noheader', '--parsable2',
            '--format=JobID,State,ExitCode,Start,End']
        started = time.time_ns()
        query = subprocess.run(accounting_argv, text=True, capture_output=True, timeout=5, check=True)
        accounting = dict(command=accounting_argv, returncode=query.returncode,
            stdout=query.stdout, stderr=query.stderr, started_unix_ns=started,
            finished_unix_ns=time.time_ns(), retained_controller=controller_ref)
    allocation = validate(controller['stdout'], job, 'terminal',
        command=str(ROOT / 'benchmark_tools/run_threadripper_shared_scaling.sh'),
        cwd=str(ROOT), time_limit=executor.TIME_LIMIT, allocation_mode='shared')
    if (allocation['scheduler_state'], allocation['scheduler_exit_code']) != ('FAILED', '1:0'):
        raise ValueError('Scheduler is not the expected terminal failure')
    if allocation['fields']['Comment'] != request_ref['sha256']:
        raise ValueError('Scheduler request digest differs')
    if accounting is not None:
        rows = [row for row in csv.reader(io.StringIO(accounting['stdout']), delimiter='|')
            if row and row[0] == str(job)]
        fields = allocation['fields']
        if rows != [[str(job), fields['JobState'], fields['ExitCode'], fields['StartTime'], fields['EndTime']]]:
            raise ValueError('Fresh accounting does not corroborate the retained controller')
    measurement = Path(run['measurement_directory'])
    session = measurement.parent.parent / 'sessions/run_17'
    paths = dict(preflight=session / 'environment_preflight.json', detail=session / 'environment_worker_evidence.json',
        lifecycle=session / 'environment_worker_lifecycle.json', result=session / 'result.json',
        wrapper=measurement.parent / 'verification.json', release=measurement / 'environment_release.json',
        go=measurement / 'go.json', aborted=measurement / 'aborted_before_native.json',
        policy=Path(executor.read(request['readiness_review'])['environment_policy']['path']))
    refs = {key: executor.record(path) for key, path in paths.items()}
    documents = {key: executor.read(ref) for key, ref in refs.items()}
    if documents['detail']['request'] != request_ref or documents['result']['request'] != request_ref:
        raise ValueError('Failure evidence belongs to another request')
    for ref in documents['preflight']['evidence']:
        executor.check(ref)
    stream_ref = executor.record(measurement / 'host_processes.jsonl')
    with Path(stream_ref['path']).open() as handle:
        rows = [json.loads(line) for line in handle]
    native_paths = [measurement / name for name in ('native.log', 'done.json', 'lineage_report.json')]
    native_paths += list(measurement.glob('point_*.json'))
    native_paths.append(Path(run['configuration']['output']))
    report = assess(run, request, documents, rows, [str(path) for path in native_paths if path.exists()])
    lookup_refs = [checked['report'] for side in ('before', 'after')
        for checked in documents['wrapper'][side]['runtime']['lookup'].values()]
    report.update(plan=plan_ref, controller=controller, terminal_accounting=accounting,
        scheduler_state='FAILED', scheduler_exit_code='1:0',
        evidence=[plan_ref, request_ref, *refs.values(), stream_ref, *lookup_refs],
        source=executor.record(__file__))
    for ref in report['evidence']:
        executor.check(ref)
    if any(path.exists() for path in native_paths):
        raise ValueError('Native artifact inventory changed during review')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--terminal-controller', type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = collect(args.terminal_controller)
    with args.output.open('x') as handle:
        json.dump(report, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')
    print(json.dumps({key: report[key] for key in ('status', 'job_id', 'native_outcome', 'resources')}, indent=2))
