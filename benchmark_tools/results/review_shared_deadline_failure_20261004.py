"""Review the retained index-20 release timeout without retrying or admitting resources."""

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


def native_artifacts(run):
    measurement = Path(run['measurement_directory'])
    paths = [measurement / name for name in ('native.log', 'done.json', 'lineage_report.json')]
    paths += list(measurement.glob('point_*.json'))
    paths.append(Path(run['configuration']['metrics']))
    present = [str(path) for path in paths if path.exists() or path.is_symlink()]
    output = Path(run['configuration']['output'])
    # The harness prepares an empty output directory before the native gate.
    if output.is_symlink() or output.exists() and (not output.is_dir() or any(output.iterdir())):
        present.append(str(output))
    return present


def assess(run, request, documents, rows, present):
    executor.expect(run, dict(index=20, method='orthohmm_satellite_v2', proteomes=4, repeat=2))
    executor.expect(request, dict(index=20, execution_scope=executor.SHARED_SCOPE))
    job = request['job_id']
    preflight, detail, marker = (documents[key] for key in ('preflight', 'detail', 'marker'))
    executor.expect(preflight, dict(index=20, job_id=job, decision='failed',
        whole_run_observer_ready=False, execution_scope=executor.SHARED_SCOPE,
        publication_deadline_expired=True, foreign_cpu_used_for_eligibility=False,
        uncontended_timing=False))
    executor.expect(marker, dict(index=20, job_id=job, wait_seconds=20, native_released=False))
    executor.expect(detail, dict(status='deadline_expired', error_type='FileNotFoundError',
        native_release_authorized=False))
    pid = documents['ready']['pid']
    if type(pid) is not int or detail['error'] != f"[Errno 2] No such file or directory: '/proc/{pid}/cgroup'":
        raise ValueError('Failure does not identify the retained parked worker')
    for key in ('release', 'wrapper'):
        executor.expect(documents[key], dict(error_type='TimeoutError',
            error='Environmental worker response deadline exceeded'))
    executor.expect(documents['release'], dict(status='environment_release_failed'))
    executor.expect(documents['wrapper'], dict(status='verified_wrapper_failed', scientific_results_admitted=False))
    executor.expect(documents['result'], dict(index=20, job_id=job, status='executor_failed',
        error='Environmental worker was not successfully joined at release', automatic_retry=False,
        next_submission_authorized=False, scientific_timings_admitted=False))
    executor.expect(documents['lifecycle'], dict(terminal=True, exit_code=1, status='failed_or_cancelled'))
    if documents['go'] != {'abort': True} or documents['aborted'] != {'status': 'observer_did_not_release_native'}:
        raise ValueError('Native gate did not retain the expected abort')
    if present:
        raise ValueError('Native artifacts exist; this is not a pre-native abort')
    requested, started, finished = (marker['requested_unix_ns'],
        preflight['observation_started_unix_ns'], preflight['observation_finished_unix_ns'])
    if (any(type(value) is not int for value in (requested, started, finished))
            or not 0 < requested <= started < requested + 20_000_000_000 <= finished):
        raise ValueError('Retained response does not cross the original deadline')
    if len(rows) != 2 or [row['index'] for row in rows] != [0, 1] or rows[0]['interval'] is not None:
        raise ValueError('Retained collector history differs')
    if rows[0]['observer_pid'] != rows[1]['observer_pid']:
        raise ValueError('Collector identity changed')
    snapshots = detail['snapshots']
    if len(snapshots) != 2 or len({s['boot_id'] for s in [*snapshots, *[r['snapshot'] for r in rows]]}) != 1:
        raise ValueError('Observations do not share one boot')
    if documents['initial'] != rows[0]:
        raise ValueError('Immutable initial sample differs from stream')
    memory, minimum = detail['available_memory_bytes'], documents['policy']['minimum_available_memory_bytes']
    if len(memory) != 2 or any(type(value) is not int or value < minimum for value in memory):
        raise ValueError('There is a separate unsafe-capacity failure')
    executor.expect(detail['process_review'], dict(process_policy_matched=True,
        execution_scope=executor.SHARED_SCOPE, foreign_cpu_used_for_eligibility=False))
    for side in ('before', 'after'):
        checked = documents['wrapper'][side]['runtime']
        executor.expect(checked, dict(status='runtime_and_lookup_checked', scientific_execution_authorized=False))
        for name in ('orthohmm', 'orthofinder'):
            executor.expect(checked['lookup'][name], dict(status='native_python_lookup_matches'))
    return dict(status='pre_native_infrastructure_failure_reviewed', index=20, job_id=job,
        method=run['method'], proteomes=4, repeat=2, native_outcome='not_started',
        resources=None, comparative_timing_eligible=False, automatic_retry=False,
        next_submission_authorized=False, scientific_timings_admitted=False,
        execution_scope=executor.SHARED_SCOPE, uncontended_timing=False,
        failure_kind='environment_response_deadline_with_parked_worker_disappearance',
        request_to_review_start_seconds=(started-requested)/1e9,
        review_to_response_seconds=(finished-started)/1e9,
        request_to_response_seconds=(finished-requested)/1e9,
        response_deadline_seconds=20, parked_worker_pid=pid,
        available_memory_bytes=memory, minimum_available_memory_bytes=minimum,
        preflight_foreign_average_cores=detail['process_review']['cpu_diagnostic']['sum_observed_foreign_average_cores'],
        retained_process_samples=len(rows),
        limitations=['Pre-native infrastructure failure, not failed OrthoHMM inference or a resource measurement.',
            'The retained timeout and later missing parked-worker cgroup do not isolate why preparation was slow.',
            'Runtime checks and safe capacity are retained; the environmental preflight remains failed.',
            'No retry, endpoint imputation or continuation authorization; explicit history resolution is still required.',
            'Shared-host contention remains accepted with unknown, potentially method-dependent distortion.'])


def collect(controller_path):
    plan_path = ROOT / 'benchmark_tools/results/threadripper_private_commands_20260928.json'
    plan_ref = executor.record(plan_path)
    run = executor.read_frozen(plan_path, executor.PRIVATE_PLAN_SHA)['runs'][20]
    work = ROOT / 'benchmarks/work/threadripper_shared_execution_20261003'
    request_ref = executor.record(work / 'request_run_20.json')
    request = executor.read(request_ref)
    job = request['job_id']
    controller_ref = executor.record(controller_path)
    controller = executor.read(controller_ref)
    executor.expect(controller, dict(command=['scontrol', 'show', 'job', str(job), '--oneliner'], returncode=0))
    allocation = validate(controller['stdout'], job, 'terminal',
        command=str(ROOT / 'benchmark_tools/run_threadripper_shared_scaling.sh'),
        cwd=str(ROOT), time_limit=executor.TIME_LIMIT, allocation_mode='shared')
    fields = allocation['fields']
    if (allocation['scheduler_state'], allocation['scheduler_exit_code']) != ('FAILED', '1:0'):
        raise ValueError('Scheduler is not the expected terminal failure')
    if fields['Comment'] != request_ref['sha256']:
        raise ValueError('Controller request digest differs')
    argv = ['sacct', '-j', str(job), '--noheader', '--parsable2', '--format=JobID,State,ExitCode,Start,End']
    started = time.time_ns()
    query = subprocess.run(argv, text=True, capture_output=True, timeout=5, check=True)
    accounting = dict(command=argv, returncode=query.returncode, stdout=query.stdout,
        stderr=query.stderr, started_unix_ns=started, finished_unix_ns=time.time_ns())
    accounts = list(csv.reader(io.StringIO(query.stdout), delimiter='|'))
    if [row for row in accounts if row and row[0] == str(job)] != [[str(job),
            fields['JobState'], fields['ExitCode'], fields['StartTime'], fields['EndTime']]]:
        raise ValueError('Fresh accounting does not corroborate retained controller')
    measurement = Path(run['measurement_directory'])
    session = measurement.parent.parent / 'sessions/run_20'
    paths = dict(preflight=session / 'environment_preflight.json', detail=session / 'environment_worker_evidence.json',
        lifecycle=session / 'environment_worker_lifecycle.json', result=session / 'result.json',
        wrapper=measurement.parent / 'verification.json', release=measurement / 'environment_release.json',
        go=measurement / 'go.json', aborted=measurement / 'aborted_before_native.json',
        marker=measurement / 'environment_review_requested.json', ready=measurement / 'ready.json',
        initial=measurement / 'preflight_initial_process_sample.json',
        policy=Path(executor.read(request['readiness_review'])['environment_policy']['path']))
    refs = {key: executor.record(path) for key, path in paths.items()}
    documents = {key: executor.read(ref) for key, ref in refs.items()}
    if any(documents[key]['request'] != request_ref for key in ('detail', 'result', 'marker')):
        raise ValueError('Evidence belongs to another request')
    if documents['detail']['marker'] != refs['marker']:
        raise ValueError('Worker reviewed another marker')
    for ref in documents['preflight']['evidence']:
        executor.check(ref)
    stream_ref = executor.record(measurement / 'host_processes.jsonl')
    with Path(stream_ref['path']).open() as handle:
        rows = [json.loads(line) for line in handle]
    report = assess(run, request, documents, rows, native_artifacts(run))
    lookup_refs = [checked['report'] for side in ('before', 'after')
        for checked in documents['wrapper'][side]['runtime']['lookup'].values()]
    report.update(plan=plan_ref, terminal_controller=controller_ref, terminal_accounting=accounting,
        scheduler_state='FAILED', scheduler_exit_code='1:0',
        evidence=[plan_ref, request_ref, controller_ref, *refs.values(), stream_ref, *lookup_refs],
        source=executor.record(__file__))
    for ref in report['evidence']:
        executor.check(ref)
    if native_artifacts(run):
        raise ValueError('Native artifact inventory changed during review')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--terminal-controller', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = collect(args.terminal_controller)
    with args.output.open('x') as handle:
        json.dump(report, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')
    print(json.dumps({key: report[key] for key in ('status', 'job_id', 'native_outcome', 'resources')}, indent=2))
