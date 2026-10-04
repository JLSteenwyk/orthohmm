"""Revalidate index 17 and resolve index 20 without overwriting history or launching jobs."""

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
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.verify_threadripper_controller import validate

RESULTS = ROOT / 'benchmark_tools/results'
WORK = ROOT / 'benchmarks/work/threadripper_shared_execution_20261003'


def put(path, data):
    save(path, data)
    return executor.record(path)


def terminal(controller_ref, job):
    controller = executor.read(controller_ref)
    fields = validate(controller['stdout'], job, 'terminal',
        command=str(ROOT / 'benchmark_tools/run_threadripper_shared_scaling.sh'),
        cwd=str(ROOT), time_limit=executor.TIME_LIMIT, allocation_mode='shared')['fields']
    argv = ['sacct', '-j', str(job), '--noheader', '--parsable2', '--format=JobID,State,ExitCode,Start,End']
    started = time.time_ns()
    query = subprocess.run(argv, text=True, capture_output=True, timeout=5, check=True)
    rows = [row for row in csv.reader(io.StringIO(query.stdout), delimiter='|') if row and row[0] == str(job)]
    if rows != [[str(job), 'FAILED', '1:0', fields['StartTime'], fields['EndTime']]]:
        raise ValueError('Fresh accounting does not corroborate the retained failure')
    return dict(command=argv, returncode=query.returncode, stdout=query.stdout,
        stderr=query.stderr, started_unix_ns=started, finished_unix_ns=time.time_ns(),
        retained_controller=controller_ref)


def resolve():
    directory = WORK / 'preparation_sync_resolution_20261004'
    paths = {index: WORK / f'resolved_session_{index:02d}_preparation_sync_20261004.json' for index in (17, 20)}
    public = RESULTS / 'threadripper_preparation_sync_history_20261004.json'
    if directory.exists() or public.exists() or any(path.exists() for path in paths.values()):
        raise FileExistsError('Retain existing history; never retry or overwrite')
    plan_path = RESULTS / 'threadripper_private_commands_20260928.json'
    plan_ref = executor.record(plan_path)
    plan = executor.read_frozen(plan_path, executor.PRIVATE_PLAN_SHA)
    validated = executor.record(RESULTS / 'threadripper_preparation_sync_validation_20261004.json')
    immutable = executor.record(RESULTS / 'threadripper_immutable_sample_revalidation_20261004.json')
    for ref, status in ((validated, 'preparation_synchronization_regression_passed'),
                        (immutable, 'immutable_initial_observation_regression_passed')):
        proof = executor.read(ref)
        executor.expect(proof, dict(status=status, returncode=0, native_runs_started=False,
            scientific_execution_authorized=False, scientific_timings_admitted=False))
        if proof['missing_required_tests'] or any(proof['counts'][key] for key in ('failures', 'errors', 'skipped')):
            raise ValueError('Captured repair regression is incomplete')
        for pin in proof['evidence']:
            executor.check(pin)
    old_ref = executor.record(RESULTS / 'threadripper_prenative_resolution_22413.json')
    if old_ref['sha256'] != 'd35501e6a969043c494a4bfb20f8a993972d5aea14f899c1c764d553f5fdf399':
        raise ValueError('Original index-17 resolution changed')
    old = executor.read(old_ref)
    original17 = executor.read(old['original_session'])
    failure_refs = {17: old['pre_native_audit'],
        20: executor.record(RESULTS / 'threadripper_shared_prenative_failure_22416.json')}
    failures = {index: executor.read(ref) for index, ref in failure_refs.items()}
    for index, job in ((17, 22413), (20, 22416)):
        executor.expect(failures[index], dict(index=index, job_id=job,
            status='pre_native_infrastructure_failure_reviewed', native_outcome='not_started',
            resources=None, comparative_timing_eligible=False, automatic_retry=False,
            next_submission_authorized=False, scientific_timings_admitted=False))
        for pin in [failures[index]['source'], *failures[index]['evidence']]:
            executor.check(pin)
    controllers = {17: original17['controller'],
        20: executor.record(RESULTS / 'threadripper_terminal_controller_22416.json')}
    accounting = {index: terminal(ref, failures[index]['job_id']) for index, ref in controllers.items()}
    directory.mkdir()
    originals, resolutions, sessions = {}, {}, {}
    for index in (17, 20):
        job, failure_ref = failures[index]['job_id'], failure_refs[index]
        accounting_ref = put(directory / f'accounting_{index}.json', accounting[index])
        if index == 17:
            original, original_ref = original17, old['original_session']
        else:
            decisions = dict(runtime='passed', environment='failed', resources='unresolved', outputs_or_failure='passed')
            reviews = {category: put(directory / f'{category}_20_review.json', dict(
                schema='threadripper_panel_review_v1', index=index, job_id=job,
                plan_sha256=plan_ref['sha256'], category=category, decision=decision,
                execution_scope=executor.SHARED_SCOPE, uncontended_timing=False,
                evidence=[failure_ref, accounting_ref], source=executor.record(__file__),
                review_reference='Retain recorded runtime checks, failed preflight and absent native endpoints; no passing environment/resource verdict.'))
                for category, decision in decisions.items()}
            original = dict(schema='threadripper_panel_session_v1', index=index, job_id=job,
                plan_sha256=plan_ref['sha256'], phase='terminal', controller=controllers[index],
                native_outcome='not_started', pre_native_audit=failure_ref, reviews=reviews)
            original_ref = put(directory / 'original_session_20.json', original)
        proof = immutable if index == 17 else validated
        evidence = [failure_ref, original_ref, *original['reviews'].values(), accounting_ref, proof]
        if index == 17:
            evidence.append(old_ref)
        resolution = dict(schema='threadripper_pre_native_failure_resolution_v1', index=index,
            job_id=job, plan_sha256=plan_ref['sha256'], execution_scope=executor.SHARED_SCOPE,
            kind='pre_native_process_stream_identity_failure' if index == 17 else 'pre_native_environment_response_deadline',
            decision='retain_excluded_attempt_and_advance', comparative_timing_eligible=False,
            automatic_retry=False, scientific_timings_admitted=False, original_session=original_ref,
            pre_native_audit=failure_ref, environment_review=original['reviews']['environment'],
            repair_validation=proof, evidence=evidence, source=executor.record(__file__),
            review_reference='Joint current-source regression validates immutable first-sample and preparation synchronization; preserve failed outcomes and advance without retry.',
            limitations=['All original failed verdicts and absent native endpoints remain unchanged.',
                'No production handoff or current runtime/protocol binding is established by resolution.'])
        if index == 17:
            resolution['supersedes_for_future_execution'] = old_ref
        originals[index] = original_ref
        resolutions[index] = put(directory / f'resolution_{index}.json', resolution)
        sessions[index] = put(paths[index], dict(original, resolution=resolutions[index]))
    prefix = [executor.record(WORK / 'resolved_session_00.json')]
    for index in range(1, 21):
        if index in sessions:
            prefix.append(sessions[index])
        else:
            summary = executor.read(executor.record(WORK / f'review_run_{index:02d}/summary.json'))
            executor.expect(summary, dict(index=index, status='shared_attempt_independently_reviewed',
                shared_host_resources_reviewed=True, original_environment_protocol_passed=True))
            prefix.append(summary['session'])
    bound = executor.bind(plan_ref, prefix, command=str(ROOT / 'benchmark_tools/run_threadripper_shared_scaling.sh'),
        cwd=str(ROOT), time_limit=executor.TIME_LIMIT, allocation_mode='shared')
    if bound['progress']['status'] != 'next_identity_requires_preflight' or bound['progress']['index'] != 21:
        raise ValueError('Actual resolved prefix does not select index 21')
    report = dict(status='preparation_sync_failures_resolved_without_retry', resolved_sessions=sessions,
        resolutions=resolutions, original_sessions=originals, prefix=prefix, progress=bound['progress'],
        evidence=bound['evidence'], plan=plan_ref, next_index=21,
        next_identity={key: plan['runs'][21][key] for key in ('method', 'proteomes', 'repeat')},
        source=executor.record(__file__), next_submission_authorized=False, native_runs_started=False,
        scientific_timings_admitted=False, limitations=['Refreshed runtime/source/protocol/readiness and a fresh capacity/handoff remain required.',
            'Historical resolution and reporting bytes remain retained; neither aborted identity is retried.'])
    put(directory / 'result.json', report)
    put(public, report)
    print(json.dumps(dict(status=report['status'], next_identity=report['next_identity'],
        resolved_sessions=sessions, result=executor.record(public)), indent=2))
    return report


if __name__ == '__main__':
    resolve()
