"""Test and retain the pre-native failure resolution; never submit or retry jobs."""

import csv
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from benchmark_tools import run_threadripper_scaling as executor
from benchmark_tools.results import continue_shared_threadripper_panel_20261003 as previous
from benchmark_tools.verify_threadripper_controller import validate

MODULES = ('test_threadripper_environment_worker', 'test_threadripper_panel_progress',
    'test_bind_threadripper_panel_history', 'test_review_shared_prenative_failure',
    'test_manage_threadripper_environment_worker', 'test_run_threadripper_scaling',
    'test_periodic_host_observer', 'test_measure_threadripper_scaling',
    'test_review_threadripper_process_stream', 'test_verify_threadripper_controller')


def resolve():
    directory = previous.WORK / 'prenative_resolution_run_17_20261004'
    resolved_path = previous.WORK / 'resolved_session_17.json'
    if directory.exists() or resolved_path.exists():
        raise FileExistsError('Retain the existing resolution; never rerun or overwrite it')
    plan_path = ROOT / 'benchmark_tools/results/threadripper_private_commands_20260928.json'
    plan_ref = executor.record(plan_path)
    plan = executor.read_frozen(plan_path, executor.PRIVATE_PLAN_SHA)
    failure_ref = executor.record(ROOT / 'benchmark_tools/results/threadripper_shared_prenative_failure_22413.json')
    if failure_ref['sha256'] != '9d643efb649b5a7c67cdda51fc6e359b5baaa89d938fce5e45c276d3186ea115':
        raise ValueError('Original independently retained failure changed')
    failure = executor.read(failure_ref)
    executor.expect(failure, dict(index=17, job_id=22413, native_outcome='not_started', resources=None,
        status='pre_native_infrastructure_failure_reviewed', next_submission_authorized=False))
    for pin in [failure['source'], *failure['evidence']]:
        executor.check(pin)
    controller_ref = executor.record(ROOT / 'benchmark_tools/results/threadripper_terminal_controller_22413.json')
    controller = executor.read(controller_ref)
    allocation = validate(controller['stdout'], 22413, 'terminal', command=str(previous.SCRIPT),
        cwd=str(ROOT), time_limit=executor.TIME_LIMIT, allocation_mode='shared')
    argv = ['sacct', '-j', '22413', '--noheader', '--parsable2', '--format=JobID,State,ExitCode,Start,End']
    started = time.time_ns()
    query = subprocess.run(argv, text=True, capture_output=True, timeout=5, check=True)
    fields = allocation['fields']
    rows = [row for row in csv.reader(io.StringIO(query.stdout), delimiter='|') if row and row[0] == '22413']
    if rows != [['22413', 'FAILED', '1:0', fields['StartTime'], fields['EndTime']]]:
        raise ValueError('Fresh accounting does not corroborate the retained abort')
    prefix, _ = previous.history(17, plan_ref)
    directory.mkdir()
    accounting_ref = previous.save(directory / 'accounting.json', dict(command=argv,
        returncode=query.returncode, stdout=query.stdout, stderr=query.stderr,
        started_unix_ns=started, finished_unix_ns=time.time_ns(), retained_controller=controller_ref))
    tests = [ROOT / f'tests/unit/{name}.py' for name in MODULES]
    sources = [ROOT / 'benchmark_tools' / name for name in (
        'threadripper_environment_worker.py', 'threadripper_panel_progress.py',
        'bind_threadripper_panel_history.py', 'manage_threadripper_environment_worker.py',
        'run_threadripper_scaling.py', 'periodic_host_observer.py', 'measure_threadripper_scaling.py',
        'review_threadripper_process_stream.py', 'verify_threadripper_controller.py')]
    tested_refs = [executor.record(path) for path in [*sources, *tests, Path(__file__)]]
    log_path, junit = directory / 'regression.log', directory / 'regression.xml'
    command = [sys.executable, '-B', '-m', 'pytest', '-q', '--tb=short',
        '--junitxml=' + str(junit), *map(str, tests)]
    started = time.time_ns()
    env = os.environ.copy()
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    with log_path.open('x') as handle:
        result = subprocess.run(command, cwd=ROOT, env=env, stdin=subprocess.DEVNULL,
            stdout=handle, stderr=subprocess.STDOUT, timeout=180)
    finished = time.time_ns()
    for pin in tested_refs:
        executor.check(pin)
    suites = list(ET.parse(junit).getroot().iter('testsuite'))
    counts = {name: sum(int(suite.attrib.get(name, '0')) for suite in suites)
        for name in ('tests', 'failures', 'errors', 'skipped')}
    passed = result.returncode == 0 and counts['tests'] > 0 and counts['failures'] == counts['errors'] == 0
    validation_ref = previous.save(directory / 'repair_validation.json', dict(
        status='immutable_initial_observation_regression_passed' if passed else 'repair_regression_failed',
        returncode=result.returncode, command=command, counts=counts,
        started_unix_ns=started, finished_unix_ns=finished,
        evidence=[*tested_refs, executor.record(log_path), executor.record(junit)],
        scientific_execution_authorized=False, scientific_timings_admitted=False,
        limitations=['Focused software regression on the current worktree, not production preflight or full runtime rebinding.',
            'A recorded repair never changes the original failed attempt or authorizes a native retry.']))
    if not passed:
        raise ValueError('Repair validation failed; preserve evidence and do not advance')
    decisions = dict(runtime='passed', environment='failed', resources='unresolved', outputs_or_failure='passed')
    reviews = {category: previous.save(directory / (category + '_review.json'), dict(
        schema='threadripper_panel_review_v1', index=17, job_id=22413,
        plan_sha256=plan_ref['sha256'], category=category, decision=decision,
        execution_scope=executor.SHARED_SCOPE, uncontended_timing=False,
        evidence=[failure_ref, accounting_ref], source=executor.record(__file__),
        review_reference='Retain the independently reviewed pre-native abort: failed preflight, absent native endpoints, recorded runtime checks and fresh terminal corroboration. No success verdict is fabricated.'))
        for category, decision in decisions.items()}
    original = dict(schema='threadripper_panel_session_v1', index=17, job_id=22413,
        plan_sha256=plan_ref['sha256'], phase='terminal', controller=controller_ref,
        native_outcome='not_started', pre_native_audit=failure_ref, reviews=reviews)
    original_ref = previous.save(directory / 'original_session.json', original)
    resolution_ref = previous.save(directory / 'resolution.json', dict(
        schema='threadripper_pre_native_failure_resolution_v1', index=17, job_id=22413,
        plan_sha256=plan_ref['sha256'], execution_scope=executor.SHARED_SCOPE,
        kind='pre_native_process_stream_identity_failure', decision='retain_excluded_attempt_and_advance',
        comparative_timing_eligible=False, automatic_retry=False, scientific_timings_admitted=False,
        original_session=original_ref, pre_native_audit=failure_ref,
        environment_review=reviews['environment'], repair_validation=validation_ref,
        evidence=[failure_ref, original_ref, *reviews.values(), accounting_ref, validation_ref, *tested_refs],
        source=executor.record(__file__),
        review_reference='Explicitly resolve the retained pre-native observer-stream race after immutable-first-sample regression. Advance without retry; no resource or environment pass is inferred.',
        limitations=['Original failed scheduler, failed preflight and missing native resources are unchanged.',
            'The aborted identity does not contribute a timing repeat; no retry or fastest-repeat selection.',
            'Current-source runtime/protocol bindings and a fresh native handoff remain required before index 18.']))
    resolved_ref = previous.save(resolved_path, dict(original, resolution=resolution_ref))
    bound = executor.bind(plan_ref, [*prefix, resolved_ref], command=str(previous.SCRIPT),
        cwd=str(ROOT), time_limit=executor.TIME_LIMIT, allocation_mode='shared')
    if bound['progress']['index'] != 18 or bound['progress']['status'] != 'next_identity_requires_preflight':
        raise ValueError('Resolved real history does not lead to the next frozen identity')
    result = dict(status='pre_native_abort_resolved_without_retry', index=17, job_id=22413,
        next_index=18, next_identity={key:plan['runs'][18][key] for key in ('method', 'proteomes', 'repeat')},
        resolved_session=resolved_ref, resolution=resolution_ref, repair_validation=validation_ref,
        progress=bound['progress'], evidence=bound['evidence'], source=executor.record(__file__),
        next_submission_authorized=False, native_runs_started=False, scientific_timings_admitted=False)
    previous.save(directory / 'result.json', result)
    return result


if __name__ == '__main__':
    result = resolve()
    print(json.dumps({key:result[key] for key in ('status', 'next_index', 'resolved_session', 'repair_validation')}, indent=2))
