"""Resolve one excluded attempt and release the next evidence-bound identity."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from benchmark_tools.results import prepare_shared_threadripper_launch_20261003 as first
from benchmark_tools import run_threadripper_scaling as executor
from benchmark_tools.derive_threadripper_resources import SCOPES, protocol

WORK, RESULTS, SCRIPT = first.WORK, first.RESULTS, first.SCRIPT
LOOKUP = RESULTS / 'threadripper_private_lookup_shared_repaired_20261003.json'
LOOKUP_SHA = '7c22f7594294295ae4ce949b37e8f6ea13aee7ce0c1610e76285740116f65b02'
PREPARATION = WORK / 'repaired_preparation.json'
RESOLVED = WORK / 'resolved_session_00.json'
FILES = {key: RESULTS / ('threadripper_shared_' + name + '_repaired_20261003.json')
    for key, name in dict(recipe='source_recipe', policy='environment_policy', readiness='readiness').items()}
RESOURCE = RESULTS / 'threadripper_resource_endpoints_shared_repaired_20261003.json'
save = first.save


def load():
    lookup = executor.read_frozen(LOOKUP, LOOKUP_SHA)
    binding = executor.read(lookup['binding'])
    if os.path.abspath(sys.executable) != binding['controller_python']['path']:
        raise ValueError('Use the frozen private controller')
    executor.check(binding['controller_python'])
    plan_path = RESULTS / 'threadripper_private_commands_20260928.json'
    plan = executor.read_frozen(plan_path, executor.PRIVATE_PLAN_SHA)
    return lookup, binding, plan, executor.record(plan_path)


def history(index, plan_ref):
    if type(index) is not int or not 1 <= index <= 27:
        raise ValueError('Continuation requires a reviewed one-attempt prefix')
    refs = [executor.record(RESOLVED)]
    for number in range(1, index):
        summary = executor.read(executor.record(WORK / f'review_run_{number:02d}/summary.json'))
        executor.expect(summary, dict(status='shared_attempt_independently_reviewed', index=number,
            shared_host_resources_reviewed=True, original_environment_protocol_passed=True))
        refs.append(summary['session'])
    bound = executor.bind(plan_ref, refs, command=str(SCRIPT), cwd=str(ROOT),
        time_limit=executor.TIME_LIMIT, allocation_mode='shared')
    expected = 'all_attempts_reviewed' if index == 27 else 'next_identity_requires_preflight'
    if bound['progress']['status'] != expected or bound['progress']['index'] != (None if index == 27 else index):
        raise ValueError('Existing prefix is live, unresolved, complete or reordered')
    return refs, bound


def resolve():
    _, _, _, plan_ref = load()
    resolution_path = WORK / 'monitoring_resolution_00.json'
    if resolution_path.exists() or RESOLVED.exists():
        raise FileExistsError('Preserve the existing resolution')
    summary_ref = executor.record(RESULTS / 'threadripper_shared_attempt_22396.json')
    summary = executor.read(summary_ref)
    executor.expect(summary, dict(status='shared_attempt_reviewed_with_environment_failure',
        index=0, job_id=22396, scheduler_state='FAILED', scheduler_exit_code='1:0',
        original_environment_protocol_passed=False, primary_resources_replayed=True,
        shared_host_resources_reviewed=False))
    original_ref = summary['session']
    original = executor.read(original_ref)
    reviewed = Path(original_ref['path']).parent
    replay_ref = executor.record(reviewed / 'environment_replay.json')
    result_ref = executor.record(ROOT / 'benchmarks/results/threadripper_scaling_v1/sessions/run_00/result.json')
    test_python = ROOT / 'benchmarks/work/release_alert_refresh_20261001/venv/bin/python'
    modules = ['test_threadripper_panel_progress', 'test_bind_threadripper_panel_history',
        'test_run_threadripper_scaling', 'test_periodic_host_observer', 'test_measure_threadripper_scaling',
        'test_review_threadripper_process_stream', 'test_verify_threadripper_controller']
    test_paths = [ROOT / ('tests/unit/' + name + '.py') for name in modules]
    command = [str(test_python), '-B', '-m', 'pytest', '-q', '--tb=short', *map(str, test_paths)]
    log_path = WORK / 'resolution_regression_00.log'
    started = time.time_ns()
    with log_path.open('x') as stream:
        result = subprocess.run(command, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT,
            env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1'), timeout=180)
    test_ref = save(WORK / 'resolution_regression_00.json', dict(command=command,
        returncode=result.returncode, started_unix_ns=started, finished_unix_ns=time.time_ns(),
        log=executor.record(log_path), source_tests=[executor.record(path) for path in test_paths]))
    if result.returncode:
        raise ValueError('Resolution regression failed; retain evidence without native retry')
    repair = [executor.record(ROOT / 'benchmark_tools' / name)
        for name in ('periodic_host_observer.py', 'measure_threadripper_scaling.py')]
    for pin in repair:
        relative = str(Path(pin['path']).relative_to(ROOT))
        if Path(pin['path']).read_bytes() != subprocess.check_output(['git', 'show', '6add4497:' + relative], cwd=ROOT):
            raise ValueError('Observer repair differs from its tested committed bytes')
    evidence = [summary_ref, original_ref, original['controller'], original['native_audit'],
        *original['reviews'].values(), replay_ref, result_ref, test_ref, executor.record(log_path),
        *repair, executor.record(RESULTS / 'THREADRIPPER_SHARED_ATTEMPT_22396.md'),
        executor.record(RESULTS / 'PUBLICATION_SHARED_HOST_AMENDMENT_20261003.md')]
    for pin in evidence:
        executor.check(pin)
    resolution_ref = save(resolution_path, dict(schema='threadripper_monitoring_failure_resolution_v1',
        index=0, job_id=22396, plan_sha256=plan_ref['sha256'],
        execution_scope=executor.SHARED_SCOPE, kind='post_native_process_cadence_failure',
        decision='retain_excluded_attempt_and_advance', comparative_timing_eligible=False,
        automatic_retry=False, scientific_timings_admitted=False,
        original_session=original_ref, native_audit=original['native_audit'],
        environment_review=original['reviews']['environment'], environment_replay=replay_ref,
        executor_result=result_ref, evidence=evidence, repair_commit='6add4497',
        source=executor.record(__file__), review_reference='Native/output/accounting checks pass; preserve the independently reproduced original cadence failure, exclude its timing and advance only after the anchored-observer repair and regression. No retry.',
        limitations=['Failed scheduler and environment verdicts are unchanged.',
            'Run 0 does not contribute an eligible repeat to a comparative timing aggregate.',
            'Resolution binds an engineering explanation, not an isolation or timing admission pass.']))
    resolved_ref = save(RESOLVED, dict(original, resolution=resolution_ref))
    _, bound = history(1, plan_ref)
    return dict(status='excluded_monitoring_failure_resolved_without_retry', session=resolved_ref,
        resolution=resolution_ref, progress=bound['progress'])


def prepare():
    if PREPARATION.exists() or any(path.exists() for path in [*FILES.values(), RESOURCE]):
        raise FileExistsError('Repaired preparation already exists')
    lookup, binding, _, plan_ref = load()
    prefix, bound = history(1, plan_ref)
    old = executor.read(executor.record(WORK / 'preparation.json'))
    old_policy = executor.read(old['policy'])
    if old_policy['process_policy']['path'] is None:
        raise ValueError('Missing original process policy')
    executor.check(old_policy['process_policy'])
    if executor.read(old_policy['process_policy'])['boot_id'] != Path('/proc/sys/kernel/random/boot_id').read_text().strip():
        raise ValueError('Process policy belongs to a previous boot')
    manifest = executor.read(binding['runtime_manifests'][0])
    by_path = {row['path']: row for row in manifest['records']}
    sources = [executor.record(path) for path in sorted(set((ROOT / 'benchmark_tools').glob('*.py')) | {SCRIPT})]
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    for pin in sources:
        path = Path(pin['path'])
        if path.read_bytes() != subprocess.check_output(['git', 'show', commit + ':' + str(path.relative_to(ROOT))], cwd=ROOT):
            raise ValueError('Uncommitted executor source')
        if path.suffix == '.py' and (by_path[pin['path']]['bytes'], by_path[pin['path']]['sha256']) != (pin['bytes'], pin['sha256']):
            raise ValueError('Runtime helper binding differs')
    recipe_ref = save(FILES['recipe'], dict(schema='threadripper_executor_recipe_v1', root=str(ROOT),
        sources=sources, source_commit=commit, execution_scope=executor.SHARED_SCOPE,
        scientific_timings_admitted=False))
    parent_path = RESULTS / 'threadripper_resource_endpoints_20260929.json'
    parent = executor.read_frozen(parent_path, 'dba9602c365184f7e802425b9116b8021175b225d1f2266c44077da924bcfcdc')
    if parent['primary_scopes'] != SCOPES:
        raise ValueError('Resource endpoint definitions changed')
    resource_ref = save(RESOURCE, dict(parent, sources=[executor.record(Path(pin['path'])) for pin in parent['sources']],
        plan=plan_ref, lookup=executor.record(LOOKUP), parent_protocol=executor.record(parent_path),
        decision='same_endpoints_prospectively_bound_to_repaired_recipe', source=executor.record(__file__),
        execution_scope=executor.SHARED_SCOPE, excluded_prefix=prefix,
        remaining_prerequisites=['fresh shared-host capacity/handoff and actual per-run reviews'],
        limitations=['Source pins refreshed before index 1; index 0 retains its original failed admission.',
            'Endpoint scopes and failure retention are unchanged.']))
    protocol(RESOURCE, resource_ref['sha256'])
    calibration_ref = executor.record(RESULTS / 'threadripper_async_calibration_audit_22380.json')
    calibration = executor.read(calibration_ref)
    if calibration['status'] != 'calibration_checks_passed' or not all(value is True for value in calibration['checks'].values()):
        raise ValueError('Retained accounting checks do not pass')
    for pin in calibration['evidence']:
        executor.check(pin)
    evidence = [executor.record(ROOT / 'benchmark_tools/PUBLICATION_GOAL_20261003.txt'),
        executor.record(RESULTS / 'PUBLICATION_SHARED_HOST_AMENDMENT_20261003.md'),
        executor.record(RESULTS / 'THREADRIPPER_SHARED_CONTINUATION_20261003.md'),
        calibration_ref, executor.record(LOOKUP), lookup['binding'], resource_ref,
        recipe_ref, *prefix, executor.record(__file__)]
    policy_ref = save(FILES['policy'], dict(old_policy, evidence=evidence,
        review_reference='Repaired observer source and explicit excluded prefix; original shared-host bounds and diagnostic roles unchanged'))
    ready_ref = save(FILES['readiness'], dict(schema='threadripper_shared_readiness_review_v1',
        decision='passed', execution_scope=executor.SHARED_SCOPE, plan_sha256=plan_ref['sha256'],
        lookup_sha256=LOOKUP_SHA, recipe_sha256=recipe_ref['sha256'], environment_policy=policy_ref,
        observer_accounting_validated=True, environment_policy_frozen=True,
        contention_annotation_required=True, isolation_required=False,
        review_reference='Retained accounting calibration, repaired source/runtime checks and reviewed excluded prefix; actual native handoff remains live-gated',
        evidence=[*evidence, policy_ref], native_handoff_established=False,
        controlled_workload_verified=False, observer_causal_slowdown_established=False,
        scientific_timings_admitted=False))
    result = dict(status='repaired_shared_execution_prepared', source=executor.record(__file__),
        source_commit=commit, lookup=executor.record(LOOKUP), binding=lookup['binding'],
        plan=plan_ref, recipe=recipe_ref, policy=policy_ref, readiness=ready_ref,
        resources=resource_ref, prefix=prefix, next_index=bound['progress']['index'],
        native_runs_started=False, scientific_timings_admitted=False)
    save(PREPARATION, result)
    return result


def launch(index):
    _, _, plan, plan_ref = load()
    prefix, _ = history(index, plan_ref)
    if index == 27:
        raise ValueError('Panel is complete; do not submit')
    preparation = executor.read(executor.record(PREPARATION))
    for key in ('lookup', 'binding', 'plan', 'recipe', 'policy', 'readiness', 'resources'):
        executor.check(preparation[key])
    run = plan['runs'][index]
    request_path = WORK / f'request_run_{index:02d}.json'
    session = Path(run['measurement_directory']).parent.parent / 'sessions' / f'run_{index:02d}'
    if request_path.exists() or session.exists() or Path(run['measurement_directory']).parent.exists():
        raise FileExistsError('Existing attempt must be observed/reviewed, never retried')
    suffix = f'{index:02d}'
    raw = first.command('submission_' + suffix, ['sbatch', '--hold', '--parsable', '--partition=gpu',
        '--chdir=' + str(ROOT), '--output=' + str(WORK / f'slurm_run_{suffix}-%j.log'),
        '--error=' + str(WORK / f'slurm_run_{suffix}-%j.log'), str(SCRIPT),
        '--request', str(request_path), '--request-sha256', 'scheduler-comment'])
    match = re.fullmatch(r'([0-9]+)(?:;[^\s]+)?\n?', raw)
    if match is None:
        raise ValueError('Cannot identify held job; preserve submission and inspect scheduler')
    job = int(match[1])
    held = first.command('held_controller_' + suffix, ['scontrol', 'show', 'job', str(job), '--oneliner'])
    pairs = re.findall(r'(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)', held)
    fields = dict(pairs)
    required = dict(JobId=str(job), JobState='PENDING', Reason='JobHeldUser', NumCPUs='64',
        NumTasks='1', MinMemoryNode='128G', TimeLimit=executor.TIME_LIMIT, Requeue='0',
        Command=str(SCRIPT), WorkDir=str(ROOT))
    required['CPUs/Task'] = '64'
    if len(fields) != len(pairs) or any(fields.get(key) != value for key, value in required.items()):
        raise ValueError('Held resources/identity differ; do not release')
    request = dict(schema='threadripper_execution_request_v1', job_id=job,
        execution_authorized=True, deployment='private_v2_20260928', execution_scope=executor.SHARED_SCOPE,
        runtime_lookup=preparation['lookup'], plan_sha256=plan_ref['sha256'], lookup_sha256=LOOKUP_SHA,
        index=index, allocation_cwd=str(ROOT), scheduler_command=str(SCRIPT),
        recipe=preparation['recipe'], history=prefix, readiness_review=preparation['readiness'],
        environment_preflight_path=str(session / 'environment_preflight.json'))
    request_ref = save(request_path, request)
    chosen, _, selected_history, sources = executor.select(request, ROOT, job)
    if chosen != run or selected_history['progress']['index'] != index:
        raise ValueError('Actual held-job selection differs from frozen next identity')
    for pin in [request_ref, *sources]:
        executor.check(pin)
    first.command('request_comment_' + suffix, ['scontrol', 'update', 'JobId=' + str(job),
        'Comment=' + request_ref['sha256']])
    pinned = first.command('pinned_controller_' + suffix, ['scontrol', 'show', 'job', str(job), '--oneliner'])
    if ('Comment=' + request_ref['sha256']) not in pinned.split():
        raise ValueError('Request digest not bound in scheduler; leave held')
    first.command('release_' + suffix, ['scontrol', 'release', str(job)])
    result = dict(status='repaired_shared_job_bound_and_released', job_id=job, index=index,
        method=run['method'], proteomes=run['proteomes'], repeat=run['repeat'], source=executor.record(__file__),
        request=request_ref, preparation=executor.record(PREPARATION), execution_scope=executor.SHARED_SCOPE,
        automatic_retry=False, uncontended_timing=False, scientific_timings_admitted=False,
        native_phase='requires_live_observation')
    save(WORK / f'launch_{suffix}.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['resolve', 'prepare', 'launch'])
    parser.add_argument('--index', type=int)
    args = parser.parse_args()
    if args.action == 'launch' and args.index is None:
        parser.error('launch requires --index')
    result = launch(args.index) if args.action == 'launch' else globals()[args.action]()
    print(json.dumps({key: result[key] for key in ('status', 'job_id', 'index', 'session', 'next_index') if key in result}, indent=2))
