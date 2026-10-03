"""Prepare shared-host evidence and release one held, explicitly bound identity."""
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
from benchmark_tools import run_threadripper_scaling as executor
from benchmark_tools.review_threadripper_process_policy import SHARED_SCOPE

RESULTS = ROOT / 'benchmark_tools/results'
WORK = ROOT / 'benchmarks/work/threadripper_shared_execution_20261003'
LOOKUP = RESULTS / 'threadripper_private_lookup_shared_20261003.json'
LOOKUP_SHA = '4d769df1272da767be3cf84e3e2a44d680e5b8a525770e9bb94a9a1b5014777f'
SCRIPT = ROOT / 'benchmark_tools/run_threadripper_shared_scaling.sh'
FILES = dict(recipe=RESULTS / 'threadripper_shared_source_recipe_20261003.json',
    process=RESULTS / 'threadripper_shared_process_policy_20261003.json',
    policy=RESULTS / 'threadripper_shared_environment_policy_20261003.json',
    readiness=RESULTS / 'threadripper_shared_readiness_20261003.json')


def save(path, data):
    with path.open('x') as stream:
        json.dump(data, stream, indent=2, sort_keys=True)
        stream.write('\n')
    return executor.record(path)


def command(name, argv):
    output = WORK / (name + '.json')
    if output.exists():
        raise FileExistsError(output)
    started = time.time_ns()
    result = subprocess.run(argv, cwd=ROOT, text=True, capture_output=True, timeout=30)
    save(output, dict(command=argv, started_unix_ns=started, finished_unix_ns=time.time_ns(),
        returncode=result.returncode, stdout=result.stdout, stderr=result.stderr))
    if result.returncode:
        raise RuntimeError(name + ': ' + result.stderr)
    return result.stdout


def load():
    lookup = executor.read_frozen(LOOKUP, LOOKUP_SHA)
    binding = executor.read(lookup['binding'])
    if os.path.abspath(sys.executable) != binding['controller_python']['path']:
        raise ValueError('Use the frozen private controller')
    executor.check(binding['controller_python'])
    plan_path = RESULTS / 'threadripper_private_commands_20260928.json'
    plan = executor.read_frozen(plan_path, executor.PRIVATE_PLAN_SHA)
    return lookup, binding, plan, executor.record(plan_path)


def prepare():
    if WORK.exists() or any(path.exists() for path in FILES.values()):
        raise FileExistsError('Shared launch preparation already exists')
    lookup, binding, plan, plan_ref = load()
    manifest = executor.read(binding['runtime_manifests'][0])
    by_path = {row['path']: row for row in manifest['records']}
    paths = sorted(set((ROOT / 'benchmark_tools').glob('*.py')) | {SCRIPT})
    sources = [executor.record(path) for path in paths]
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    for ref in sources:
        relative = str(Path(ref['path']).relative_to(ROOT))
        if Path(ref['path']).read_bytes() != subprocess.check_output(['git', 'show', commit + ':' + relative], cwd=ROOT):
            raise ValueError('Uncommitted execution source: ' + relative)
        if relative.endswith('.py'):
            row = by_path[ref['path']]
            if (row['bytes'], row['sha256']) != (ref['bytes'], ref['sha256']):
                raise ValueError('Stale runtime helper: ' + relative)
    calibration_path = RESULTS / 'threadripper_async_calibration_audit_22380.json'
    calibration_ref = executor.record(calibration_path)
    calibration = executor.read(calibration_ref)
    if (calibration['status'] != 'calibration_checks_passed' or calibration['job_id'] != 22380
            or not calibration['checks'] or not all(value is True for value in calibration['checks'].values())):
        raise ValueError('Retained accounting calibration does not pass')
    for ref in calibration['evidence']:
        executor.check(ref)
    evidence = [executor.record(RESULTS / 'PUBLICATION_SHARED_HOST_AMENDMENT_20261003.md'),
        executor.record(ROOT / 'benchmark_tools/PUBLICATION_GOAL_20261003.txt'),
        executor.record(RESULTS / 'THREADRIPPER_SHARED_HOST_EXECUTION_20261003.md'),
        executor.record(RESULTS / 'threadripper_resource_endpoints_20260929.json'),
        calibration_ref, executor.record(LOOKUP), lookup['binding'], executor.record(Path(__file__))]
    for ref in evidence:
        executor.check(ref)
    WORK.mkdir()
    recipe_ref = save(FILES['recipe'], dict(schema='threadripper_executor_recipe_v1', root=str(ROOT),
        sources=sources, source_commit=commit, execution_scope=SHARED_SCOPE,
        scientific_timings_admitted=False))
    boot = Path('/proc/sys/kernel/random/boot_id').read_text().strip()
    process_ref = save(FILES['process'], dict(schema='threadripper_process_policy_v3',
        execution_scope=SHARED_SCOPE, boot_id=boot,
        review_reference='User-authorized shared-host scope, 2026-10-03',
        outside_activity='observed_diagnostic_not_isolation_review'))
    policy_ref = save(FILES['policy'], dict(schema='threadripper_environment_policy_v3',
        decision='reviewed', host='bizon', execution_scope=SHARED_SCOPE,
        plan_sha256=plan_ref['sha256'], process_policy=process_ref, configuration_files=[],
        foreign_cpu_role='diagnostic_only', native_pressure_role='diagnostic_only',
        preflight_pressure_role='diagnostic_only', maximum_foreign_average_cores=0.,
        maximum_pressure_percent=dict(cpu=100., memory=100., io=100.),
        maximum_sample_period_s=35., maximum_pressure_sample_period_s=1.5,
        minimum_available_memory_bytes=128 * 1024**3,
        review_reference='Prospective shared-host resource/annotation policy; no isolation claim',
        evidence=evidence, scientific_timings_admitted=False,
        limitations=['All CPU/PSI magnitudes are diagnostics; limits are not outside-work eligibility.',
            'Outside read errors/churn limit contention estimates; native accounting remains separately validated.',
            '128 GiB available memory is checked at launch, not a guarantee of future host capacity.']))
    ready_ref = save(FILES['readiness'], dict(schema='threadripper_shared_readiness_review_v1',
        decision='passed', execution_scope=SHARED_SCOPE, plan_sha256=plan_ref['sha256'],
        lookup_sha256=LOOKUP_SHA, recipe_sha256=recipe_ref['sha256'], environment_policy=policy_ref,
        observer_accounting_validated=True, environment_policy_frozen=True,
        contention_annotation_required=True, isolation_required=False,
        review_reference='Retained 22380 accounting/cadence checks plus current source/private-runtime binding; actual environmental release is validated by the parked-worker guard',
        evidence=[*evidence, recipe_ref, process_ref, policy_ref],
        native_handoff_established=False, controlled_workload_verified=False,
        observer_causal_slowdown_established=False, scientific_timings_admitted=False))
    preparation = dict(status='shared_source_policy_accounting_readiness_prepared',
        source=executor.record(Path(__file__)), source_commit=commit, recipe=recipe_ref,
        policy=policy_ref, readiness=ready_ref, lookup=executor.record(LOOKUP),
        binding=lookup['binding'], plan=plan_ref, sources=len(sources),
        inherited_accounting_checks=calibration['checks'], native_runs_started=False,
        limitations=['Readiness is scoped accounting/source preparedness, not an environmental handoff or isolated performance pass.'])
    save(WORK / 'preparation.json', preparation)
    return preparation


def launch(index):
    lookup, binding, plan, plan_ref = load()
    preparation = executor.read(executor.record(WORK / 'preparation.json'))
    for key in ('recipe', 'policy', 'readiness', 'lookup', 'binding', 'plan'):
        executor.check(preparation[key])
    if index != 0:
        raise ValueError('Later identities require independently reviewed history; this first-launch source does not create it')
    run = plan['runs'][index]
    request_path = WORK / f'request_run_{index:02d}.json'
    session = Path(run['measurement_directory']).parent.parent / 'sessions' / f'run_{index:02d}'
    if request_path.exists() or session.exists() or Path(run['measurement_directory']).parent.exists():
        raise FileExistsError('Existing request/session/native attempt must be reviewed, not retried')
    raw = command('submission_00', ['sbatch', '--hold', '--parsable', '--partition=gpu',
        '--chdir=' + str(ROOT), '--output=' + str(WORK / 'slurm_run_00-%j.log'),
        '--error=' + str(WORK / 'slurm_run_00-%j.log'), str(SCRIPT),
        '--request', str(request_path), '--request-sha256', 'scheduler-comment'])
    match = re.fullmatch(r'([0-9]+)(?:;[^\s]+)?\n?', raw)
    if match is None:
        raise ValueError('Cannot identify submitted held job; preserve submission before any further action')
    job = int(match[1])
    held = command('held_controller_00', ['scontrol', 'show', 'job', str(job), '--oneliner'])
    fields = dict(re.findall(r'(?<!\S)([A-Za-z][^\s=]*)=([^\s]+)', held))
    required = dict(JobId=str(job), JobState='PENDING', Reason='JobHeldUser', NumCPUs='64',
        NumTasks='1', MinMemoryNode='128G', TimeLimit='1-02:00:00', Requeue='0',
        Command=str(SCRIPT), WorkDir=str(ROOT))
    required['CPUs/Task'] = '64'
    if any(fields.get(key) != value for key, value in required.items()):
        raise ValueError('Held allocation identity/resources differ; do not release')
    request = dict(schema='threadripper_execution_request_v1', job_id=job,
        execution_authorized=True, deployment='private_v2_20260928', execution_scope=SHARED_SCOPE,
        runtime_lookup=preparation['lookup'], plan_sha256=plan_ref['sha256'], lookup_sha256=LOOKUP_SHA,
        index=index, allocation_cwd=str(ROOT), scheduler_command=str(SCRIPT),
        recipe=preparation['recipe'], history=[], readiness_review=preparation['readiness'],
        environment_preflight_path=str(session / 'environment_preflight.json'))
    request_ref = save(request_path, request)
    chosen, _, history, sources = executor.select(request, ROOT, job)
    if chosen != run or history['progress']['index'] != index:
        raise ValueError('Actual job selection differs from frozen identity')
    for ref in [request_ref, *sources]:
        executor.check(ref)
    command('request_comment_00', ['scontrol', 'update', 'JobId=' + str(job),
                                  'Comment=' + request_ref['sha256']])
    pinned = command('pinned_controller_00', ['scontrol', 'show', 'job', str(job), '--oneliner'])
    if ('Comment=' + request_ref['sha256']) not in pinned.split():
        raise ValueError('Scheduler did not preserve actual request hash; leave held')
    command('release_00', ['scontrol', 'release', str(job)])
    report = dict(status='held_job_request_bound_and_released', job_id=job, index=index,
        method=run['method'], proteomes=run['proteomes'], repeat=run['repeat'],
        source=executor.record(Path(__file__)), request=request_ref, preparation=executor.record(WORK / 'preparation.json'),
        execution_scope=SHARED_SCOPE, automatic_retry=False, uncontended_timing=False,
        scientific_timings_admitted=False, native_phase='requires_live_observation')
    save(WORK / 'launch_00.json', report)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'launch'])
    args = parser.parse_args()
    result = prepare() if args.action == 'prepare' else launch(0)
    print(json.dumps({key: result[key] for key in ('status', 'job_id', 'index', 'sources') if key in result}, indent=2))
