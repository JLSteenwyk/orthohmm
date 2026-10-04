"""Capture a running repaired panel identity without admitting a final timing."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from benchmark_tools.results import continue_shared_prepared_panel_20261004 as launcher
from benchmark_tools import run_threadripper_scaling as executor
from benchmark_tools.verify_threadripper_controller import validate


def capture(index):
    _, _, plan, _ = launcher.load()
    if type(index) is not int or not 21 <= index < len(plan['runs']):
        raise ValueError('Invalid continuation index')
    launch_ref = executor.record(launcher.WORK / f'launch_{index:02d}.json')
    launched = executor.read(launch_ref)
    job = launched['job_id']
    executor.expect(launched, dict(index=index, execution_scope=executor.SHARED_SCOPE,
        status='repaired_shared_job_bound_and_released'))
    output = launcher.RESULTS / f'threadripper_shared_live_{job}.json'
    if output.exists():
        raise FileExistsError(output)
    run = plan['runs'][index]
    directory = Path(run['measurement_directory'])
    session = directory.parent.parent / 'sessions' / f'run_{index:02d}'
    preflight_ref = executor.record(session / 'environment_preflight.json')
    preflight = executor.read(preflight_ref)
    executor.expect(preflight, dict(index=index, job_id=job, decision='passed',
        execution_scope=executor.SHARED_SCOPE, uncontended_timing=False))
    ready_ref = executor.record(directory / 'ready.json')
    ready = executor.read(ready_ref)
    if ready['placement']['affinity'] != list(range(32)) or 'job_' + str(job) + '/step_0/' not in Path(f"/proc/{ready['pid']}/cgroup").read_text():
        raise ValueError('Selected native worker is not live with frozen placement')
    request = executor.read(launched['request'])
    policy_ref = executor.read(request['readiness_review'])['environment_policy']
    handoff = launcher.observer_handoff(session, directory, launched['request'], policy_ref, index, job)
    refs = [*handoff['evidence'], launch_ref, launched['request'], launched['preparation'], preflight_ref, ready_ref,
        executor.record(directory / 'preflight_initial_process_sample.json'),
        *[executor.record(session / name) for name in ('started.json', 'lookup_checks/checked_01.json')],
        *[executor.record(directory / name) for name in ('command.json', 'go.json',
            'environment_release.json', 'release_budget.json')]]
    argv = ['scontrol', 'show', 'job', str(job), '--oneliner']
    started = time.time_ns()
    result = subprocess.run(argv, capture_output=True, text=True, timeout=5, check=True)
    finished = time.time_ns()
    allocation = validate(result.stdout, job, 'running', command=str(launcher.SCRIPT),
        cwd=str(ROOT), time_limit=executor.TIME_LIMIT, allocation_mode='shared')
    if allocation['fields'].get('Comment') != launched['request']['sha256']:
        raise ValueError('Live request digest differs')
    samples = []
    with (directory / 'host_processes.jsonl').open() as stream:
        for line in stream:
            samples.append(json.loads(line))
            if len(samples) == 2:
                break
    gap = None if len(samples) < 2 else samples[1]['snapshot']['started_monotonic_s'] - samples[0]['snapshot']['started_monotonic_s']
    for pin in refs:
        executor.check(pin)
    data = dict(status='repaired_shared_native_observed_running', index=index, job_id=job,
        method=run['method'], proteomes=run['proteomes'], repeat=run['repeat'],
        source=executor.record(__file__), evidence=refs,
        controller=dict(command=argv, returncode=result.returncode, stdout=result.stdout,
            stderr=result.stderr, started_unix_ns=started, finished_unix_ns=finished),
        prepared_handoff=handoff, native_placement=ready['placement'], preflight={key: preflight[key] for key in
            ('available_memory_bytes', 'observed_foreign_average_cores', 'background_competition_recorded')},
        first_observed_process_start_gap_s=gap,
        native_log_excerpt=(directory / 'native.log').read_text()[-2048:],
        execution_scope=executor.SHARED_SCOPE, uncontended_timing=False,
        terminal_review_completed=False, scientific_timings_admitted=False,
        next_submission_authorized=False,
        limitations=['Live snapshot only; obtain a fresh scheduler query on every continuation.',
            'An early sampling interval does not prove whole-run monitoring validity.',
            'Background diagnostics do not quantify causal slowdown or correct timings.',
            'The excluded index-0 measurement and index-17/index-20 pre-native aborts remain excluded; no native attempt is retried.'])
    return launcher.save(output, data)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index', type=int, required=True)
    args = parser.parse_args()
    print(json.dumps(capture(args.index), indent=2))
