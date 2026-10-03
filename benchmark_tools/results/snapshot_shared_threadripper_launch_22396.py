"""Retain a live launch observation, not a terminal timing or admission review."""
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.verify_threadripper_controller import validate


def capture():
    work = ROOT / 'benchmarks/work/threadripper_shared_execution_20261003'
    session = ROOT / 'benchmarks/results/threadripper_scaling_v1/sessions/run_00'
    measurement = ROOT / 'benchmarks/results/threadripper_scaling_v1/run_00/measurement'
    output = ROOT / 'benchmark_tools/results/threadripper_shared_launch_22396.json'
    if output.exists():
        raise FileExistsError(output)
    paths = [work / name for name in ('preparation.json', 'submission_00.json',
        'held_controller_00.json', 'request_run_00.json', 'request_comment_00.json',
        'pinned_controller_00.json', 'release_00.json', 'launch_00.json')]
    paths += [session / name for name in ('started.json', 'environment_preflight.json',
        'environment_worker_evidence.json', 'lookup_checks/checked_01.json')]
    paths += [measurement / name for name in ('command.json', 'ready.json', 'go.json',
        'environment_release.json', 'release_budget.json')]
    refs = [record(path) for path in paths]
    launch = json.loads((work / 'launch_00.json').read_text())
    preflight = json.loads((session / 'environment_preflight.json').read_text())
    ready = json.loads((measurement / 'ready.json').read_text())
    expected = dict(job_id=22396, index=0, execution_scope='shared_host_matched_resources')
    if any(launch.get(key) != value or preflight.get(key) != value
           for key, value in expected.items()) or preflight['decision'] != 'passed':
        raise ValueError('Launch and successful preflight identities differ')
    if ready['placement']['affinity'] != list(range(32)):
        raise ValueError('Native affinity differs')
    argv = ['scontrol', 'show', 'job', '22396', '--oneliner']
    started = time.time_ns()
    result = subprocess.run(argv, capture_output=True, text=True, timeout=5, check=True)
    finished = time.time_ns()
    allocation = validate(result.stdout, 22396, 'running',
        command=str(ROOT / 'benchmark_tools/run_threadripper_shared_scaling.sh'),
        cwd=str(ROOT), time_limit='1-02:00:00', allocation_mode='shared')
    if allocation['fields']['Comment'] != launch['request']['sha256']:
        raise ValueError('Live scheduler request binding differs')
    if 'job_22396/step_0/' not in Path(f"/proc/{ready['pid']}/cgroup").read_text():
        raise ValueError('Recorded native worker is not live in the selected step')
    for ref in refs:
        check(ref)
    data = dict(status='shared_native_released_job_observed_running', **expected,
        method=launch['method'], proteomes=launch['proteomes'], repeat=launch['repeat'],
        source=record(__file__), evidence=refs, controller=dict(command=argv,
            returncode=result.returncode, stdout=result.stdout, stderr=result.stderr,
            started_unix_ns=started, finished_unix_ns=finished),
        native_placement=ready['placement'],
        preflight={key: preflight[key] for key in ('decision', 'available_memory_bytes',
            'observed_foreign_average_cores', 'background_competition_recorded',
            'foreign_cpu_used_for_eligibility', 'preflight_pressure_limits_used')},
        scientific_timings_admitted=False, terminal_review_completed=False,
        next_submission_authorized=False, uncontended_timing=False,
        limitations=['Live launch snapshot only; future state requires a fresh scheduler query.',
            'Background CPU estimate is preflight diagnostic, not whole-run contention or a slowdown correction.',
            'Native completion, post-runtime identity, whole-run environment/resource replay and output validation remain pending.',
            'Legacy helper quiet-host wording does not supersede the explicit shared-host policy.'])
    with output.open('x') as stream:
        json.dump(data, stream, indent=2, sort_keys=True)
        stream.write('\n')
    return record(output)


if __name__ == '__main__':
    print(json.dumps(capture(), indent=2))
