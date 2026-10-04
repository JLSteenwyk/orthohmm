"""Retain a terminal scheduler record; never launch, retry or admit native work."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from benchmark_tools.results import continue_shared_prepared_panel_20261004 as panel
from benchmark_tools.verify_threadripper_controller import validate


def capture(index):
    if type(index) is not int or not 21 <= index < 27:
        raise ValueError('Require a prepared-panel identity')
    launch_ref = panel.executor.record(panel.WORK / f'launch_{index:02d}.json')
    launch = panel.executor.read(launch_ref)
    panel.executor.expect(launch, dict(index=index,
        status='repaired_shared_job_bound_and_released', execution_scope=panel.executor.SHARED_SCOPE,
        automatic_retry=False, uncontended_timing=False))
    job = launch['job_id']
    request_ref = launch['request']
    request = panel.executor.read(request_ref)
    panel.executor.expect(request, dict(index=index, job_id=job,
        allocation_cwd=str(ROOT), scheduler_command=str(panel.SCRIPT)))
    output = panel.WORK / f'terminal_controller_{job}.json'
    if output.exists() or output.is_symlink():
        raise FileExistsError('Existing terminal evidence must not be overwritten')
    argv = ['scontrol', 'show', 'job', str(job), '--oneliner']
    started = time.time_ns()
    result = subprocess.run(argv, capture_output=True, text=True, timeout=5, check=True)
    finished = time.time_ns()
    allocation = validate(result.stdout, job, 'terminal', command=str(panel.SCRIPT),
        cwd=str(ROOT), time_limit=panel.executor.TIME_LIMIT, allocation_mode='shared')
    if allocation['fields'].get('Comment') != request_ref['sha256']:
        raise ValueError('Terminal controller request digest differs')
    for pin in (launch_ref, request_ref):
        panel.executor.check(pin)
    controller = dict(command=argv, returncode=result.returncode, stdout=result.stdout,
        stderr=result.stderr, started_unix_ns=started, finished_unix_ns=finished,
        capture_source=panel.executor.record(__file__), launch=launch_ref, request=request_ref,
        native_review_completed=False, next_submission_authorized=False,
        scientific_timings_admitted=False)
    return panel.save(output, controller)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index', required=True, type=int)
    args = parser.parse_args()
    print(json.dumps(capture(args.index), indent=2))
