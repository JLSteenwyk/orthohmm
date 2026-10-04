"""Capture joint preparation-barrier and retained collector regression; never launch jobs."""

import json
import os
from pathlib import Path
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save

MODULES = ('test_threadripper_environment_worker', 'test_manage_threadripper_environment_worker',
    'test_run_threadripper_scaling', 'test_threadripper_panel_progress',
    'test_bind_threadripper_panel_history', 'test_periodic_host_observer',
    'test_measure_threadripper_scaling', 'test_review_threadripper_process_stream',
    'test_verify_threadripper_controller', 'test_review_shared_prenative_failure',
    'test_shared_deadline_failure_review')
HELPERS = ('threadripper_environment_worker', 'manage_threadripper_environment_worker',
    'run_threadripper_scaling', 'threadripper_panel_progress', 'bind_threadripper_panel_history',
    'periodic_host_observer', 'measure_threadripper_scaling',
    'review_threadripper_process_stream', 'verify_threadripper_controller')
REQUIRED_TESTS = ('test_preparation_barrier_precedes_release_and_retains_bound_receipt',
    'test_slow_preparation_finishes_before_release_window_begins',
    'test_preparation_timeout_reaps_child_without_starting_native_work',
    'test_dead_prepared_worker_cannot_release_native_work',
    'test_bound_deadline_abort_keeps_failed_verdict_and_missing_resources',
    'test_actual_timeout_has_no_native_resources',
    'test_bound_pre_native_abort_keeps_failure_and_missing_resources')


def validate():
    work = ROOT / 'benchmarks/work/threadripper_shared_execution_20261003/preparation_sync_validation_20261004'
    public = ROOT / 'benchmark_tools/results/threadripper_preparation_sync_validation_20261004.json'
    immutable = ROOT / 'benchmark_tools/results/threadripper_immutable_sample_revalidation_20261004.json'
    if work.exists() or public.exists() or immutable.exists():
        raise FileExistsError('Retain the existing validation; no overwrite or silent retry')
    tests = [ROOT / f'tests/unit/{name}.py' for name in MODULES]
    sources = [ROOT / f'benchmark_tools/{name}.py' for name in HELPERS]
    tested = [record(path) for path in [*sources, *tests, Path(__file__)]]
    work.mkdir()
    log, junit = work / 'regression.log', work / 'regression.xml'
    command = [sys.executable, '-B', '-m', 'pytest', '-q', '--tb=short',
        '--junitxml=' + str(junit), *map(str, tests)]
    env = os.environ.copy()
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    started = time.time_ns()
    with log.open('x') as handle:
        result = subprocess.run(command, cwd=ROOT, env=env, stdin=subprocess.DEVNULL,
            stdout=handle, stderr=subprocess.STDOUT, timeout=180)
    finished = time.time_ns()
    for pin in tested:
        check(pin)
    tree = ET.parse(junit).getroot()
    suites = list(tree.iter('testsuite'))
    counts = {key: sum(int(suite.attrib.get(key, 0)) for suite in suites)
        for key in ('tests', 'failures', 'errors', 'skipped')}
    names = {case.attrib['name'].split('[')[0] for case in tree.iter('testcase')}
    missing = sorted(set(REQUIRED_TESTS) - names)
    passed = result.returncode == 0 and counts['tests'] > 0 and not any(
        counts[key] for key in ('failures', 'errors', 'skipped')) and not missing
    report = dict(status='preparation_synchronization_regression_passed' if passed else 'repair_regression_failed',
        returncode=result.returncode, command=command, counts=counts,
        required_tests=list(REQUIRED_TESTS), missing_required_tests=missing,
        started_unix_ns=started, finished_unix_ns=finished,
        evidence=[*tested, record(log), record(junit)],
        scientific_execution_authorized=False, scientific_timings_admitted=False,
        native_runs_started=False, automatic_retry=False,
        limitations=['Current source-bound software regression, not a production handoff or refreshed runtime binding.',
            'Preparation precedes collector startup; response, parked-worker and monitoring bounds are unchanged.',
            'Original failed attempts remain failed, excluded and without native resource endpoints.',
            'Historical source proofs need explicit revalidation/retention, not silent mutation, before history can advance.'])
    save(work / 'validation.json', report)
    save(public, report)
    if not passed:
        raise ValueError('Repair regression failed; retain evidence and do not advance')
    save(immutable, dict(report, status='immutable_initial_observation_regression_passed',
        joint_validation=record(public), purpose='Revalidate retained first-sample controls after preparation-barrier/source changes'))
    print(json.dumps(dict(validation=record(public), immutable_sample_revalidation=record(immutable), counts=counts), indent=2))
    return report


if __name__ == '__main__':
    validate()
