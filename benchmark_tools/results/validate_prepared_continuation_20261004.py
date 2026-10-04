"""Retain a source-pinned continuation regression without launching scientific work."""
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


def validate():
    work = ROOT / 'benchmarks/work/threadripper_shared_execution_20261003/prepared_continuation_validation_20261004'
    public = ROOT / 'benchmark_tools/results/threadripper_prepared_continuation_validation_20261004.json'
    if work.exists() or public.exists():
        raise FileExistsError('Existing validation must be retained')
    test = ROOT / 'tests/unit/test_continue_shared_prepared_panel.py'
    helpers = [ROOT / f'benchmark_tools/results/{name}_shared_prepared_panel_20261004.py'
        for name in ('continue', 'review', 'snapshot')]
    top = [ROOT / f'benchmark_tools/{name}.py' for name in (
        'run_threadripper_scaling', 'manage_threadripper_environment_worker',
        'threadripper_environment_worker', 'bind_threadripper_panel_history',
        'threadripper_panel_progress')]
    sources = [record(path) for path in [test, *helpers, *top, Path(__file__)]]
    work.mkdir()
    log, junit = work / 'regression.log', work / 'regression.xml'
    command = [sys.executable, '-B', '-m', 'pytest', '-q', '--tb=short',
        '--junitxml=' + str(junit), str(test)]
    env = os.environ.copy()
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    started = time.time_ns()
    with log.open('x') as stream:
        result = subprocess.run(command, cwd=ROOT, env=env, stdin=subprocess.DEVNULL,
            stdout=stream, stderr=subprocess.STDOUT, timeout=120)
    finished = time.time_ns()
    for pin in sources:
        check(pin)
    tree = ET.parse(junit).getroot()
    counts = {key: sum(int(row.attrib.get(key, 0)) for row in tree.iter('testsuite'))
        for key in ('tests', 'failures', 'errors', 'skipped')}
    passed = result.returncode == 0 and counts['tests'] == 54 and not any(
        counts[key] for key in ('failures', 'errors', 'skipped'))
    report = dict(status='prepared_continuation_regression_passed' if passed else 'continuation_regression_failed',
        command=command, returncode=result.returncode, counts=counts,
        started_unix_ns=started, finished_unix_ns=finished,
        evidence=[*sources, record(log), record(junit)], native_runs_started=False,
        scientific_execution_authorized=False, scientific_timings_admitted=False,
        limitations=['Software regression only; production preparation/handoff remains to be observed.',
            'Existing 655-test kernel regression and runtime refresh are retained, not repeated.',
            'Both aborts remain excluded; resource, response and freshness bounds are unchanged.'])
    save(work / 'validation.json', report)
    save(public, report)
    if not passed:
        raise ValueError('Continuation regression did not pass; preserve failure evidence')
    print(json.dumps(dict(validation=record(public), counts=counts), indent=2))
    return report


if __name__ == '__main__':
    validate()
