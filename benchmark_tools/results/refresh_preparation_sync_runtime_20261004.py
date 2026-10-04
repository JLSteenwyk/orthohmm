"""Refresh only the five tested preparation/history helper identities; no inference."""

import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from benchmark_tools.bind_private_timing_runtime import bind
from benchmark_tools.check_threadripper_runtime import check_manifests
from benchmark_tools.inspect_native_python_lookup import compare_lookup
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.probe_dgx_step_separation import save

WORK = ROOT / 'benchmarks/work/threadripper_private_binding_preparation_sync_20261004'
OLD = ROOT / 'benchmark_tools/results/threadripper_private_lookup_prenative_20261004.json'
OLD_SHA = '37901e9a28db2c717096dd1e92de7318d3ad52c528392e35cb595191eea9d1ed'
PUBLIC = ROOT / 'benchmark_tools/results/threadripper_private_lookup_preparation_sync_20261004.json'
CHANGED = {'threadripper_environment_worker.py', 'manage_threadripper_environment_worker.py',
           'run_threadripper_scaling.py', 'threadripper_panel_progress.py', 'bind_threadripper_panel_history.py'}


def refresh():
    if PUBLIC.exists():
        raise FileExistsError(PUBLIC)
    old = read_frozen(OLD, OLD_SHA)
    prior_binding = read_frozen(Path(old['binding']['path']), old['binding']['sha256'])
    print('Inventorying the tested helper repair', flush=True)
    binding = bind(ROOT, WORK)
    binding_ref = record(WORK / 'binding.json')
    for key in ('baseline', 'controller_python', 'command_plan', 'baseline_paths', 'retired_roots'):
        if binding[key] != prior_binding[key]:
            raise ValueError('Scientific/private binding changed: ' + key)
    if binding['runtime_specs'][1:] != prior_binding['runtime_specs'][1:]:
        raise ValueError('Private runtime manifests changed')
    manifests = [read_frozen(Path(value['runtime_specs'][0][0]), value['runtime_specs'][0][1])
        for value in (prior_binding, binding)]
    before, after = [{row['path']:row for row in manifest['records']} for manifest in manifests]
    changes = dict(added=sorted(after.keys()-before.keys()), removed=sorted(before.keys()-after.keys()),
        changed=sorted(path for path in before.keys() & after.keys() if before[path] != after[path]))
    expected = sorted(str(ROOT / 'benchmark_tools' / name) for name in CHANGED)
    if changes != dict(added=[], removed=[], changed=expected):
        raise ValueError('Runtime has differences beyond the five tested helpers: ' + json.dumps(changes))
    repair_ref = record(ROOT / 'benchmark_tools/results/threadripper_preparation_sync_validation_20261004.json')
    repair = read_frozen(Path(repair_ref['path']), repair_ref['sha256'])
    repaired = {pin['path']:pin for pin in repair['evidence']}
    for path in expected:
        if record(path) != repaired[path]:
            raise ValueError('Helper differs from its captured successful regression')
    delta = dict(previous=prior_binding['runtime_manifests'][0], current=binding['runtime_manifests'][0],
        changes=changes, non_helper_drift=False, repair_validation=repair_ref)
    save(WORK / 'helper_delta.json', delta)
    print('Exactly five helper changes; scientific/private identities preserved', flush=True)
    baseline = binding['baseline']
    inspector = ROOT / 'benchmark_tools/inspect_native_python_lookup.py'
    if record(inspector) != old['source']:
        raise ValueError('Native lookup inspector changed')
    checks = {}
    started = time.monotonic()
    checks['before'] = check_manifests(binding['runtime_specs'])
    checks['before_wall_s'] = time.monotonic()-started
    command = [binding['controller_python']['path'], '-B', str(inspector),
        '--baseline', baseline['path'], '--baseline-sha256', baseline['sha256'],
        '--binding', binding_ref['path'], '--binding-sha256', binding_ref['sha256'],
        '--output', str(WORK / 'lookup_current')]
    env = os.environ.copy()
    for name in ('PYTHONHOME', 'LD_PRELOAD', 'LD_LIBRARY_PATH', 'LD_AUDIT'):
        env.pop(name, None)
    env.update(PYTHONNOUSERSITE='1', PYTHONHASHSEED='0', PYTHONDONTWRITEBYTECODE='1',
        PYTHONPYCACHEPREFIX=str(WORK / 'controller_cache'))
    with (WORK / 'lookup_process.log').open('x') as handle:
        process = subprocess.run(command, env=env, stdin=subprocess.DEVNULL,
            stdout=handle, stderr=subprocess.STDOUT, timeout=240)
    save(WORK / 'lookup_process.json', dict(command=command, exit_code=process.returncode))
    if process.returncode or (WORK / 'controller_cache').exists():
        raise ValueError('Private-controller startup failed or wrote bytecode')
    interpreters = {}
    evidence = [record(OLD), binding_ref, baseline, record(inspector), repair_ref]
    for name in ('orthohmm', 'orthofinder'):
        prior_ref = old['interpreters'][name]['reports'][-1]
        observed_ref = record(WORK / 'lookup_current' / (name + '.json'))
        previous = read_frozen(Path(prior_ref['path']), prior_ref['sha256'])
        observed = read_frozen(Path(observed_ref['path']), observed_ref['sha256'])
        interpreters[name] = dict(reports=[prior_ref, observed_ref], modules=len(observed['modules']),
            comparison=compare_lookup(previous, observed))
        evidence += [prior_ref, observed_ref]
    started = time.monotonic()
    checks['after'] = check_manifests(binding['runtime_specs'])
    checks['after_wall_s'] = time.monotonic()-started
    for pin in evidence:
        check(pin)
    report = dict(status='native_lookup_repeated_identity_match', baseline=baseline, binding=binding_ref,
        interpreters=interpreters, source=record(inspector), supersedes=record(OLD),
        scientific_execution_authorized=False,
        current_validation=dict(source=record(__file__), runtime_checks=checks,
            delta=record(WORK / 'helper_delta.json'), process=record(WORK / 'lookup_process.json'),
            log=record(WORK / 'lookup_process.log')),
        limitations=old['limitations'])
    save(WORK / 'lookup.json', report)
    save(PUBLIC, report)
    print(json.dumps(dict(lookup=record(PUBLIC), records=binding['records'],
        modules={name:value['modules'] for name,value in interpreters.items()}), indent=2), flush=True)
    return report


if __name__ == '__main__':
    refresh()
