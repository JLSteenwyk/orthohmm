import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import prepare_pressure_overhead_panel as prepare
from benchmark_tools import run_dgx_frontier_overhead as module
from benchmark_tools import probe_native_pressure, audit_dgx_pressure
from benchmark_tools.prepare_frontier_overhead_panel import relocate, ROOT
from tests.unit.test_run_dgx_frontier_overhead import save

RESULTS = Path(__file__).resolve().parents[2] / 'benchmark_tools/results'
PARENT = RESULTS / 'dgx_frontier_overhead_plan_20260918.json'
PLAN = RESULTS / 'dgx_pressure_overhead_plan_20260919.json'


def test_plan_only_changes_observation_and_output_locations():
    original = json.loads(PARENT.read_text())
    plan = json.loads(PLAN.read_text())
    assert prepare.build(PARENT) == plan
    assert hashlib.sha256(PLAN.read_bytes()).hexdigest() == module.PRESSURE_PLAN_SHA
    expected = relocate(original['runs'], ROOT + '/frontier_overhead_v1', ROOT + '/pressure_frontier_overhead_v1')
    assert plan['runs'] == expected
    for key in ('runtime_manifests', 'core_commit', 'order', 'resource_plan', 'engineering_budget', 'native_timeout_s'):
        assert plan[key] == original[key]
    assert plan['native_pressure'] is True
    assert plan['execution_authorized'] is False
    assert plan['scientific_timings_admitted'] is False


@pytest.fixture
def pinned(tmp_path):
    paths = [Path(module.__file__), Path(module.periodic.__code__.co_filename),
             Path(module.boundary.__code__.co_filename), PLAN,
             Path(probe_native_pressure.__file__), Path(audit_dgx_pressure.__file__)]
    recipe = tmp_path / 'recipe.json'
    recipe_sha = save(recipe, dict(records=[dict(kind='file', path=str(p.resolve()),
                                                sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths]))
    auth = tmp_path / 'auth.json'
    auth_sha = save(auth, dict(purpose='native_pressure_frontier_incremental_overhead', execution_authorized=True,
        scientific_execution_authorized=False, plan_sha256=module.PRESSURE_PLAN_SHA,
        recipe_sha256=recipe_sha, allowed_indices=list(range(18))))
    return [PLAN, auth, auth_sha, recipe, recipe_sha]


@pytest.mark.parametrize('index', range(18))
def test_pressure_selection(pinned, index):
    plan, auth, row = module.select(*pinned, index, module.PRESSURE_PLAN_SHA)
    assert row == plan['runs'][index]
    assert auth['scientific_execution_authorized'] is False


@pytest.mark.parametrize('fault', ['unknown_plan', 'old_plan', 'old_purpose', 'missing_probe'])
def test_cross_plan_or_unpinned_pressure_rejected(pinned, fault):
    pin = module.PRESSURE_PLAN_SHA
    if fault == 'unknown_plan':
        pin = '0' * 64
    elif fault == 'old_plan':
        pin = module.PLAN_SHA
    else:
        auth = json.loads(pinned[1].read_text())
        if fault == 'old_purpose':
            auth['purpose'] = 'native_frontier_incremental_overhead'
        else:
            recipe = json.loads(pinned[3].read_text())
            recipe['records'].pop()
            pinned[4] = save(pinned[3], recipe)
            auth['recipe_sha256'] = pinned[4]
        pinned[2] = save(pinned[1], auth)
    with pytest.raises(ValueError):
        module.select(*pinned, 0, pin)


@pytest.mark.parametrize('index', [0, 1])
def test_launch_enables_pressure_and_retains_plan_binding(pinned, tmp_path, monkeypatch, index):
    plan, auth, row = module.select(*pinned, index, module.PRESSURE_PLAN_SHA)
    directory = tmp_path / 'run'
    row['run']['measurement_directory'] = str(directory / 'measurement')
    monkeypatch.setattr(module, 'select', lambda *a: (plan, auth, row))
    monkeypatch.setattr(module.os, 'uname', lambda: SimpleNamespace(nodename='spark-7ff0'))
    monkeypatch.setattr(module.os, 'chdir', lambda p: None)
    monkeypatch.setattr(module.sys, 'dont_write_bytecode', True)
    for key,value in dict(SLURM_JOB_ID='123',SLURM_CPUS_PER_TASK='20',SLURM_MEM_PER_NODE='98304',PYTHONHASHSEED='0').items():
        monkeypatch.setenv(key,value)
    for key in ('LD_PRELOAD','LD_LIBRARY_PATH','LD_AUDIT'):
        monkeypatch.delenv(key,raising=False)
    exists = Path.exists
    monkeypatch.setattr(Path,'exists',lambda p: False if str(p)=='/etc/ld.so.preload' or str(p).startswith(plan['cache_directory']) else exists(p))
    def measure(*args, **kwargs):
        collector=args[4]
        assert collector.func is (module.boundary if row['mode']=='boundary' else module.periodic)
        assert collector.keywords == dict(native_pressure=True)
        directory.mkdir()
        return dict(status='command_failed')
    monkeypatch.setattr(module,'measure_run',measure)
    assert module.launch(*pinned,index,module.PRESSURE_PLAN_SHA)['status']=='command_failed'
    result=json.loads((directory/'overhead_task.json').read_text())
    assert result['plan_sha256']==module.PRESSURE_PLAN_SHA
    assert not result['scientific_timings_admitted']
