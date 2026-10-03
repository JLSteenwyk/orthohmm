"""Reconcile retained resource scopes with the actually executed source recipe."""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from benchmark_tools.results import prepare_shared_threadripper_launch_20261003 as launch
from benchmark_tools import run_threadripper_scaling as executor
from benchmark_tools.derive_threadripper_resources import SCOPES, protocol


def bind():
    parent_path = launch.RESULTS / 'threadripper_resource_endpoints_20260929.json'
    parent_sha = 'dba9602c365184f7e802425b9116b8021175b225d1f2266c44077da924bcfcdc'
    parent = executor.read_frozen(parent_path, parent_sha)
    if parent['primary_scopes'] != SCOPES:
        raise ValueError('Prospective endpoint definitions differ')
    prepared_ref = executor.record(launch.WORK / 'preparation.json')
    prepared = executor.read(prepared_ref)
    recipe = executor.read(prepared['recipe'])
    source_pins = {pin['path']: pin for pin in recipe['sources']}
    sources = [executor.record(Path(pin['path'])) for pin in parent['sources']]
    for pin in sources:
        executor.check(pin)
        if source_pins[pin['path']] != pin:
            raise ValueError('Resource source differs from actually executed recipe')
    output = launch.RESULTS / 'threadripper_resource_endpoints_shared_binding_20261003.json'
    if output.exists():
        raise FileExistsError(output)
    value = dict(parent, sources=sources, plan=prepared['plan'], lookup=prepared['lookup'],
        parent_protocol=executor.record(parent_path), executed_preparation=prepared_ref,
        source=executor.record(__file__), execution_scope=executor.SHARED_SCOPE,
        decision='source_binding_reconciled_without_endpoint_change',
        reconciliation_timing='after_run_00_before_any_additional_native_submission',
        remaining_prerequisites=['actual per-run runtime/resource/environment/output review',
            'fresh capacity and shared-host launch preflight'],
        limitations=['The original source inventory predates asynchronous host observation; its historical bytes are unchanged.',
            'This reconciles source pins against the actually executed recipe, not a prospective claim for new code.',
            'Primary/secondary scopes, failure retention and unavailable endpoints are unchanged.',
            'The original run-00 cadence failure remains failed; this binding does not admit timings or authorize continuation.'])
    ref = launch.save(output, value)
    protocol(output, ref['sha256'])
    return ref


if __name__ == '__main__':
    print(json.dumps(bind(), indent=2))
