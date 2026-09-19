"""Derive fresh pressure-enabled overhead tasks without changing native work."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.prepare_frontier_overhead_panel import relocate, ROOT
from benchmark_tools.prepare_ob_candidate_neighborhood import record

PARENT_SHA = '58e482e6f123bc82de5e98df3c93f4bf0ca3a47d2302ffcaa47290c7d1e20685'


def build(parent):
    original = read_pinned(parent, PARENT_SHA)
    plan = relocate(original, ROOT + '/frontier_overhead_v1', ROOT + '/pressure_frontier_overhead_v1')
    plan.update(status='prospective_native_pressure_frontier_overhead_plan',
                derived_from=record(parent), source=record(__file__), native_pressure=True,
                cache_directory=ROOT + '/pressure_frontier_overhead_v1',
                purpose='native_pressure_frontier_incremental_overhead')
    plan['limitations'] += [
        'Both arms collect native-step PSI; pressure is diagnostic, not an interference exclusion threshold.',
        'Entire new panel, not selective replacements for historical failed or missing pairs.',
        'Terminal scheduler records must be collected from the controller before expiry.',
        'Frozen historical results and the original scientific timing admission flags are unchanged.']
    return plan


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    plan = build(args.parent.resolve())
    with args.output.open('x') as stream:
        json.dump(plan, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
