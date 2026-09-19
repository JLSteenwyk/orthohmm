"""Derive a fresh paired panel without changing frozen native workloads."""

import argparse
import hashlib
import json
from pathlib import Path

from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.prepare_frontier_overhead_panel import relocate, ROOT
from benchmark_tools.prepare_ob_candidate_neighborhood import record

PARENT_SHA = "b644e165dbf4d0beabf1cf4d9b6c314de522e3ebd1b91598ebebea99094c8fff"


def build(parent, protocol):
    original = read_pinned(parent, PARENT_SHA)
    if [row["index"] for row in original["runs"]] != list(range(18)):
        raise ValueError("Require complete frozen18-task parent")
    plan = relocate(original, ROOT + "/pressure_frontier_overhead_v2", ROOT + "/dual_collector_overhead_v1")
    plan.update(status="prospective_dual_collector_overhead_plan",
        purpose="dual_collector_incremental_overhead", derived_from=record(parent), source=record(__file__),
        protocol_sha256=hashlib.sha256(protocol.read_bytes()).hexdigest(),
        periodic_collector="benchmark_tools.measure_native_dual_bracket_step.measure",
        boundary_collector="benchmark_tools.measure_frontier_boundary_step.measure",
        native_pressure=True, cache_directory=ROOT + "/dual_collector_overhead_v1",
        execution_authorized=False, scientific_timings_admitted=False)
    plan["limitations"] += ["Fresh complete panel; no selected replacement pairs.",
        "Dual periodic versus pressure-enabled boundary collector; native commands unchanged.",
        "Numerical overhead budgets and environmental uncertainty must be reported separately."]
    return plan


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("parent", "protocol", "output"):
        parser.add_argument("--"+name, type=Path, required=True)
    args = parser.parse_args()
    result = build(args.parent.resolve(), args.protocol.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
