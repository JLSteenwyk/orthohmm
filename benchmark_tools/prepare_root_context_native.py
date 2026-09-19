"""Derive native root-context diagnostics without changing frozen method commands."""

import argparse
import hashlib
import json
from pathlib import Path

from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.prepare_frontier_overhead_panel import relocate, METHODS

ROOT = Path("/home/jlsteenwyk/projects/orthohmm-publication")
PARENT_SHA = "714ef04458904e1c01526d171d0b2fde658ac135bd416af49ab107dc5ed14bfe"
PROTOCOL_SHA = "d89c15df3648dd791b0a897ef1692952fcc58b4a376c3b1c8cf1dae496e46015"


def build(parent_path, protocol_path):
    parent = read_pinned(parent_path, PARENT_SHA)
    if hashlib.sha256(protocol_path.read_bytes()).hexdigest() != PROTOCOL_SHA:
        raise ValueError("Native root-context protocol changed")
    if (len(parent["runs"]) != 3 or [row["method"] for row in parent["runs"]] != list(METHODS)
            or parent["native_timeout_s"] != 900 or parent["interval_s"] != 1.):
        raise ValueError("Unexpected parent diagnostic inventory")
    runs = []
    for index, task in enumerate(parent["runs"]):
        if (task["index"] != index or task["run"]["proteomes"] != 4
                or task["run"]["dataset"]["proteins"] != 73266):
            raise ValueError("Unexpected frozen native task")
        run = relocate(task["run"], str(ROOT / "lineage_native_v1" / f"run_{index:02d}"),
                       str(ROOT / "root_context_native_v1" / f"run_{index:02d}"))
        runs.append(dict(task, run=run))
    return dict(parent, status="prospective_root_context_native_diagnostic", runs=runs,
        parent_plan_sha256=PARENT_SHA, protocol_sha256=PROTOCOL_SHA,
        collector=dict(module="benchmark_tools.measure_native_root_context", entry="measure_native_run",
            native_report="lineage_report.json", supplementary_report="root_context_report.json",
            monitor_host=True, host_interval_s=30.),
        allocation=dict(node="spark-7ff0", cpus=20, memory_gib=96, exclusive=True,
                        time_limit_s=3600, requeue=False, sequential_steps=True),
        failure_policy="stop_after_failure_retain_unrun", execution_authorized=False,
        scientific_timings_admitted=False, publication_ready=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    folder = Path(__file__).resolve().parent / "results"
    plan = build(folder / "dgx_lineage_native_plan_20260919.json", folder / "ROOT_CONTEXT_NATIVE_PROTOCOL_20260919.md")
    with args.output.open("x") as stream:
        json.dump(plan, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
