"""Freeze paired periodic lineage/root-context tasks without changing native work."""

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.prepare_frontier_overhead_panel import relocate, METHODS
from benchmark_tools.prepare_root_context_native import ROOT

PARENT_SHA = "9b99cd810aaf7e040cda0242dd9b6d1da82bcb231c22657b8dc0108c4c4178ad"
PROTOCOL_SHA = "079ea021eb9f180143f87a1099d6a852c423f21554f95403b81367970a79220e"
PRIOR_AUDIT_SHA = "9b2a243a13cd83ab431025a7f24d5b166672d50cbaa581697a2eba2557e54b23"
OUTPUT_ROOT = ROOT / "root_context_overhead_v1"
BLOCKS = ((0, 1, 2), (1, 2, 0), (2, 0, 1))


def build(parent_path, protocol_path):
    parent = read_pinned(parent_path, PARENT_SHA)
    if hashlib.sha256(protocol_path.read_bytes()).hexdigest() != PROTOCOL_SHA:
        raise ValueError("Root-context overhead protocol differs")
    if ([row["method"] for row in parent["runs"]] != list(METHODS)
            or parent["native_timeout_s"] != 900 or parent["interval_s"] != 1.):
        raise ValueError("Unexpected parent native inventory")
    for index, task in enumerate(parent["runs"]):
        if (task["index"] != index or task["run"]["proteomes"] != 4
                or task["run"]["dataset"]["proteins"] != 73266):
            raise ValueError("Unexpected parent native task")
    runs = []
    for block, methods in enumerate(BLOCKS):
        arms = ("root_context", "lineage") if block == 1 else ("lineage", "root_context")
        for position, method_index in enumerate(methods):
            task = parent["runs"][method_index]
            pair = block * 3 + position
            for arm in arms:
                index = len(runs)
                run = relocate(task["run"], str(ROOT / "root_context_native_v1" / f"run_{method_index:02d}"),
                               str(OUTPUT_ROOT / f"run_{index:02d}"))
                runs.append(dict(index=index, block=block, pair=pair, arm=arm,
                    method=task["method"], native_parent_index=method_index, run=run))
    shared = {key: deepcopy(parent[key]) for key in (
        "baseline_sha256", "core_commit", "enumerator", "environment_overrides", "environment_paths",
        "interval_s", "launcher_python", "native_timeout_s", "order", "runtime_manifests", "unset_environment")}
    return dict(shared, status="prospective_root_context_incremental_overhead", schema=1,
        parent_plan_sha256=PARENT_SHA, protocol_sha256=PROTOCOL_SHA, prior_native_audit_sha256=PRIOR_AUDIT_SHA,
        runs=runs, output_root=str(OUTPUT_ROOT), failure_policy="stop_after_failure_retain_unrun",
        collectors=dict(lineage=dict(module="benchmark_tools.measure_native_lineage_step", entry="measure",
                                     native_report="lineage_report.json", supplementary_report=None),
            root_context=dict(module="benchmark_tools.measure_native_root_context", entry="measure_native_run",
                              native_report="lineage_report.json", supplementary_report="root_context_report.json")),
        monitor_host=True, host_interval_s=30.,
        allocation=dict(node="spark-7ff0", partition="spark", cpus=20, memory_gib=96, exclusive=True,
                        time_limit_s=18000, requeue=False, sequential_steps=True),
        waiting_session=dict(remote_timeout_s=18120, local_timeout_s=18150, remote_kill_grace_s=10),
        overhead_statistic="root_context_native_wall_s / lineage_native_wall_s - 1",
        engineering_budget=dict(every_pair_max=0.10, per_method_median_max=0.05, required_pairs_per_method=3),
        execution_authorized=False, scientific_timings_admitted=False, publication_ready=False,
        limitations=["Both arms are periodic; only supplementary root-context collection differs.",
            "Fixed order is not randomized or perfectly balanced; elapsed-time changes are descriptive.",
            "No flags removed, partial-method medians, selective retries or timing corrections.",
            "Passing engineering budgets does not establish total overhead, isolation or scientific timing eligibility."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    results = Path(__file__).resolve().parent / "results"
    plan = build(results / "dgx_root_context_native_plan_20260919.json",
                 results / "ROOT_CONTEXT_OVERHEAD_PROTOCOL_20260919.md")
    with args.output.open("x") as stream:
        json.dump(plan, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
