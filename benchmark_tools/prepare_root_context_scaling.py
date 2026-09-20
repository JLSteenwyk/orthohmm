"""Prepare all original scaling identities with fresh paths, without authorizing execution."""

import argparse
from collections import Counter
from copy import deepcopy
import gzip
import json
from pathlib import Path

from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.prepare_frontier_overhead_panel import relocate
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_root_context_native import ROOT

PROTOCOL_SHA = "2bccc38c7fd830539e1214af070d782cf670485ad6b2e129bea56c9367be0419"
PLAN_SHA = "7096348236f8372ef7f6ad12a3e829eae3dfee812b17c3b0bdd582480538ff1b"
ORDER_SHA = "7972eb5e1224c0f14dc09a36b26d38f777fddfad3b6075cb5d20999a55a5a451"
OVERHEAD_PLAN_SHA = "02475e4a4290664f776d430b041bd65ccfa640221d07f95b1709b53678a91873"
OVERHEAD_AUDIT_SHA = "6da8441923e51a404fb2e4fec4c7302593fc1a1bc7e338b8996b9d2d2f24641a"
OUTPUT_ROOT = ROOT / "scaling_root_context_v1"
METHOD_IDS = ("orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_3_1_5_full")
AMENDMENT_SHA = "74eefc57da826b4767cfc76019c740cc0f7f91ac708253dce74c00a8dd13c649"
PREPARED_V1_SHA = "65e0f850f32d09e0049e7e55c8700637588d3b7857aa59a2a70eb942d6c5b348"


def assemble(original, order, overhead, audit):
    if (type(audit["validated_tasks"]) is not int or audit["validated_tasks"] != 18 or audit["issues"] != []
            or audit["comparison"]["engineering_budget_passed"] is not True
            or audit["scientific_timings_admitted"] is not False):
        raise ValueError("Require complete validated engineering evidence, not timing admission")
    shared_keys = ("core_commit", "environment_overrides", "environment_paths", "unset_environment")
    if any(original[k] != overhead[k] for k in shared_keys):
        raise ValueError("Native environment differs from tested collector panel")
    runs = original["runs"]
    expected = Counter((method, size, repeat) for method in METHOD_IDS
                       for size in (4, 8, 12) for repeat in range(3))
    if (len(runs) != 27 or [r["index"] for r in runs] != list(range(27))
            or any(type(r[k]) is not int for r in runs for k in ("index", "proteomes", "repeat"))
            or Counter((r["method"], r["proteomes"], r["repeat"]) for r in runs) != expected):
        raise ValueError("Require all 27 original method/size/repeat identities")
    for run in runs:
        matching = [d for d in order["datasets"] if d["input_directory"] == run["dataset"]["input_directory"]]
        if (len(matching) != 1 or sorted(matching[0]["inputs_in_native_order"], key=lambda x:x["path"])
                != sorted(run["dataset"]["inputs"], key=lambda x:x["path"])):
            raise ValueError("Native order and input bytes differ")
    if str(original["output_root"]) != str(ROOT / "scaling_native_v1"):
        raise ValueError("Unexpected original output root")
    tasks = [dict(index=r["index"], method=r["method"], proteomes=r["proteomes"], repeat=r["repeat"],
                  run=relocate(r, str(ROOT / "scaling_native_v1"), str(OUTPUT_ROOT))) for r in runs]
    shared = {k:deepcopy(overhead[k]) for k in (*shared_keys, "baseline_sha256", "enumerator",
              "interval_s", "launcher_python", "runtime_manifests")}
    return dict(shared, schema=1, status="prospective_root_context_scaling_prepared_not_authorized",
        protocol_sha256=PROTOCOL_SHA, original_plan_sha256=PLAN_SHA, original_order_sha256=ORDER_SHA,
        overhead_plan_sha256=OVERHEAD_PLAN_SHA, overhead_audit_sha256=OVERHEAD_AUDIT_SHA,
        runs=tasks, orders=deepcopy(order["datasets"]), output_root=str(OUTPUT_ROOT),
        collector=deepcopy(overhead["collectors"]["root_context"]),
        monitor_host=True, host_interval_s=30., native_timeout_s=85800,
        allocation=dict(node="spark-7ff0", partition="spark", cpus=20, memory_gib=96,
                        exclusive=True, time_limit_s=86400, requeue=False, max_concurrent_runs=1),
        waiting_session=dict(remote_timeout_s=86520, local_timeout_s=86550, remote_kill_grace_s=10),
        failure_policy=dict(native="retain_then_continue_only_after_terminal_and_environment_checks",
            infrastructure_or_policy="pause_submissions_retain_unrun", automatic_retry=False),
        reporting=dict(required_repeats_per_method_size=3, partial_medians=False,
                       retain_all_flags=True, overhead_subtraction=False),
        environment_policy=dict(status="unresolved_requires_separate_freeze", service_change_authorized=False),
        remaining_gates=["User decision and frozen service/workload policy", "Tested launcher and independent audit",
            "Full source recipe and deployment verification", "Environmental evidence and contamination handling"],
        execution_authorized=False, scientific_timings_admitted=False, publication_ready=False)


def build(results):
    protocol = results / "DGX_SCALING_REPLACEMENT_PROTOCOL_20260920.md"
    source = record(protocol)
    if source["sha256"] != PROTOCOL_SHA:
        raise ValueError("Prospective protocol changed")
    original = read_pinned(results / "dgx_scaling_commands_20260917.json", PLAN_SHA)
    order = read_pinned(results / "dgx_native_input_order_20260917.json", ORDER_SHA)
    overhead = read_pinned(results / "dgx_root_context_overhead_plan_20260919.json", OVERHEAD_PLAN_SHA)
    path = results / "root_context_overhead_audit_22022_20260920.json.gz"
    audit_source = record(path)
    if audit_source["sha256"] != OVERHEAD_AUDIT_SHA:
        raise ValueError("Overhead audit hash differs")
    with gzip.open(path, "rt") as stream:
        audit = json.load(stream)
    plan = assemble(original, order, overhead, audit)
    check(source)
    check(audit_source)
    return plan


def build_long_run(results):
    plan = build(results)
    retained = read_pinned(results / "dgx_root_context_scaling_plan_20260920.json", PREPARED_V1_SHA)
    if plan != retained:
        raise ValueError("Original prepared plan does not reproduce")
    amendment = record(results / "DGX_SCALING_LONG_RUN_AMENDMENT_20260920.md")
    if amendment["sha256"] != AMENDMENT_SHA:
        raise ValueError("Long-run amendment changed")
    plan.update(status="prospective_root_context_scaling_v2_prepared_not_authorized",
        prepared_v1_sha256=PREPARED_V1_SHA, long_run_amendment_sha256=AMENDMENT_SHA)
    plan["collector"].update(module="benchmark_tools.measure_scaling_root_context", entry="measure")
    plan["host_observation_scope"] = "counter snapshots, not exhaustive process/GPU/device-I/O inventory"
    plan["remaining_gates"].extend(["Matching long-run replay and composed lifecycle validation",
        "Long-run observer memory/overhead assessment; prior four-proteome results are not a general bound"])
    check(amendment)
    return plan


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--long-run-collector", action="store_true")
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    plan = (build_long_run if args.long_run_collector else build)(args.results.resolve())
    with args.output.open("x") as stream:
        json.dump(plan, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
