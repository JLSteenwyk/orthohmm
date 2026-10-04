"""Freeze thirteen missing cost identities; this does not authorize a launch."""

import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from benchmark_tools.run_native_factorial_cost import IDENTITIES, MEMORY, CORE_COMMIT, read, validate_plan
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.snapshot_orthohmm_input_order import snapshot

PINS = dict(
    orthobench=("benchmark_tools/results/orthobench_factorial_prepared_20260916.json", "5c325f4d77865e0c7571fe4bb4df0be0959977a49e1d22169f192518700f9382", "fasta_inputs"),
    qfo_corrected=("benchmarks/results/qfo_corrected_factorial_v1/manifest.json", "d8385c50426e690afd6d32f3c5302e678de6977c9841451d013201b0f75b564a", "input_fastas"))


def prepare(output, panel, runtime_lookup=None, supersedes=None, retained_failure=None):
    if output.exists() or output.is_symlink() or panel.exists() or panel.is_symlink():
        raise FileExistsError("Require a new plan and absent new factorial panel")
    datasets, evidence = {}, []
    for label, (name, sha, field) in PINS.items():
        ref = record(ROOT / name)
        if ref["sha256"] != sha:
            raise ValueError("Original factorial input manifest differs")
        datasets[label] = read(ref)[field]
        evidence.append(ref)
        for pin in datasets[label]:
            check(pin)
    old_lookup = record(ROOT / "benchmark_tools/results/threadripper_private_lookup_preparation_sync_20261004.json")
    if old_lookup["sha256"] != "7e3ac02b46cf1ead426ec19889fc5bd23e6a0e01dbd3737bf374c3ae3180ed1d":
        raise ValueError("Retained runtime lookup differs")
    lookup = old_lookup if runtime_lookup is None else runtime_lookup
    receipt = read(lookup)
    old = read(old_lookup)
    if runtime_lookup is not None:
        prior_binding, current_binding = read(old["binding"]), read(receipt["binding"])
        if (receipt.get("status") != "native_lookup_repeated_identity_match" or receipt["baseline"] != old["baseline"]
                or receipt["source"] != old["source"]
                or any(prior_binding[k] != current_binding[k] for k in
                    ("baseline", "controller_python", "command_plan", "baseline_paths", "retired_roots"))
                or prior_binding["runtime_specs"][1:] != current_binding["runtime_specs"][1:]):
            raise ValueError("New lookup changed the frozen scientific/private deployment")
        evidence.extend([lookup, old_lookup, receipt["current_validation"]["delta"]])
    if (supersedes is None) != (retained_failure is None):
        raise ValueError("A replacement must retain both the old plan and failed attempt")
    if supersedes is not None:
        old_plan, failure = read(supersedes), read(retained_failure)
        validate_plan(old_plan)
        if (failure.get("plan") != supersedes or failure.get("status") != "factorial_attempt_failed_retained"
                or failure.get("index") != 0 or failure.get("wrapper", {}).get("status") != "verified_wrapper_failed"
                or failure["wrapper"].get("error") != "Runtime inventory changed"):
            raise ValueError("Replacement is not the retained runtime-preflight failure")
        evidence.extend([supersedes, retained_failure, record(ROOT / "benchmark_tools/NATIVE_FACTORIAL_RUNTIME_REPAIR_22426.md")])
    baseline = receipt["baseline"]
    frozen = read(baseline)
    if frozen["core_commit"] != CORE_COMMIT:
        raise ValueError("Wrong frozen scientific commit")
    runtime = dict(records=[dict(path=r["absolute_path"], bytes=r["bytes"], sha256=r["sha256"])
                           for r in frozen["core_sources"]])
    orders, probes = {}, {}
    for label, inputs in datasets.items():
        # Empty-name fixture checks filesystem enumeration, never runs inference.
        with tempfile.TemporaryDirectory(prefix="factorial_order_", dir=panel.parent) as directory:
            directory = Path(directory)
            for name in sorted(Path(r["path"]).name for r in inputs):
                (directory / name).touch(exist_ok=False)
            rows = [record(p) for p in directory.iterdir()]
            probes[label] = snapshot(frozen["core_root"], dict(datasets=[dict(proteomes=len(inputs),
                                      input_directory=str(directory), inputs=rows)]), runtime)
            orders[label] = probes[label]["datasets"][0]["native_order"]
    for name in ("NATIVE_FACTORIAL_COST_PROTOCOL_20261004.md", "results/native_factorial_adapter_execution_20261004.json",
                 "results/THREADRIPPER_ASYNC_CALIBRATION_RESULT_22380.md", "results/threadripper_shared_panel_snapshot_20261004_v27/panel.json",
                 "NATIVE_FACTORIAL_EXECUTION_PROTOCOL_20261004.md"):
        evidence.append(record(ROOT / "benchmark_tools" / name))
    helpers = [record(p) for p in sorted((ROOT / "benchmark_tools").glob("*.py"))]
    helpers.append(record(ROOT / "benchmark_tools/run_native_factorial_cost.sh"))
    runs = []
    for i, (dataset, cell) in enumerate(IDENTITIES):
        inputs = datasets[dataset]
        directories = {str(Path(r["path"]).parent) for r in inputs}
        if len(directories) != 1:
            raise ValueError("Input files do not share a dataset directory")
        runs.append(dict(index=i, dataset=dataset, cell=cell, repeat=0, inputs=inputs,
            input_directory=directories.pop(), native_order=orders[dataset],
            input_creation_order=sorted(Path(r["path"]).name for r in inputs),
            output_root=str(panel / f"run_{i:02d}"), proteomes=12 if i < 6 else 78, genes=251378 if i < 6 else 984137))
    plan = dict(schema="native_factorial_cost_plan_v1", root=str(ROOT), panel_root=str(panel),
        core_commit=CORE_COMMIT, source_commit=subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        baseline=baseline, runtime_lookup=lookup, runs=runs, execution_scope="shared_host_matched_resources",
        automatic_retry=False, scientific_execution_authorized=False,
        resources=dict(native_cpu_ids=list(range(32)), slurm_slots=64, memory_bytes=MEMORY, timeout_s=85800,
            sample_period_s=1., host_period_s=30., minimum_available_memory_bytes=MEMORY),
        helper_sources=helpers, evidence=evidence, filename_enumeration_probes=probes,
        limitations=["Plan preparation is not execution authorization; each held Slurm job requires an explicit bound request.",
            "Descriptive one-attempt costs, not replicated causal component effects or isolated efficiency ranks.",
            "Fresh persistent input copies differ in storage scope/order from historical tmpfs runs; no exact-output equality is assumed."])
    if supersedes is not None:
        plan.update(supersedes_plan=supersedes, retained_preflight_failure=retained_failure,
                    replacement_reason="Known packaging-only stale inventory; original attempt retained and no native inference previously ran")
    validate_plan(plan)
    with output.open("x") as handle:
        json.dump(plan, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return plan


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--panel-root", type=Path, required=True)
    parser.add_argument("--runtime-lookup", type=Path)
    parser.add_argument("--runtime-lookup-sha256")
    parser.add_argument("--supersedes-plan", type=Path)
    parser.add_argument("--retained-preflight-failure", type=Path)
    args = parser.parse_args()
    runtime_ref = record(args.runtime_lookup) if args.runtime_lookup else None
    if bool(args.runtime_lookup) != bool(args.runtime_lookup_sha256) or runtime_ref and runtime_ref["sha256"] != args.runtime_lookup_sha256:
        parser.error("A refreshed lookup requires its explicit SHA256")
    prepare(args.output.absolute(), args.panel_root.absolute(), runtime_lookup=runtime_ref,
        supersedes=record(args.supersedes_plan) if args.supersedes_plan else None,
        retained_failure=record(args.retained_preflight_failure) if args.retained_preflight_failure else None)
