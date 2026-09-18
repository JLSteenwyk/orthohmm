"""Freeze fresh corrected-input OrthoHMM and OrthoFinder primary commands."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_simulation_methods import commands, CORE_COMMIT
from benchmark_tools.prepare_scaling_commands import native_orthohmm
from benchmark_tools.run_simulation_methods import read_frozen, verify_environment, execution_environment

BASELINE_SHA = "bf677728cf9312edd9755d9ac1ceefbce3abb8c3848d43413631dc00ca00eb2f"
STAGING_SHA = "07a890eb816944f946d46039559a6046c2a3664f033eaa0f30f35bb25b9c9ab8"
INVENTORY_SHA = "c86c2d6337928938a5a682de0761cc0d293143a1fc5760fe39e5e0c8065e97ca"


def primary_commands(input_directory, output, baseline, cpus=32):
    if type(cpus) is not int or cpus < 1:
        raise ValueError("Positive integer CPU count required")
    configs = commands({"input": str(input_directory)}, output, Path(baseline["core_root"]),
                       Path(baseline["tool_entrypoints"]["orthohmm_python"]["absolute_path"]),
                       Path(baseline["tool_entrypoints"]["orthofinder"]["absolute_path"]))
    result = {}
    for method in ("orthohmm_high_sensitivity", "orthofinder_full"):
        config = configs[method]
        for flag in (("-t", "-a") if method == "orthofinder_full" else ("--cpu",)):
            if config["argv"].count(flag) != 1:
                raise ValueError("Ambiguous resource option")
            config["argv"][config["argv"].index(flag) + 1] = str(cpus)
        result[method] = {**config, "native_argv": config["argv"] if method == "orthofinder_full" else native_orthohmm(config),
                          "cwd": baseline["core_root"], "search_reuse": False,
                          "accuracy_admitted": False}
    return result


def prepare(root, output, destination, runtime_path, runtime_sha):
    if output.exists() or destination.exists():
        raise FileExistsError("Require new output root and command manifest")
    results = root / "benchmark_tools/results"
    stage_path = results / "qfo_corrected_staging_manifest_20260918.json"
    inventory_path = results / "qfo_corrected_staged_inventory_20260918.json"
    baseline_path = results / "publication_variable_native_methods_20260916.json"
    stage = read_frozen(stage_path, STAGING_SHA)
    inventory = read_frozen(inventory_path, INVENTORY_SHA)
    baseline = read_frozen(baseline_path, BASELINE_SHA)
    runtime = read_frozen(runtime_path, runtime_sha)
    if (inventory["status"] != "corrected_staged_inventory_verified_pending_execution_freeze"
            or inventory["total_sequences"] != 984137 or inventory["proteomes"] != 78
            or inventory["inputs"][0]["sha256"] != STAGING_SHA or baseline["core_commit"] != CORE_COMMIT
            or runtime["status"] != "current_resolution_observed"
            or runtime["baseline"] != record(baseline_path)):
        raise ValueError("Changed staged inputs, core or runtime provenance")
    inputs = stage["input_fastas"]
    if [r["file"] for r in inventory["files"]] != inputs:
        raise ValueError("Staging and direct inventory disagree")
    parents = {Path(r["path"]).parent for r in inputs}
    if len(parents) != 1:
        raise ValueError("Require a single staged input directory")
    input_directory = next(iter(parents))
    expected = {Path(r["path"]).name for r in inputs} | {"staging_manifest.json"}
    if {p.name for p in input_directory.iterdir()} != expected:
        raise ValueError("Corrected input directory changed")
    checked = [record(stage_path), record(inventory_path), record(baseline_path), record(runtime_path),
               *inventory["inputs"], *inputs]
    for item in checked:
        check(item)
    verify_environment(baseline)
    _, resolved = execution_environment(baseline)
    for tool in runtime["child_tools"].values():
        if tool["status"] != "resolved":
            raise ValueError("Missing OrthoFinder child tool")
        check(tool["file"])
    for item in runtime["package_sources"]:
        check(item)
    configs = primary_commands(input_directory, output, baseline)
    report = {"status": "corrected_primary_commands_frozen_unrun", "source": record(__file__),
              "inputs": checked, "baseline": record(baseline_path), "runtime": record(runtime_path),
              "input_directory": str(input_directory), "output_root": str(output),
              "core_commit": CORE_COMMIT, "environment_overrides": baseline["environment_overrides"],
              "prepend_path": baseline["prepend_path"], "resolved_outer_executables": resolved,
              "methods": configs, "resource_plan": {"node": "bizon", "cpus": 32, "memory_gib": 192,
                  "time_limit_hours": 72, "execution": "sequential native runs; shared host, not dedicated timings"},
              "helper_sources": [record(Path(__file__).with_name(name)) for name in
                  ("prepare_simulation_methods.py", "prepare_scaling_commands.py", "benchmark_production.py", "run_simulation_methods.py")],
              "required_remaining_rows": ["orthohmm_satellite_v2", "orthofinder_sequence_only", "sonicparanoid",
                  "proteinortho", "fastoma", "orthomcl_1_4"],
              "required_factorial_cells": [f"p{p}_c{c}_r{r}" for p in (0, 1) for c in (0, 1) for r in (0, 1)],
              "execution_authorized": False, "accuracy_evaluated": False,
              "remaining_gates": ["Pinned executor/runner with preflight and postflight environment and input checks.",
                  "Fresh per-method directories, OrthoFinder input copies, no old groups or search checkpoints.",
                  "Native output, complete gene coverage and graph/tree validation before conversion/scoring.",
                  "Freeze remaining comparators and corrected-input factorial; these two runs do not complete the protocol.",
                  "Reuse corrected HMM evidence only after checkpoint provenance/identity verification; no original-release reuse.",
                  "Runtime/memory descriptive on shared host; DGX timing experiment remains separate."]}
    for item in checked:
        check(item)
    with destination.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "output-root", "manifest", "runtime"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--runtime-sha256", required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.output_root.resolve(), args.manifest.resolve(),
            args.runtime.resolve(), args.runtime_sha256)
