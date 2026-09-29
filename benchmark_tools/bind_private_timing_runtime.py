"""Inventory a prospective private deployment without authorizing timing."""

import argparse
import json
from pathlib import Path

from benchmark_tools.build_private_timing_environment import record, write
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.snapshot_runtime_trees import inventory


def root_partition(roots):
    retired, retained = [], []
    for value in sorted(set(roots)):
        path = Path(value)
        if any(path == base or path.is_relative_to(base) for base in (
                Path("/home/bizon/anaconda3"), Path("/home/bizon/.local/lib/python3.10/site-packages"),
                Path("/home/bizon/.cache/matplotlib"))):
            retired.append(value)
        else:
            retained.append(value)
    return retained, retired


def bind(repo, output):
    if output.exists():
        raise FileExistsError(output)
    results = repo / "benchmark_tools/results"
    previous = results / "threadripper_runtime_binding_v4_20260928.json"
    prior = read_frozen(previous, "c6e9225f43604adea6b29adbf5c1f51d3e72ff2f99583f44e42da1939c4cffd6")
    amendment_path = results / "threadripper_private_deployment_20260928.json"
    amendment = json.loads(amendment_path.read_text())
    private_path = results / "threadripper_private_trees_20260928.json"
    private = json.loads(private_path.read_text())
    private_inventory = read_frozen(Path(private["inventory"]["path"]), private["inventory"]["sha256"])
    baseline = read_frozen(Path(amendment["baseline"]["path"]), amendment["baseline"]["sha256"])
    read_frozen(Path(amendment["command_plan"]["path"]), amendment["command_plan"]["sha256"])
    historical = [read_frozen(Path(p), sha) for p, sha in prior["runtime_specs"]]
    roots, retired = root_partition([p for item in historical for p in item["roots"]])
    roots = sorted(set(roots) | {str(p) for p in (repo / "benchmark_tools").glob("*.py")})
    output.mkdir(parents=True)
    write(output / "started.json", dict(source=record(Path(__file__)), previous=record(previous),
          amendment=record(amendment_path), private=record(private_path), roots=roots, retired_roots=retired))
    native = inventory(roots)
    write(output / "native_os_helpers.json", native)
    current_private = inventory(private_inventory["roots"])
    if current_private != private_inventory:
        write(output / "private_drift.json", current_private)
        raise ValueError("Private runtime changed since validated snapshot")
    all_rows = native["records"] + current_private["records"]
    if len({r["path"] for r in all_rows}) != len(all_rows):
        raise ValueError("Overlapping inventory scopes")
    by_path = {r["path"]: r for r in all_rows}
    checked = []
    for ref in (baseline["core_sources"] + baseline["adapter_sources"]
                + baseline["orthofinder_distribution"] + list(baseline["tool_entrypoints"].values())):
        actual = by_path.get(ref["absolute_path"])
        if actual is None:
            raise ValueError("Missing baseline file: " + ref["absolute_path"])
        size = actual.get("bytes", actual.get("target_bytes"))
        sha = actual.get("sha256", actual.get("target_sha256"))
        if (size, sha) != (ref["bytes"], ref["sha256"]):
            raise ValueError("Baseline file changed: " + ref["absolute_path"])
        checked.append(ref["absolute_path"])
    manifests = [record(output / "native_os_helpers.json"), private["inventory"]]
    controller = repo / "benchmarks/work/threadripper_private_controller_20260928/venv/bin/python"
    result = dict(status="prospective_private_timing_runtime_bound", source=record(Path(__file__)),
        supersedes=record(previous), amendment=record(amendment_path), baseline=amendment["baseline"],
        command_plan=amendment["command_plan"], controller_python=record(controller),
        runtime_specs=[[r["path"], r["sha256"]] for r in manifests], runtime_manifests=manifests,
        generator_source=record(repo / "benchmark_tools/snapshot_runtime_trees.py"),
        records=len(all_rows), baseline_paths=sorted(set(checked)), retired_roots=retired,
        external_symlinks=native["external_symlinks"], scientific_execution_authorized=False,
        limitations=["Retired shared/user Python roots require private lookup validation, not presumed absence of use.",
                     "Fresh OS/helper inventory is a new baseline, not proof of equality to the historical runtime.",
                     "Repeated native/controller lookup, collector integration and timing eligibility remain unvalidated."])
    write(output / "binding.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(bind(args.repo.resolve(), args.output.resolve())["status"])
