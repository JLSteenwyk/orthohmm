"""Package exact scaling inputs and tracked Python workflows for a dedicated host."""

import argparse
import json
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_scaling_commands import INPUT_SHA
from benchmark_tools.prepare_scaling_inputs import planned_runs
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_simulation_tree_mode_control import METHOD_SHA


def portable_inputs(inputs):
    if inputs["planned_runs"] != planned_runs() or len(inputs["ordered_proteomes"]) != 12:
        raise ValueError("Changed transfer panel")
    files, names = [], set()
    for row in inputs["ordered_proteomes"]:
        item = row["input"]
        name = Path(item["path"]).name
        if not name or name in names or name in {".", ".."}:
            raise ValueError("Ambiguous portable input basename")
        names.add(name)
        files.append({"path": "inputs/" + name, "bytes": item["bytes"], "sha256": item["sha256"],
                      "proteins": row["proteins"], "sequence_characters": row["sequence_characters"]})
    datasets = []
    if [d["proteomes"] for d in inputs["datasets"]] != [4, 8, 12]:
        raise ValueError("Incomplete nested transfer datasets")
    for dataset in inputs["datasets"]:
        chosen = files[:dataset["proteomes"]]
        source = {Path(r["path"]).name: (r["bytes"], r["sha256"]) for r in dataset["inputs"]}
        expected = {PurePosixPath(r["path"]).name: (r["bytes"], r["sha256"]) for r in chosen}
        if source != expected or len(dataset["inputs"]) != len(chosen):
            raise ValueError("Nested dataset differs from frozen proteome order")
        if (sum(r["proteins"] for r in chosen) != dataset["proteins"] or
                sum(r["sequence_characters"] for r in chosen) != dataset["sequence_characters"]):
            raise ValueError("Portable input count mismatch")
        datasets.append({"proteomes": dataset["proteomes"], "inputs": [r["path"] for r in chosen],
                         "proteins": dataset["proteins"], "sequence_characters": dataset["sequence_characters"]})
    return files, datasets


def copy_verified(source, destination, expected):
    before = record(source)
    if (before["bytes"], before["sha256"]) != (expected["bytes"], expected["sha256"]):
        raise ValueError("Transfer source changed")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    shutil.copy2(source, destination)
    after = record(destination)
    if (after["bytes"], after["sha256"]) != (expected["bytes"], expected["sha256"]):
        raise ValueError("Transfer copy differs")
    check(before)


def prepare(root, output):
    if output.exists():
        raise FileExistsError(output)
    inputs = read_frozen(root / "benchmark_tools/results/publication_scaling_inputs_20260916.json", INPUT_SHA)
    baseline = read_frozen(root / "benchmark_tools/results/publication_variable_native_methods_20260916.json", METHOD_SHA)
    files, datasets = portable_inputs(inputs)
    workflow_paths = subprocess.check_output(["git", "-C", str(root), "ls-files", "benchmark_tools/*.py"], text=True).splitlines()
    workflow_paths = [p for p in workflow_paths if len(PurePosixPath(p).parts) == 2]
    if not workflow_paths:
        raise ValueError("No tracked workflow sources")
    subprocess.run(["git", "-C", str(root), "diff", "--exit-code", "HEAD", "--", *workflow_paths], check=True)
    revision = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    report = {"status": "preparing_transfer", "source": record(__file__), "input_manifest_sha256": INPUT_SHA,
              "baseline_manifest_sha256": METHOD_SHA, "workflow_revision": revision,
              "core_revision": baseline["core_commit"], "inputs": files, "datasets": datasets,
              "planned_runs": inputs["planned_runs"], "resource_plan": inputs["resource_plan"],
              "reference_package_inventories": baseline["environments"], "workflow_files": [],
              "destination_validated": False, "inference_started": False,
              "requirements": ["Destination host/access, working directory and scheduler must be supplied separately.",
                  "Use the frozen core git checkout; no machine-specific native binaries are transferred.",
                  "Record destination-specific native builds and profile/output checks; the existing AVX2 builder is not portable to ARM.",
                  "Install and verify OrthoFinder3.1.5 and the frozen external tool versions; do not reuse local executable paths.",
                  "Regenerate absolute commands and runtime manifests at destination without changing settings, run order or inputs.",
                  "Validate actual CPU/memory limits, dedicated workload, output semantics and instrumentation overhead before timing.",
                  "Transfer preparation is neither a dependency lock nor validated destination execution."]}
    output.mkdir(parents=True)
    try:
        for item, original in zip(files, inputs["ordered_proteomes"]):
            copy_verified(Path(original["input"]["path"]), output / item["path"], item)
        for relative in workflow_paths:
            source = root / relative
            target = "workflow/" + relative
            item = record(source)
            copy_verified(source, output / target, item)
            report["workflow_files"].append({**item, "path": target})
        report["status"] = "scaling_transfer_prepared_not_sent"
    except Exception as error:
        report.update(status="transfer_preparation_failed", error=str(error))
        raise
    finally:
        (output / "manifest.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.output.resolve())
