"""Admit the complete prespecified tree-perturbation panel without scoring."""

import argparse
import csv
import io
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import Phylo
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.prepare_species_tree_robustness import topology
from benchmark_tools.recover_factorial_postflight import PREPARED_HASH, ENVIRONMENT_HASH
from benchmark_tools.run_orthobench_factorial_cell import select_cell, verify_prepared
from benchmark_tools.run_simulation_methods import read_frozen, verify_environment
from benchmark_tools.run_species_tree_control import (
    TREE_MANIFEST_SHA, CONTROL_VALIDATION_SHA, PERTURBATIONS, perturbation_cell,
)
from benchmark_tools.validate_factorial_native import validate_native_cell
from benchmark_tools.validate_simulation_outputs import verify_process
from benchmark_tools.validate_species_tree_control import validate as validate_control


def completed_tasks(accounting):
    rows = list(csv.DictReader(io.StringIO(accounting), delimiter="|"))
    expected = {f"21299_{i}" for i in range(6)}
    selected = [row for row in rows if row["JobID"] in expected]
    if len(selected) != 6 or {row["JobID"] for row in selected} != expected:
        raise ValueError("Missing or duplicate perturbation tasks")
    if any(row["State"] != "COMPLETED" or row["ExitCode"] != "0:0" for row in selected):
        raise ValueError("All six perturbations must complete successfully before admission")
    if len({row["JobIDRaw"] for row in selected}) != 6:
        raise ValueError("Duplicate raw scheduler identities")
    return {int(row["JobID"].split("_")[1]): row for row in selected}


def check_status(status, preflight, postflight, cell, tree, index, scheduler, launcher):
    if (status["status"] != "finished_pending_native_validation" or status["failed_methods"] != []
            or status["accuracy_evaluated"] is not False or status["native_outputs_validated"] is not False):
        raise ValueError("Perturbation is not successfully completed and unscored")
    if status["dataset"] != cell["label"] or set(status["methods"]) != {cell["label"]}:
        raise ValueError("Perturbation method inventory differs")
    if status["provenance"] != preflight:
        raise ValueError("Execution provenance differs from preflight")
    expected = {"job_id": scheduler["JobIDRaw"], "array_job_id": "21299", "array_task_id": str(index),
                "perturbation_index": index, "tree_variant": tree, "cell": cell, "cwd": str(launcher)}
    if any(preflight[key] != value for key, value in expected.items()):
        raise ValueError("Perturbation identity, command or working directory differs")
    if preflight["control_admission"]["status"] != "equivalent":
        raise ValueError("Supplied-control admission failed")
    if postflight != {"status": "complete_pending_native_validation", "accuracy_evaluated": False,
                      "baseline_source_unchanged": True, "cell": cell}:
        raise ValueError("Perturbation postflight did not pass")


def check_topology(source, supplied, observed, expected_distance):
    taxa = lambda tree: [leaf.name for leaf in tree.get_terminals()]
    names = taxa(source)
    if len(names) != len(set(names)) or any(name is None for name in names):
        raise ValueError("Invalid source taxa")
    if sorted(taxa(supplied)) != sorted(names) or sorted(taxa(observed)) != sorted(names):
        raise ValueError("Perturbation taxon inventory changed")
    distance = len(topology(source) ^ topology(supplied))
    if distance != expected_distance or topology(supplied) != topology(observed):
        raise ValueError("Perturbation topology or rooted RF distance differs")
    return distance


def validate(root):
    root = root.resolve()
    if Path.cwd().resolve() != root:
        raise ValueError("Run from the original repository verification directory")
    accounting = subprocess.check_output(["sacct", "-j", "21299", "--parsable2",
                                         "--format=JobID,JobIDRaw,State,ExitCode,Elapsed"], text=True)
    tasks = completed_tasks(accounting)
    results = root / "benchmark_tools/results"
    prepared_path = results / "orthobench_factorial_prepared_20260916.json"
    env_path = results / "publication_variable_native_methods_20260916.json"
    trees_path = results / "ob_species_tree_robustness_prepared_20260916.json"
    prepared = read_frozen(prepared_path, PREPARED_HASH)
    environment = read_frozen(env_path, ENVIRONMENT_HASH)
    trees = read_frozen(trees_path, TREE_MANIFEST_SHA)
    frozen_control = read_frozen(results / "ob_supplied_tree_control_validation_20260916.json", CONTROL_VALIDATION_SHA)
    control = validate_control(root)
    if control["status"] != "equivalent" or frozen_control["status"] != "equivalent":
        raise ValueError("Supplied-control equivalence is not established")
    baseline = control["baseline_validation"]
    if baseline["species_tree"] != trees["source_tree"]:
        raise ValueError("Prepared trees use a different source")
    original, _, launcher = select_cell(prepared, 3)
    executor = root / "benchmarks/work/publication_ob_tree_perturbations_v1"
    revision = subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip()
    expected_revision = subprocess.check_output(["git", "-C", str(root), "rev-parse", "cfe0a09^{commit}"], text=True).strip()
    if revision != expected_revision:
        raise ValueError("Perturbation executor revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    verify_prepared(prepared, original, launcher, 21161)
    verify_environment(environment)
    source_tree = Phylo.read(trees["source_tree"]["path"], "newick")
    reports, tracked = {}, []
    for index, label in enumerate(PERTURBATIONS):
        matches = [tree for tree in trees["variants"] if tree["label"] == label]
        if len(matches) != 1:
            raise ValueError("Missing or duplicate prepared tree")
        tree = matches[0]
        verify_file(Path(tree["tree"]["path"]), tree["tree"])
        output = root / "benchmarks/results/ob_species_tree_perturbations_v1" / label
        cell = perturbation_cell(original, tree, output, index)
        records = [file_provenance(output / name) for name in ("preflight.json", "postflight.json", "execution/status.json")]
        preflight, postflight, status = [json.loads(Path(record["path"]).read_text()) for record in records]
        check_status(status, preflight, postflight, cell, tree, index, tasks[index], launcher)
        for key, path in {"source": executor / "benchmark_tools/run_species_tree_control.py",
                          "tree_manifest": trees_path, "prepared_manifest": prepared_path, "environment_manifest": env_path}.items():
            if preflight[key] != file_provenance(path):
                raise ValueError("Input or executor provenance differs: " + key)
        for key in ("native_manifest", "native_metrics", "species_tree", "partition"):
            if baseline[key] != preflight["baseline_native_admission"][key]:
                raise ValueError("Original native checkpoint source changed")
        method = {"argv": cell["argv"], "output": str(output / "output"), "metrics": str(output / "metrics.json")}
        artifacts = verify_process(method, status["methods"][cell["label"]])
        native = validate_native_cell(prepared, environment, cell, output, launcher,
                                      {"scheduler": tasks[index], "verified_artifacts": len(artifacts)})
        distance = check_topology(source_tree, Phylo.read(tree["tree"]["path"], "newick"),
                                  Phylo.read(native["species_tree"]["path"], "newick"), tree["rooted_rf_clade_distance"])
        reports[label] = {"native_validation": native, "execution": records, "cell": cell,
                          "rooted_rf_clade_distance": distance, "scheduler": tasks[index]}
        tracked.extend([*records, tree["tree"]])
    for record in tracked:
        verify_file(Path(record["path"]), record)
    verify_prepared(prepared, original, launcher, 21161)
    verify_environment(environment)
    return {"status": "native_panel_validated_unscored", "accuracy_evaluated": False, "variants": reports,
            "control_validation": control, "verifier": file_provenance(Path(__file__)),
            "limitations": "Incremental supplied-tree perturbations with reused raw gene trees; not end-to-end timing or posterior uncertainty."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = validate(args.root)
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
