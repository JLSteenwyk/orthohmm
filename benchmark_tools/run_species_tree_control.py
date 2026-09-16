"""Run the unchanged supplied-tree control before tree perturbation experiments."""

import argparse
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.recover_factorial_postflight import PREPARED_HASH, ENVIRONMENT_HASH
from benchmark_tools.run_orthobench_factorial_cell import select_cell, verify_prepared
from benchmark_tools.run_simulation_methods import read_frozen, verify_environment, execution_environment, execute
from benchmark_tools.validate_factorial_native import validate

TREE_MANIFEST_SHA = "4ff0566c2ecac9e5afe91fee722a29312e2f82a35a771d592aa407dca3e0e8f5"
CONTROL_VALIDATION_SHA = "f80b5e45ebba88c0855a387fe489748e01f0d2bf399a981deb56519f745d7fa1"
PERTURBATIONS = tuple(f"nni{level}_{index}" for level in (1, 2) for index in range(3))


def control_cell(original, tree, output):
    if original["label"] != "p1_c1_r1" or tree["label"] != "supplied_control" or tree["rooted_rf_clade_distance"] != 0:
        raise ValueError("Only the unchanged supplied-tree control is admitted")
    return supplied_cell(original, tree, output)


def perturbation_cell(original, tree, output, index):
    if isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < 6:
        raise ValueError("Unknown perturbation index")
    if original["label"] != "p1_c1_r1" or tree["label"] != PERTURBATIONS[index] or tree["rooted_rf_clade_distance"] != (2 if index < 3 else 4):
        raise ValueError("Tree does not match the fixed perturbation inventory")
    return supplied_cell(original, tree, output)


def supplied_cell(original, tree, output):
    args = list(original["argv"])
    if args.count("--species-tree-mode") != 1 or args[args.index("--species-tree-mode") + 1] != "infer":
        raise ValueError("Expected the frozen inferred-tree baseline")
    if "--species-tree" in args or "--checkpoint-source" in args or args.count("--membership-constraints") != 1:
        raise ValueError("Unexpected baseline tree, checkpoint or membership arguments")
    source = args[args.index("--output-directory") + 1]
    args[args.index("--species-tree-mode") + 1] = "supplied"
    args[args.index("--output-directory") + 1] = str(output / "output")
    args[args.index("--json") + 1] = str(output / "metrics.json")
    args.extend(["--species-tree", tree["tree"]["path"], "--checkpoint-source", source])
    return {**original, "label": "p1_c1_r1_" + tree["label"], "argv": args,
            "checkpoint_source": source, "species_tree": tree["tree"],
            "prediction": str(output / "output/orthohmm_phylogeny/orthohmm_root_hogs.tsv")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--perturbation-index", type=int, choices=range(6))
    args = parser.parse_args()
    root = args.root.resolve()
    if Path.cwd().resolve() != root:
        raise ValueError("Start from the original repository verification directory")
    results = root / "benchmark_tools/results"
    tree_path = results / "ob_species_tree_robustness_prepared_20260916.json"
    trees = read_frozen(tree_path, TREE_MANIFEST_SHA)
    if trees["status"] != "trees_prepared_not_reconciled" or trees["accuracy_evaluated"] is not False:
        raise ValueError("Unexpected tree preparation status")
    baseline = validate(root, 3)
    if baseline["species_tree"] != trees["source_tree"]:
        raise ValueError("Prepared trees use a different source tree")
    perturbation = args.perturbation_index
    control_admission = None
    if perturbation is not None:
        from benchmark_tools.validate_species_tree_control import validate as validate_control
        frozen_admission = read_frozen(results / "ob_supplied_tree_control_validation_20260916.json", CONTROL_VALIDATION_SHA)
        control_admission = validate_control(root)
        if frozen_admission["status"] != "equivalent" or control_admission["status"] != "equivalent":
            raise ValueError("Supplied-tree baseline equivalence is not established")
    label = "supplied_control" if perturbation is None else PERTURBATIONS[perturbation]
    candidates = [t for t in trees["variants"] if t["label"] == label]
    if len(candidates) != 1:
        raise ValueError("Missing or duplicate control tree")
    tree = candidates[0]
    verify_file(Path(tree["tree"]["path"]), tree["tree"])
    prepared_path = results / "orthobench_factorial_prepared_20260916.json"
    environment_path = results / "publication_variable_native_methods_20260916.json"
    prepared = read_frozen(prepared_path, PREPARED_HASH)
    environment = read_frozen(environment_path, ENVIRONMENT_HASH)
    original, _, launcher = select_cell(prepared, 3)
    verify_prepared(prepared, original, launcher, 21161)
    verify_environment(environment)
    output = root / "benchmarks/results/ob_supplied_tree_control_v1"
    if perturbation is not None:
        output = root / "benchmarks/results/ob_species_tree_perturbations_v1" / label
    if output.exists():
        raise FileExistsError(output)
    cell = control_cell(original, tree, output) if perturbation is None else perturbation_cell(original, tree, output, perturbation)
    env, resolved = execution_environment(environment)
    env["PYTHONPATH"] = str(launcher)
    env.update(prepared["environment_overrides"])
    provenance = {"source": file_provenance(Path(__file__)), "tree_manifest": file_provenance(tree_path),
                  "prepared_manifest": file_provenance(prepared_path), "environment_manifest": file_provenance(environment_path),
                  "resolved_tools": resolved, "baseline_native_admission": baseline, "cell": cell,
                  "job_id": os.environ.get("SLURM_JOB_ID"), "cwd": str(launcher),
                  "scope": "Incremental supplied-tree control with reused raw gene trees; no accuracy evaluation"}
    if perturbation is not None:
        provenance.update(perturbation_index=perturbation, tree_variant=tree,
                          control_admission=control_admission,
                          array_job_id=os.environ.get("SLURM_ARRAY_JOB_ID"),
                          array_task_id=os.environ.get("SLURM_ARRAY_TASK_ID"),
                          scope="Prespecified rooted topology perturbation; incremental cached reconciliation; unscored")
    output.mkdir(parents=True)
    (output / "preflight.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    argv = cell["argv"]
    method = {"argv": argv, "output": str(output / "output"), "metrics": str(output / "metrics.json")}
    inputs = {"status": "ready", "inputs": [{**r, "absolute_path": r["path"]} for r in prepared["fasta_inputs"]]}
    try:
        os.chdir(launcher)
        result = execute({"label": cell["label"], "methods": {cell["label"]: method}},
                         [cell["label"]], env, output / "execution", inputs, provenance)
    finally:
        os.chdir(root)
    verify_prepared(prepared, original, launcher, 21161)
    verify_environment(environment)
    verify_file(Path(tree["tree"]["path"]), tree["tree"])
    after = validate(root, 3)
    if after != baseline:
        raise ValueError("Checkpoint source changed during copied replay")
    if result.get("failed_methods"):
        raise SystemExit("Supplied-tree replay failed; evidence preserved")
    (output / "postflight.json").write_text(json.dumps({"status": "complete_pending_native_equivalence" if perturbation is None else "complete_pending_native_validation",
        "accuracy_evaluated": False, "baseline_source_unchanged": True, "cell": cell}, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
