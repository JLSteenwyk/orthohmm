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


def control_cell(original, tree, output):
    if original["label"] != "p1_c1_r1" or tree["label"] != "supplied_control" or tree["rooted_rf_clade_distance"] != 0:
        raise ValueError("Only the unchanged supplied-tree control is admitted")
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
    return {**original, "label": "p1_c1_r1_supplied_control", "argv": args,
            "checkpoint_source": source, "species_tree": tree["tree"],
            "prediction": str(output / "output/orthohmm_phylogeny/orthohmm_root_hogs.tsv")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
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
    candidates = [t for t in trees["variants"] if t["label"] == "supplied_control"]
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
    if output.exists():
        raise FileExistsError(output)
    cell = control_cell(original, tree, output)
    env, resolved = execution_environment(environment)
    env["PYTHONPATH"] = str(launcher)
    env.update(prepared["environment_overrides"])
    provenance = {"source": file_provenance(Path(__file__)), "tree_manifest": file_provenance(tree_path),
                  "prepared_manifest": file_provenance(prepared_path), "environment_manifest": file_provenance(environment_path),
                  "resolved_tools": resolved, "baseline_native_admission": baseline, "cell": cell,
                  "job_id": os.environ.get("SLURM_JOB_ID"), "cwd": str(launcher),
                  "scope": "Incremental supplied-tree control with reused raw gene trees; no accuracy evaluation"}
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
    (output / "postflight.json").write_text(json.dumps({"status": "complete_pending_native_equivalence",
        "accuracy_evaluated": False, "baseline_source_unchanged": True, "cell": cell}, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
