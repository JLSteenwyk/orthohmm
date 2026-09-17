"""Run four prespecified candidate variants through frozen inferred phylogeny."""

import argparse
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.recover_factorial_postflight import PREPARED_HASH, ENVIRONMENT_HASH
from benchmark_tools.run_orthobench_factorial_cell import select_cell, verify_prepared
from benchmark_tools.run_simulation_methods import read_frozen, verify_environment, execution_environment, execute
from benchmark_tools.validate_factorial_native import validate

ADMISSION_SHA = "38a4cc8c16e8e45800c88af618f9b5a9b3d77ab6b3d0e5b79fa555b88e849514"
VARIANTS = ("norm_low", "norm_high", "margin_low", "margin_high")
CPM_VARIANTS = ("cpm_low", "cpm_high")
CPM_ADMISSION_SHA = "5acd1c56fe72e6267913a1170efca7f5f189960785422b7194f6a5fcc545a5fd"


def select_variant_arm(admission, index, cpm=False):
    labels = CPM_VARIANTS if cpm else VARIANTS
    status = "cpm_candidates_verified_unscored" if cpm else "candidate_neighborhood_preparation_verified_unscored"
    if (type(index) is not int or not 0 <= index < len(labels) or admission["status"] != status
            or admission["accuracy_evaluated"] is not False
            or [row["label"] for row in admission["arms"]] != ["control", *labels]):
        raise ValueError("Invalid or unadmitted candidate panel")
    arm = admission["arms"][index + 1]
    if cpm:
        arm = {**arm, "partition": arm["candidate_partition"], "constraints": arm["membership_constraints"]}
    return arm


def variant_cell(original, arm, output):
    if original["label"] != "p1_c1_r1" or arm["label"] not in (*VARIANTS, *CPM_VARIANTS):
        raise ValueError("Unexpected baseline or candidate variant")
    argv = list(original["argv"])
    flags = ("--candidate-clusters", "--membership-constraints", "--output-directory", "--json", "--species-tree-mode")
    if any(argv.count(flag) != 1 for flag in flags):
        raise ValueError("Missing or repeated baseline argument")
    if (argv[argv.index("--species-tree-mode") + 1] != "infer"
            or "--species-tree" in argv or "--checkpoint-source" in argv):
        raise ValueError("Expected independent inferred-tree baseline")
    checkpoint = argv[argv.index("--output-directory") + 1]
    for flag, value in (("--candidate-clusters", arm["partition"]["path"]),
                        ("--membership-constraints", arm["constraints"]["path"]),
                        ("--output-directory", str(output / "output")),
                        ("--json", str(output / "metrics.json"))):
        argv[argv.index(flag) + 1] = value
    argv.extend(["--checkpoint-source", checkpoint])
    return {**original, "label": "candidate_" + arm["label"], "argv": argv,
        "candidate_partition": arm["partition"]["path"], "checkpoint_source": checkpoint,
        "prediction": str(output / "output/orthohmm_phylogeny/orthohmm_root_hogs.tsv")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--index", required=True, type=int, choices=range(4))
    parser.add_argument("--cpm", action="store_true", help="Use the separately admitted CPM-specific candidate panel")
    args = parser.parse_args()
    root = args.root.resolve()
    if Path.cwd().resolve() != root:
        raise ValueError("Start from original repository for environment verification")
    results = root / "benchmark_tools/results"
    admission_path = results / "ob_candidate_neighborhood_verified_20260916.json"
    if args.cpm:
        admission_path = results / "ob_cpm_candidates_verified_20260916.json"
    admission = read_frozen(admission_path, CPM_ADMISSION_SHA if args.cpm else ADMISSION_SHA)
    arm = select_variant_arm(admission, args.index, args.cpm)
    for item in admission["provenance_checked"]:
        check(item)
    prepared_path = results / "orthobench_factorial_prepared_20260916.json"
    environment_path = results / "publication_variable_native_methods_20260916.json"
    prepared = read_frozen(prepared_path, PREPARED_HASH)
    environment = read_frozen(environment_path, ENVIRONMENT_HASH)
    original, _, launcher = select_cell(prepared, 3)
    verify_prepared(prepared, original, launcher, 21161)
    verify_environment(environment)
    baseline = validate(root, 3)
    output = root / "benchmarks/results/ob_candidate_neighborhood_phylogeny_v1" / arm["label"]
    if args.cpm:
        output = root / "benchmarks/results/ob_cpm_phylogeny_v1" / arm["label"]
    if output.exists():
        raise FileExistsError(output)
    cell = variant_cell(original, arm, output)
    env, resolved = execution_environment(environment)
    env["PYTHONPATH"] = str(launcher)
    env.update(prepared["environment_overrides"])
    provenance = {"source": record(__file__), "admission": record(admission_path),
        "prepared": record(prepared_path), "environment": record(environment_path),
        "helpers": [record(Path(__file__).with_name(name)) for name in
            ("run_orthobench_factorial_cell.py", "run_simulation_methods.py", "validate_factorial_native.py")],
        "baseline_native_admission": baseline, "arm": arm, "cell": cell,
        "resolved_tools": resolved, "cwd": str(launcher),
        "job_id": os.environ.get("SLURM_JOB_ID"), "array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"),
        "array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
        "scope": "Unscored fixed candidate-threshold variant; inferred species tree; validated raw-tree reuse; incremental shared-node timing"}
    if args.cpm:
        provenance["scope"] = "Unscored prespecified CPM variant with its own HMM seeds/candidates; inferred species tree; validated raw-tree reuse; incremental shared-node timing"
    output.mkdir(parents=True)
    (output / "preflight.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    method = {"argv": cell["argv"], "output": str(output / "output"), "metrics": str(output / "metrics.json")}
    inputs = {"status": "ready", "inputs": [{**r, "absolute_path": r["path"]} for r in prepared["fasta_inputs"]]}
    try:
        os.chdir(launcher)
        execution = execute({"label": cell["label"], "methods": {cell["label"]: method}},
            [cell["label"]], env, output / "execution", inputs, provenance)
    finally:
        os.chdir(root)
    verify_prepared(prepared, original, launcher, 21161)
    verify_environment(environment)
    for item in admission["provenance_checked"]:
        check(item)
    if validate(root, 3) != baseline:
        raise ValueError("Baseline checkpoint source changed")
    if execution.get("failed_methods"):
        raise SystemExit("Candidate variant failed; evidence preserved without retry")
    (output / "postflight.json").write_text(json.dumps({"status": "complete_pending_native_validation",
        "accuracy_evaluated": False, "baseline_source_unchanged": True, "cell": cell},
        indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
