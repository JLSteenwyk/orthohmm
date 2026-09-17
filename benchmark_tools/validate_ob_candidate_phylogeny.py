"""Admit native outputs for fixed threshold or CPM phylogeny variants."""

import argparse
from copy import deepcopy
import csv
import io
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.recover_factorial_postflight import PREPARED_HASH, ENVIRONMENT_HASH
from benchmark_tools.run_ob_candidate_neighborhood import (
    ADMISSION_SHA, CPM_ADMISSION_SHA, VARIANTS, CPM_VARIANTS, select_variant_arm, variant_cell,
)
from benchmark_tools.run_orthobench_factorial_cell import select_cell, verify_prepared
from benchmark_tools.run_simulation_methods import read_frozen, verify_environment
from benchmark_tools.validate_factorial_native import validate as validate_baseline, validate_native_cell
from benchmark_tools.validate_simulation_outputs import verify_process


def completed_task(accounting, index, cpm=False):
    rows = list(csv.DictReader(io.StringIO(accounting), delimiter="|"))
    job = "21324" if cpm else "21316"
    selected = [row for row in rows if row["JobID"] == f"{job}_{index}"]
    if len(selected) != 1 or selected[0]["State"] != "COMPLETED" or selected[0]["ExitCode"] != "0:0":
        raise ValueError("Variant has not completed successfully")
    return selected[0]


def adapted_prepared(prepared, arm):
    if arm["label"] not in (*VARIANTS, *CPM_VARIANTS):
        raise ValueError("Unknown candidate variant")
    result = deepcopy(prepared)
    target = result["candidate_arms"]["p1_c1"]
    target["candidate_partition"] = arm["partition"]
    target["membership_constraints"] = arm["constraints"]
    if arm["label"] in CPM_VARIANTS:
        target["seed_partition"] = arm["seed_partition"]
    return result


def check_execution(status, preflight, postflight, cell, arm, scheduler, index, launcher, cpm=False):
    if (status["status"] != "finished_pending_native_validation" or status["failed_methods"] != []
            or status["accuracy_evaluated"] is not False or status["native_outputs_validated"] is not False
            or status["dataset"] != cell["label"] or set(status["methods"]) != {cell["label"]}
            or status["provenance"] != preflight):
        raise ValueError("Execution is not successful, unscored and provenance-consistent")
    expected = {"job_id": scheduler["JobIDRaw"], "array_job_id": "21324" if cpm else "21316", "array_task_id": str(index),
                "cell": cell, "arm": arm, "cwd": str(launcher)}
    if any(preflight[key] != value for key, value in expected.items()):
        raise ValueError("Variant identity or command differs")
    if postflight != {"status": "complete_pending_native_validation", "accuracy_evaluated": False,
                      "baseline_source_unchanged": True, "cell": cell}:
        raise ValueError("Variant postflight failed")


def check_baseline(current, prior):
    # Validator paths differ between the frozen executor and this admission.
    keys = ("cell", "status", "native_manifest", "native_metrics", "species_tree", "partition",
            "membership", "native_outputs_validated", "accuracy_evaluated")
    if any(current[key] != prior[key] for key in keys):
        raise ValueError("Baseline native source changed")


def validate(root, index, cpm=False):
    if type(index) is not int or index not in range(2 if cpm else 4):
        raise ValueError("Invalid variant index")
    if Path.cwd().resolve() != root:
        raise ValueError("Run from original environment verification directory")
    accounting = subprocess.check_output(["sacct", "-j", "21324" if cpm else "21316", "--parsable2",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = completed_task(accounting, index, cpm)
    results = root / "benchmark_tools/results"
    admission_path = results / "ob_candidate_neighborhood_verified_20260916.json"
    if cpm:
        admission_path = results / "ob_cpm_candidates_verified_20260916.json"
    admission = read_frozen(admission_path, CPM_ADMISSION_SHA if cpm else ADMISSION_SHA)
    arm = select_variant_arm(admission, index, cpm)
    for item in admission["provenance_checked"]:
        check(item)
    prepared_path = results / "orthobench_factorial_prepared_20260916.json"
    environment_path = results / "publication_variable_native_methods_20260916.json"
    prepared = read_frozen(prepared_path, PREPARED_HASH)
    environment = read_frozen(environment_path, ENVIRONMENT_HASH)
    original, _, launcher = select_cell(prepared, 3)
    verify_prepared(prepared, original, launcher, 21161)
    verify_environment(environment)
    baseline = validate_baseline(root, 3)
    executor = root / "benchmarks/work/publication_ob_candidate_phylogeny_v1"
    if cpm:
        executor = root / "benchmarks/work/publication_ob_cpm_phylogeny_v1"
    revision = subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip()
    expected = subprocess.check_output(["git", "-C", str(root), "rev-parse",
        "7e8c3e1^{commit}" if cpm else "e8e86d8^{commit}"], text=True).strip()
    if revision != expected:
        raise ValueError("Changed frozen executor revision")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    output = root / "benchmarks/results/ob_candidate_neighborhood_phylogeny_v1" / arm["label"]
    if cpm:
        output = root / "benchmarks/results/ob_cpm_phylogeny_v1" / arm["label"]
    cell = variant_cell(original, arm, output)
    records = [record(output / name) for name in ("preflight.json", "postflight.json", "execution/status.json")]
    preflight, postflight, status = [json.loads(Path(item["path"]).read_text()) for item in records]
    check_execution(status, preflight, postflight, cell, arm, scheduler, index, launcher, cpm)
    for key, path in {"source": executor / "benchmark_tools/run_ob_candidate_neighborhood.py",
                      "admission": admission_path, "prepared": prepared_path, "environment": environment_path}.items():
        if preflight[key] != record(path):
            raise ValueError("Input/executor provenance differs: " + key)
    for item in preflight["helpers"]:
        check(item)
    check_baseline(baseline, preflight["baseline_native_admission"])
    method = {"argv": cell["argv"], "output": str(output / "output"), "metrics": str(output / "metrics.json")}
    artifacts = verify_process(method, status["methods"][cell["label"]])
    native = validate_native_cell(adapted_prepared(prepared, arm), environment, cell, output, launcher,
        {"scheduler": scheduler, "verified_artifacts": len(artifacts), "execution": records})
    for item in [*records, *admission["provenance_checked"]]:
        check(item)
    verify_environment(environment)
    return {"status": "candidate_variant_native_validated_unscored", "accuracy_evaluated": False,
        "variant": arm["label"], "scheduler": scheduler, "native_validation": native,
        "cell": cell, "admission": record(admission_path), "source": record(__file__),
        "helpers": [record(Path(__file__).with_name(name)) for name in
            ("run_ob_candidate_neighborhood.py", "validate_factorial_native.py", "validate_factorial_partition.py",
             "validate_simulation_outputs.py")],
        "limitations": "Native group-output admission, not independent tree reconstruction or accuracy evaluation; incremental cached shared-node timing."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--index", required=True, type=int, choices=range(4))
    parser.add_argument("--cpm", action="store_true", help="Admit the separately frozen two-arm CPM panel")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = validate(args.root.resolve(), args.index, args.cpm)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
