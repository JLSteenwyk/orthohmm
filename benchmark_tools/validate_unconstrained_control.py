"""Validate completed OrthoBench unconstrained diagnostic without reading labels."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.recover_factorial_postflight import PREPARED_HASH, ENVIRONMENT_HASH
from benchmark_tools.run_orthobench_factorial_cell import select_cell, unconstrained_cell, verify_prepared
from benchmark_tools.run_simulation_methods import read_frozen, verify_environment, execution_environment
from benchmark_tools.validate_factorial_native import validate_native_cell
from benchmark_tools.validate_simulation_outputs import verify_process
from benchmark_tools.verify_ygob_validation import require_completed_job

EXECUTOR_COMMIT = "d52add73035ebe3472e3e1cb95c38a9421b09dee"


def check_status(status, cell):
    if (status.get("status") != "finished_pending_native_validation" or status.get("failed_methods") != []
            or status.get("accuracy_evaluated") is not False or status.get("native_outputs_validated") is not False):
        raise ValueError("Diagnostic inference has not completed successfully and unscored")
    if status["provenance"].get("slurm_job_id") != "21290" or status["provenance"].get("slurm_array_task_id") is not None:
        raise ValueError("Wrong diagnostic scheduler identity")
    if status["dataset"] != cell["label"] or set(status["methods"]) != {cell["label"]}:
        raise ValueError("Wrong diagnostic method inventory")
    if status["provenance"].get("parent_cell") != "p1_c1_r1" or status["provenance"].get("omitted_membership_constraints") != cell["omitted_membership_constraints"]:
        raise ValueError("Wrong diagnostic treatment provenance")


def validate(root):
    root = root.resolve()
    if Path.cwd().resolve() != root:
        raise ValueError("Run from the original repository verification directory")
    accounting = subprocess.check_output(["sacct", "-j", "21290", "--parsable2",
                                         "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 21290)
    prepared_path = root / "benchmark_tools/results/orthobench_factorial_prepared_20260916.json"
    environment_path = root / "benchmark_tools/results/publication_variable_native_methods_20260916.json"
    prepared, environment = read_frozen(prepared_path, PREPARED_HASH), read_frozen(environment_path, ENVIRONMENT_HASH)
    original, output, launcher = select_cell(prepared, 3)
    cell = unconstrained_cell(original, output)
    executor = root / "benchmarks/work/publication_ob_unconstrained_executor_v2"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR_COMMIT:
        raise ValueError("Diagnostic executor revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    status_path = output / "execution" / cell["label"] / "status.json"
    status_record = file_provenance(status_path)
    status = json.loads(status_path.read_text())
    check_status(status, cell)
    provenance = status["provenance"]
    for key, path in {"prepared_manifest": prepared_path, "environment_manifest": environment_path,
                      "executor": executor / "benchmark_tools/run_orthobench_factorial_cell.py",
                      "execution_helper": executor / "benchmark_tools/run_simulation_methods.py",
                      "launcher": Path(prepared["launcher"]["path"])}.items():
        if provenance[key] != file_provenance(path):
            raise ValueError("Diagnostic source provenance changed: " + key)
    if provenance["cwd"] != str(launcher):
        raise ValueError("Wrong inference working directory")
    expected_inputs = {"status": "ready", "inputs": [{**r, "absolute_path": r["path"]} for r in prepared["fasta_inputs"]]}
    if status["verified_inputs"] != expected_inputs:
        raise ValueError("Diagnostic input inventory differs")
    verify_prepared(prepared, cell, launcher, 21161)
    verify_environment(environment)
    _, resolved = execution_environment(environment)
    if provenance["resolved_tools"] != resolved:
        raise ValueError("Diagnostic tool resolution changed")
    argv = cell["argv"]
    config = {"argv": argv, "output": argv[argv.index("--output-directory") + 1],
              "metrics": argv[argv.index("--json") + 1]}
    artifacts = verify_process(config, status["methods"][cell["label"]])
    integrity = {"scheduler": scheduler, "execution": status_record, "verified_artifacts": len(artifacts),
                 "executor_commit": EXECUTOR_COMMIT, "accuracy_evaluated": False}
    result = validate_native_cell(prepared, environment, cell, output, launcher, integrity)
    verify_file(status_path, status_record)
    result["verifier"] = file_provenance(Path(__file__))
    result["diagnostic"] = "Unconstrained satellite replay, compared with p1_c1_r1; not a ninth factorial cell"
    return result


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
