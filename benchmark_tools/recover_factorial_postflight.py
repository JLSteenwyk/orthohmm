"""Audit preserved array 21248 outputs after the known cwd postflight defect.

This does not repair scheduler history, validate native semantics, or score.
"""

import argparse
import csv
import io
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.run_orthobench_factorial_cell import select_cell, verify_prepared
from benchmark_tools.run_simulation_methods import read_frozen, verify_environment, execution_environment
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.validate_simulation_outputs import verify_process

COMMIT = "9a8630197401d97e4fcc131f4939433bd3b89aac"
PREPARED_HASH = "5c325f4d77865e0c7571fe4bb4df0be0959977a49e1d22169f192518700f9382"
ENVIRONMENT_HASH = "bf677728cf9312edd9755d9ac1ceefbce3abb8c3848d43413631dc00ca00eb2f"


def require_known_failure(accounting, index, status, log):
    rows = list(csv.DictReader(io.StringIO(accounting), delimiter="|"))
    matches = [r for r in rows if r["JobID"] == f"21248_{index}"]
    if len(matches) != 1 or matches[0]["State"] != "FAILED" or matches[0]["ExitCode"] != "1:0":
        raise ValueError("Expected unique terminal known-failure task")
    row = matches[0]
    if (status["provenance"]["slurm_job_id"] != row["JobIDRaw"] or
            status["provenance"]["slurm_array_task_id"] != str(index)):
        raise ValueError("Scheduler identity mismatch")
    if (status.get("status") != "finished_pending_native_validation" or
            status.get("failed_methods") != [] or status.get("accuracy_evaluated") is not False or
            status.get("native_outputs_validated") is not False):
        raise ValueError("Unexpected inference status")
    if (log.count("Traceback (most recent call last):") != 1 or
            not log.rstrip().endswith("ValueError: Package inventory changed: orthofinder") or
            'line 117, in main\n    verify_environment(environment)' not in log):
        raise ValueError("Not the known postflight traceback")
    return row


def audit(root, index):
    if Path.cwd().resolve() != root:
        raise ValueError("Run from the original repository verification directory")
    results = root / "benchmark_tools/results"
    prepared_path = results / "orthobench_factorial_prepared_20260916.json"
    environment_path = results / "publication_variable_native_methods_20260916.json"
    prepared = read_frozen(prepared_path, PREPARED_HASH)
    environment = read_frozen(environment_path, ENVIRONMENT_HASH)
    cell, output, launcher = select_cell(prepared, index)
    status_path = output / "execution" / cell["label"] / "status.json"
    status_record = file_provenance(status_path)
    status = json.loads(status_path.read_text())
    log_path = root / f"benchmarks/work/publication_ob_factorial_logs_v1/cell_21248_{index}.log"
    log_record = file_provenance(log_path)
    accounting = subprocess.check_output(["sacct", "-j", "21248", "--parsable2",
                                          "--format=JobID,JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_known_failure(accounting, index, status, log_path.read_text())
    executor = root / "benchmarks/work/publication_ob_factorial_execution_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != COMMIT:
        raise ValueError("Original executor revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    provenance = status["provenance"]
    expected = {"prepared_manifest": prepared_path, "environment_manifest": environment_path,
                "executor": executor / "benchmark_tools/run_orthobench_factorial_cell.py",
                "execution_helper": executor / "benchmark_tools/run_simulation_methods.py",
                "launcher": Path(prepared["launcher"]["path"])}
    for key, path in expected.items():
        if provenance[key] != file_provenance(path):
            raise ValueError(f"Execution provenance mismatch: {key}")
    if provenance["cwd"] != str(launcher) or status["dataset"] != cell["label"]:
        raise ValueError("Execution context mismatch")
    inputs = {"status": "ready", "inputs": [{**r, "absolute_path": r["path"]} for r in prepared["fasta_inputs"]]}
    if status["verified_inputs"] != inputs or set(status["methods"]) != {cell["label"]}:
        raise ValueError("Input or method inventory mismatch")
    verify_prepared(prepared, cell, launcher, 21161)
    verify_environment(environment)
    _, resolved = execution_environment(environment)
    if provenance["resolved_tools"] != resolved:
        raise ValueError("Executable resolution changed")
    argv = cell["argv"]
    config = {"argv": argv, "output": argv[argv.index("--output-directory") + 1],
              "metrics": argv[argv.index("--json") + 1]}
    artifacts = verify_process(config, status["methods"][cell["label"]])
    verify_file(status_path, status_record)
    verify_file(log_path, log_record)
    return {"schema_version": 1, "cell": cell["label"], "status": "postflight_integrity_verified",
            "scheduler": scheduler, "accounting_raw": accounting, "original_status": status_record,
            "original_batch_log": log_record, "executor_commit": COMMIT,
            "verification_cwd": str(root), "verifier": file_provenance(Path(__file__)),
            "prepared_manifest": file_provenance(prepared_path),
            "environment_manifest": file_provenance(environment_path),
            "verified_artifact_count": len(artifacts), "native_outputs_validated": False,
            "accuracy_evaluated": False, "scoring_admitted": False,
            "limitation": "Original scheduler failure retained; native semantics and conversion gates remain required."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--index", type=int, choices=range(4), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = audit(args.root.resolve(), args.index)
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"{report['cell']}: {report['verified_artifact_count']} artifacts verified; not admitted for scoring")


if __name__ == "__main__":
    main()
