"""Execute the frozen DIAMOND control, preserving phase failures and provenance."""

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.prepare_sequence_search_control import search_command
from benchmark_tools.run_simulation_methods import read_frozen

MANIFEST_SHA = "9edb1bb8232e114e12a6199412dd45a31d740c8a1dcba36fcf05d2c58b23095b"


def verify_plan(report):
    if report["status"] != "prepared_not_searched" or report["accuracy_evaluated"] is not False:
        raise ValueError("Expected unscored frozen search preparation")
    records = [report[k] for k in ("source", "frozen_factorial", "queries", "gene_metadata", "diamond")]
    for record in [*records, *report["inputs"]]:
        verify_file(Path(record["path"]), record)
    binary = report["diamond"]["path"]
    if subprocess.check_output([binary, "version"], text=True).strip() != report["diamond_version"]:
        raise ValueError("DIAMOND version changed")
    if report["diamond_version"] != "diamond version 2.1.11" or len(report["searches"]) != 12:
        raise ValueError("Unexpected frozen search panel")
    root = Path(report["queries"]["path"]).parent
    for index, (cell, source) in enumerate(zip(report["searches"], report["inputs"], strict=True)):
        directory = root / f"target_{index:02d}"
        db, hits = directory / "target", directory / "hits.tsv"
        expected = {"index": index, "target_fasta": source, "output": str(hits),
                    "makedb": [binary, "makedb", "--in", source["path"], "--db", str(db), "--threads", "32"],
                    "search": search_command(binary, report["queries"]["path"], db, hits)}
        if cell != expected:
            raise ValueError("Search command or target order differs from preparation")
    return root


def run_phase(command, directory, name, env):
    stdout, timing = directory / (name + ".log"), directory / (name + ".time.log")
    if stdout.exists() or timing.exists():
        raise FileExistsError("Phase evidence already exists")
    started = time.monotonic()
    with stdout.open("x") as handle:
        completed = subprocess.run(["/usr/bin/time", "-v", "-o", str(timing), *command],
                                   cwd=directory, env=env, stdout=handle, stderr=subprocess.STDOUT)
    return {"argv": command, "exit_code": completed.returncode, "wall_s": time.monotonic() - started,
            "log": file_provenance(stdout), "gnu_time": file_provenance(timing)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    report = read_frozen(args.manifest, MANIFEST_SHA)
    root = verify_plan(report)
    status_path = root / "execution.json"
    if status_path.exists() or any(any(Path(cell["output"]).parent.iterdir()) for cell in report["searches"]):
        raise FileExistsError("Search output/evidence directories are not pristine")
    if args.check_only:
        print("Verified 12 frozen searches; no execution")
        return
    env = os.environ.copy()
    overrides = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    env.update(overrides)
    status = {"schema_version": 1, "status": "running", "accuracy_evaluated": False,
              "started_at": datetime.now(timezone.utc).isoformat(), "job_id": os.environ.get("SLURM_JOB_ID"),
              "manifest": file_provenance(args.manifest), "executor": file_provenance(Path(__file__)),
              "environment_overrides": overrides, "targets": [],
              "runtime_kind": "shared-machine database/search phases; GNU time RSS not simultaneous process-tree RSS"}
    with status_path.open("x") as handle:
        json.dump(status, handle, indent=2, sort_keys=True)
    try:
        for cell in report["searches"]:
            target = {"index": cell["index"], "status": "running", "phases": {}}
            status["targets"].append(target)
            directory = Path(cell["output"]).parent
            for phase in ("makedb", "search"):
                target["phases"][phase] = run_phase(cell[phase], directory, phase, env)
                status_path.write_text(json.dumps(status, indent=2, sort_keys=True) + "\n")
                if target["phases"][phase]["exit_code"] != 0:
                    target["status"] = "failed"
                    raise RuntimeError(f"Target {cell['index']} {phase} failed")
            target["hits"] = file_provenance(Path(cell["output"]))
            target["database"] = file_provenance(directory / "target.dmnd")
            target["status"] = "complete_pending_numeric_validation"
            status_path.write_text(json.dumps(status, indent=2, sort_keys=True) + "\n")
        verify_file(args.manifest, status["manifest"])
        verify_plan(report)
        status["status"] = "complete_pending_numeric_validation"
    except Exception as error:
        status.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        status["finished_at"] = datetime.now(timezone.utc).isoformat()
        status_path.write_text(json.dumps(status, indent=2, sort_keys=True) + "\n")
    print("All 12 searches completed; numeric conversion and scoring remain separate")


if __name__ == "__main__":
    main()
