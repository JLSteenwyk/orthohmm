"""Execute the frozen 78-target corrected QfO sequence-search diagnostic."""

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_qfo_sequence_search_control import search_plan, validate_inputs, STAGING_SHA, INVENTORY_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_sequence_search_control import run_phase
from benchmark_tools.run_simulation_methods import read_frozen

MANIFEST_SHA = "fcef19d9c745c5806197aa219dec1c89bd5d90082a35f24dadafa3422f461436"


def verify_plan(manifest_path):
    report = read_frozen(manifest_path, MANIFEST_SHA)
    root = manifest_path.parent
    if (report["status"] != "corrected_qfo_search_control_prepared_unrun"
            or report["accuracy_evaluated"] is not False or report["execution_authorized"] is not False
            or report["genes"] != 984137 or report["proteomes"] != 78):
        raise ValueError("Require unrun corrected search preparation")
    for item in [report["source"], report["queries"], report["gene_metadata"], report["diamond"],
                 *report["checked_records"], *report["helpers"]]:
        check(item)
    stage_record, inventory_record = report["checked_records"][:2]
    stage = read_frozen(Path(stage_record["path"]), STAGING_SHA)
    inventory = read_frozen(Path(inventory_record["path"]), INVENTORY_SHA)
    if validate_inputs(stage, inventory) != report["inputs"]:
        raise ValueError("Corrected input order differs")
    if (report["queries"]["path"] != str(root / "queries.fasta")
            or report["gene_metadata"]["path"] != str(root / "gene_metadata.json")):
        raise ValueError("Prepared files are outside the search root")
    binary = Path(report["diamond"]["path"])
    version = subprocess.check_output([str(binary), "version"], text=True).strip()
    if version != report["diamond_version"] or version != "diamond version 2.1.11":
        raise ValueError("Search binary version changed")
    if report["searches"] != search_plan(report["inputs"], binary, root / "queries.fasta", root):
        raise ValueError("Search commands differ from frozen settings")
    return report


def run(manifest_path, check_only=False):
    report = verify_plan(manifest_path)
    status_path = manifest_path.parent / "execution.json"
    if status_path.exists() or status_path.is_symlink() or any(
            any(Path(cell["output"]).parent.iterdir()) for cell in report["searches"]):
        raise FileExistsError("Require pristine target and execution paths")
    if check_only:
        return {"status": "corrected_search_preflight_verified", "targets": len(report["searches"]),
                "accuracy_evaluated": False}
    if (not os.environ.get("SLURM_JOB_ID", "").isdigit()
            or os.environ.get("SLURM_CPUS_PER_TASK") != "32" or os.environ.get("SLURM_JOB_NODELIST") != "bizon"):
        raise ValueError("Require scheduled 32-CPU workstation execution")
    overrides = {"PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    env = {**os.environ, **overrides}
    status = {"status": "running", "job_id": os.environ["SLURM_JOB_ID"], "node": "bizon",
        "started_at": datetime.now(timezone.utc).isoformat(), "manifest": record(manifest_path),
        "executor": record(__file__), "phase_helper": record(Path(__file__).with_name("run_sequence_search_control.py")),
        "time_binary": record(Path("/usr/bin/time")), "environment_overrides": overrides,
        "targets": [], "accuracy_evaluated": False, "numeric_validated": False, "publication_ready": False,
        "runtime_kind": "Shared-host database/search phases; GNU-time RSS is not simultaneous process-tree RSS"}
    with status_path.open("x") as stream:
        json.dump(status, stream, indent=2, sort_keys=True)
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
            target.update(hits=record(Path(cell["output"])), database=record(directory / "target.dmnd"),
                          status="complete_pending_numeric_validation")
            status_path.write_text(json.dumps(status, indent=2, sort_keys=True) + "\n")
        verify_plan(manifest_path)
        for key in ("manifest", "executor", "phase_helper", "time_binary"):
            check(status[key])
        for target in status["targets"]:
            for key in ("hits", "database"):
                check(target[key])
            for phase in target["phases"].values():
                check(phase["log"])
                check(phase["gnu_time"])
        status["status"] = "complete_pending_numeric_validation"
    except BaseException as error:
        status.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        status["finished_at"] = datetime.now(timezone.utc).isoformat()
        status_path.write_text(json.dumps(status, indent=2, sort_keys=True) + "\n")
    return status


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    result = run(args.manifest.resolve(), args.check_only)
    print(result["status"])
