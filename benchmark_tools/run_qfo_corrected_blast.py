"""Execute fresh legacy formatdb/BLAST stages; leave inference admission separate."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.prepare_qfo_corrected_orthomcl import commands, require_parity
from benchmark_tools.run_qfo_corrected_sonic import verify_corrected_inputs
from benchmark_tools.prepare_qfo_corrected_proteinortho import PRIMARY_SHA
from benchmark_tools.snapshot_runtime_trees import verify as verify_tree

PLAN_SHA = "3106012dc42c053d42e7f9d8d08532168d5826aa6812bab4b4c232f900e3a8ff"
RUNTIME_SHA = "9ce46b5329f34384b03980b62fbfe9522f244dc227ce6506071ce7000b154e42"


def environment():
    return {"HOME": str(Path.home()), "USER": "bizon", "LOGNAME": "bizon",
            "PATH": "/usr/bin:/bin", "LANG": "C", "LC_ALL": "C"}


def verify(plan_path, runtime_path):
    plan = read_frozen(plan_path, PLAN_SHA)
    if plan["status"] != "corrected_orthomcl_inputs_prepared_pending_execution_freeze" or plan["execution_authorized"] is not False:
        raise ValueError("Unexpected prepared input state")
    require_parity(plan["independent_parity"])
    for item in [plan["source"], *plan["checked_records"], *plan["prepared_inputs"]]:
        check(item)
    primary = read_frozen(Path(plan["primary_manifest"]["path"]), PRIMARY_SHA)
    verify_corrected_inputs(primary, plan["input_fastas"])
    root = Path(plan["output_root"])
    expected = commands(root / "work", plan["tools"]["blastall"]["path"], plan["tools"]["formatdb"]["path"])
    if expected != plan["search_commands"]:
        raise ValueError("Changed legacy search command")
    for path in (root / "work/.ncbirc", Path.home() / ".ncbirc", Path("/etc/.ncbirc"),
                 Path("/etc/ncbi.ini"), Path("/etc/ncbi"), Path("/etc/ld.so.preload")):
        if path.exists() or path.is_symlink():
            raise ValueError("Unreviewed configuration or loader override: " + str(path))
    verify_tree(read_frozen(runtime_path, RUNTIME_SHA))
    return plan


def run(plan_path, runtime_path, check_only=False):
    if not check_only and (os.environ.get("SLURM_CPUS_PER_TASK") != "180" or not os.environ.get("SLURM_JOB_ID")):
        raise ValueError("Require scheduled 180-CPU allocation")
    plan = verify(plan_path, runtime_path)
    root = Path(plan["output_root"])
    work = root / "work"
    execution = root / "search_execution"
    if execution.exists() or {p.name for p in work.iterdir()} != {"all.fa", "all.gg"}:
        raise FileExistsError("Search outputs already present; no implicit resume")
    if check_only:
        return {"status": "preflight_passed_no_search"}
    execution.mkdir(exist_ok=False)
    status = execution / "status.json"
    report = {"status": "preparing", "source": record(__file__), "plan": record(plan_path),
              "runtime": record(runtime_path), "job_id": os.environ["SLURM_JOB_ID"],
              "node": os.uname().nodename, "cwd": str(work), "environment": environment(),
              "commands": plan["search_commands"], "stages": {}, "accuracy_admitted": False,
              "search_admitted": False,
              "limitations": ["Shared-host stage timing, not matched end-to-end timing.",
                              "Zero exit does not exclude sequence-specific BLAST errors; query/log audit required.",
                              "Database, BPO, clustering and final-group conversion require separate admission."]}
    try:
        for stage in ("formatdb", "blast"):
            report.update(status="running_" + stage)
            started = time.time()
            report["stages"][stage] = {"started_epoch": started}
            status.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
            with (execution / (stage + ".log")).open("xb") as log:
                result = subprocess.run(["/usr/bin/time", "-v", "-o", str(execution / (stage + ".time.txt")),
                                         *plan["search_commands"][stage]], cwd=work, env=environment(),
                                        stdout=log, stderr=subprocess.STDOUT)
            report["stages"][stage].update(exit_code=result.returncode, finished_epoch=time.time(),
                log=record(execution / (stage + ".log")), timing=record(execution / (stage + ".time.txt")))
            if result.returncode:
                raise RuntimeError(f"{stage} exited {result.returncode}")
            verify(plan_path, runtime_path)
            if stage == "formatdb":
                database = [work / ("all.fa." + suffix) for suffix in ("phr", "pin", "psq")]
                if any(not p.is_file() or not p.stat().st_size for p in database):
                    raise ValueError("Incomplete formatted database")
                report["database_files"] = [record(p) for p in sorted(work.glob("all.fa.*")) if p.is_file()]
            else:
                for item in report["database_files"]:
                    check(item)
                partial = work / "all.blast.partial"
                if not partial.is_file() or not partial.stat().st_size:
                    raise ValueError("Empty BLAST output")
                final = work / "all.blast"
                if final.exists():
                    raise FileExistsError(final)
                partial.rename(final)
                report["blast_output"] = record(final)
        report["status"] = "search_exited_zero_pending_query_and_database_admission"
    except Exception as exc:
        report.update(status="failed", error=str(exc), finished_epoch=time.time())
        raise
    finally:
        status.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--runtime", type=Path, required=True)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    print(json.dumps({"status": run(args.plan.resolve(), args.runtime.resolve(), args.check_only)["status"]}))
