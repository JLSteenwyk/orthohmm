"""Admit a complete corrected QfO search panel before numeric hit validation."""

import argparse
import json
import math
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_qfo_sequence_search_control import verify_plan, MANIFEST_SHA
from benchmark_tools.verify_ygob_validation import require_completed_job

EXECUTOR = "b66d3f374224b89f782c98a4b29a9c33b6cefdc7"


def validate_execution(report, plan, scheduler, manifest, executor):
    if ((scheduler["State"], scheduler["ExitCode"], scheduler["NodeList"], scheduler["AllocCPUS"])
            != ("COMPLETED", "0:0", "bizon", "32") or scheduler["ReqMem"] not in ("192G", "192Gn")):
        raise ValueError("Require successful scheduled 32-CPU/192-GiB search")
    expected_env = {"PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    if (report["status"] != "complete_pending_numeric_validation" or report["job_id"] != scheduler["JobIDRaw"]
            or report["node"] != "bizon" or report["manifest"] != manifest
            or report["executor"] != record(executor / "benchmark_tools/run_qfo_sequence_search_control.py")
            or report["phase_helper"] != record(executor / "benchmark_tools/run_sequence_search_control.py")
            or report["time_binary"] != record(Path("/usr/bin/time"))
            or report["environment_overrides"] != expected_env
            or report["accuracy_evaluated"] is not False or report["numeric_validated"] is not False
            or report["publication_ready"] is not False or len(report["targets"]) != 78
            or len(plan["searches"]) != 78):
        raise ValueError("Search execution is incomplete or provenance differs")
    records = [report[k] for k in ("manifest", "executor", "phase_helper", "time_binary")]
    for expected, target in zip(plan["searches"], report["targets"], strict=True):
        directory = Path(expected["output"]).parent
        if (type(target["index"]) is not int or target["index"] != expected["index"]
                or target["status"] != "complete_pending_numeric_validation"
                or set(target["phases"]) != {"makedb", "search"}
                or target["hits"]["path"] != expected["output"]
                or target["database"]["path"] != str(directory / "target.dmnd")):
            raise ValueError("Target inventory, order or completion differs")
        records.extend([target["hits"], target["database"]])
        for name in ("makedb", "search"):
            phase = target["phases"][name]
            if (phase["argv"] != expected[name] or type(phase["exit_code"]) is not int or phase["exit_code"] != 0
                    or type(phase["wall_s"]) not in (int, float) or not math.isfinite(phase["wall_s"]) or phase["wall_s"] < 0
                    or phase["log"]["path"] != str(directory / (name + ".log"))
                    or phase["gnu_time"]["path"] != str(directory / (name + ".time.log"))):
                raise ValueError("Search phase command, exit or evidence differs")
            records.extend([phase["log"], phase["gnu_time"]])
    return records


def admit(root, job, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS,ReqMem"], text=True)
    scheduler = require_completed_job(accounting, job)
    executor = root / "benchmarks/work/publication_qfo_sequence_search_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR:
        raise ValueError("Frozen search executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    source = Path(__file__).with_name("run_qfo_sequence_search_control.py")
    if record(source)["sha256"] != record(executor / "benchmark_tools/run_qfo_sequence_search_control.py")["sha256"]:
        raise ValueError("Search plan verifier differs from frozen runner")
    manifest = root / "benchmarks/work/qfo_sequence_search_control_v1/manifest.json"
    plan = verify_plan(manifest)
    execution_path = manifest.parent / "execution.json"
    execution_record = record(execution_path)
    execution = json.loads(execution_path.read_text())
    records = validate_execution(execution, plan, scheduler, record(manifest), executor)
    checked = [execution_record, *records, *plan["checked_records"], plan["queries"], plan["gene_metadata"]]
    for item in checked:
        check(item)
    verify_plan(manifest)
    for item in checked:
        check(item)
    result = {"status": "corrected_qfo_search_panel_admitted_pending_numeric_validation",
        "source": record(__file__), "executor_commit": EXECUTOR, "manifest": record(manifest),
        "manifest_sha256": MANIFEST_SHA, "execution": execution_record, "scheduler": scheduler,
        "accounting": accounting, "targets": execution["targets"], "genes": plan["genes"],
        "proteomes": plan["proteomes"], "checked_records": checked,
        "numeric_validated": False, "accuracy_evaluated": False, "publication_ready": False,
        "limitations": ["Validates complete panel execution and retained file identities, not hit content or biological completeness.",
            "Numeric validation must check IDs, lengths, target ownership, duplicate pairs and score fields before graph replay.",
            "Phase wall and GNU-time logs are retained; shared-host cost is not matched dedicated timing evidence."]}
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--job", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.job, args.output.absolute())
