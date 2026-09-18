"""Bind terminal corrected BLAST provenance to database and full-table audits."""

import argparse
import json
import math
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_qfo_corrected_blast import verify, environment
from benchmark_tools.audit_orthomcl_database import audit as audit_database
from benchmark_tools.audit_orthomcl_search_table import audit as audit_table
from benchmark_tools.verify_ygob_validation import require_completed_job

EXECUTOR = "b0c7128d464e2abe5bf1fd3646c0378625bbb0de"


def validate_execution(plan, report, scheduler, plan_record, runtime_record, source):
    if (scheduler["State"] != "COMPLETED" or scheduler["ExitCode"] != "0:0"
            or scheduler["NodeList"] != "bizon" or scheduler["AllocCPUS"] != "180"
            or scheduler["ReqMem"] != "900G" or report["job_id"] != scheduler["JobIDRaw"]):
        raise ValueError("Wrong BLAST scheduler identity/allocation")
    if (report["status"] != "search_exited_zero_pending_query_and_database_admission"
            or report["search_admitted"] is not False or report["accuracy_admitted"] is not False):
        raise ValueError("Require successful unadmitted BLAST execution")
    work = Path(plan["output_root"]) / "work"
    execution = Path(plan["output_root"]) / "search_execution"
    if (report["plan"] != plan_record or report["runtime"] != runtime_record
            or report["source"] != source or report["commands"] != plan["search_commands"]
            or report["environment"] != environment() or report["cwd"] != str(work)
            or report["node"] != "bizon"):
        raise ValueError("Changed BLAST execution provenance")
    if set(report["stages"]) != {"formatdb", "blast"}:
        raise ValueError("Wrong native stage set")
    last_finish = 0
    checked = []
    for name in ("formatdb", "blast"):
        stage = report["stages"][name]
        start, finish = stage["started_epoch"], stage["finished_epoch"]
        if (type(stage["exit_code"]) is not int or stage["exit_code"] != 0
                or any(type(v) not in (int, float) or not math.isfinite(v) for v in (start, finish))
                or not finish >= start > 0 or start < last_finish):
            raise ValueError("Invalid native stage status/timestamps")
        last_finish = finish
        for key, filename in (("log", name + ".log"), ("timing", name + ".time.txt")):
            item = stage[key]
            if (item["path"] != str(execution / filename) or type(item["bytes"]) is not int
                    or item["bytes"] < (1 if key == "timing" else 0)):
                raise ValueError("Wrong native log/timing artifact")
            checked.append(item)
    expected_paths = [str(work / ("all.fa." + suffix)) for suffix in ("phr", "pin", "psq")]
    if [item["path"] for item in report["database_files"]] != expected_paths:
        raise ValueError("Wrong formatted database inventory")
    if (report["blast_output"]["path"] != str(work / "all.blast")
            or report["blast_output"]["bytes"] <= 0):
        raise ValueError("Wrong final BLAST table")
    return [*checked, *report["database_files"], report["blast_output"]]


def admit(root, job, output):
    if output.exists():
        raise FileExistsError(output)
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS,ReqMem"], text=True)
    scheduler = require_completed_job(accounting, job)
    executor = root / "benchmarks/work/publication_qfo_corrected_blast_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR:
        raise ValueError("Changed frozen BLAST executor")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    plan_path = root / "benchmark_tools/results/qfo_corrected_orthomcl_prepared_20260918.json"
    runtime_path = root / "benchmark_tools/results/qfo_corrected_legacy_blast_runtime_20260918.json"
    plan = verify(plan_path, runtime_path)
    path = Path(plan["output_root"]) / "search_execution/status.json"
    execution_record = record(path)
    execution = json.loads(path.read_text())
    source = record(executor / "benchmark_tools/run_qfo_corrected_blast.py")
    native = validate_execution(plan, execution, scheduler, record(plan_path), record(runtime_path), source)
    checked = [execution_record, record(plan_path), record(runtime_path), source, plan["source"],
               *plan["checked_records"], *plan["prepared_inputs"], *native,
               *[record(Path(__file__).with_name(name)) for name in (
                   "run_qfo_corrected_blast.py", "audit_orthomcl_database.py", "audit_orthomcl_search_table.py")]]
    work = Path(plan["output_root"]) / "work"
    if (work / "all.blast.partial").exists():
        raise ValueError("Unexpected partial output alongside final BLAST table")
    if (work / "formatdb.log").is_file():
        checked.append(record(work / "formatdb.log"))
    for item in checked:
        check(item)
    output.mkdir(parents=True, exist_ok=False)
    report = {"status": "running_validation", "source": record(__file__), "scheduler": scheduler,
              "accounting": accounting, "execution_report": execution_record, "checked_records": checked,
              "search_admitted": False, "accuracy_admitted": False, "publication_ready": False,
              "downstream_execution_authorized": False}
    try:
        database = audit_database(work / "all.fa", runtime_path, output / "database")
        report["database_audit"] = record(output / "database/report.json")
        if (database["status"] != "database_exact_sequence_parity_verified"
                or database["content"]["input_sequences"] != 984137):
            raise ValueError("Corrected database parity requires review")
        table = audit_table(work / "all.blast", work / "all.fa",
                            Path(execution["stages"]["blast"]["log"]["path"]), output / "table.json")
        report["table_audit"] = record(output / "table.json")
        if table["content"]["input_proteins"] != 984137:
            raise ValueError("Wrong corrected table input universe")
        if verify(plan_path, runtime_path) != plan:
            raise ValueError("Changed prepared search provenance")
        for item in [*checked, report["source"], report["database_audit"], report["table_audit"],
                     *database["outputs"], *database["checked_records"], *table["checked_records"]]:
            check(item)
        report.update(status="corrected_orthomcl_search_evidence_verified", search_admitted=True,
                      query_coverage=table["content"], database_content=database["content"])
        report["limitations"] = [
            "Search evidence admission preserves logged failed queries; it does not assert success for every input query.",
            "No-hit queries are not automatically failures, and incoming hits do not repair outgoing search failures.",
            "BPO/index validation, final OrthoMCL groups and reference impact/scoring remain separate steps.",
            "No automatic downstream execution is authorized by this report.",
            "Shared-host stage accounting is not matched-resource end-to-end timing."]
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with (output / "report.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = admit(args.root.resolve(), args.job, args.output.resolve())
    print(json.dumps({"status": result["status"]}))
