"""Run the two retained post-completion audits once, without rerunning inference."""

import argparse
import csv
import io
import json
import os
from pathlib import Path
import signal
import subprocess
import sys

from benchmark_tools.run_integrated_full_job import record, save
from benchmark_tools.run_integrated_publication_workflow import check

JOBS = [22377, 22378]


def protocol(path, digest):
    if record(path)["sha256"] != digest:
        raise ValueError("Completion-audit protocol differs")
    value = json.loads(path.read_text())
    if (value["schema"] != "publication_completion_audits_v1" or value["jobs"] != JOBS
            or value["production_timing"] is not False or type(value["attempts_per_stage"]) is not int
            or value["attempts_per_stage"] != 1
            or value["stage_timeout_seconds"] != 7200
            or [v["name"] for v in value["stages"]] != ["restored_orthobench", "observer_calibration"]):
        raise ValueError("Completion-audit scope differs")
    if [v["parent_job"] for v in value["stages"]] != JOBS or not Path(value["output"]).is_absolute():
        raise ValueError("Completion-audit identities differ")
    for pin in [value["interpreter"], value["submission_script"], *value["sources"]]:
        check(pin)
    if record(__file__) not in value["sources"]:
        raise ValueError("Audit runner not frozen")
    return value


def scheduler_success(raw, job):
    reader = csv.DictReader(io.StringIO(raw), delimiter="|")
    fields = reader.fieldnames or []
    if len(set(fields)) != len(fields) or not {"JobIDRaw", "State", "ExitCode"} <= set(fields):
        return False
    values = list(reader)
    if any(not all(isinstance(row.get(key), str) for key in ("JobIDRaw", "State", "ExitCode")) for row in values):
        return False
    rows = [row for row in values if row["JobIDRaw"] == str(job) or row["JobIDRaw"].startswith(str(job) + ".")]
    names = [row["JobIDRaw"] for row in rows]
    if len(set(names)) != len(names) or not {str(job), f"{job}.batch"} <= set(names):
        return False
    return all(row["State"] == "COMPLETED" and row["ExitCode"] == "0:0" for row in rows)


def execute(stage, directory, timeout=7200, *, parent_success=False):
    path = directory / (stage["name"] + ".log")
    result = dict(name=stage["name"], command=stage["command"], returncode=None,
                  timed_out=False, result=None, semantic_check_passed=False,
                  parent_scheduler_success=parent_success)
    child = None
    try:
        if Path(stage["result_path"]).exists():
            raise FileExistsError("Existing audit result; do not repeat")
        with path.open("x") as log:
            child = subprocess.Popen(stage["command"], stdout=log, stderr=subprocess.STDOUT,
                                     start_new_session=True)
            try:
                result["returncode"] = child.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                result["timed_out"] = True
                os.killpg(child.pid, signal.SIGKILL)
                result["returncode"] = child.wait()
        output = Path(stage["result_path"])
        if output.is_file():
            pin = record(output)
            value = json.loads(output.read_text())
            check(pin)
            result["result"] = pin
            result["reported_status"] = value.get("status")
            if value.get("job_id") != stage["parent_job"]:
                raise ValueError("Audit result belongs to another job")
            if stage["name"] == "restored_orthobench":
                result["reproduction_equal"] = value.get("reproduction_equal")
                result["semantic_check_passed"] = (value.get("status") == "restored_archive_orthobench_independently_audited"
                                                   and value.get("reproduction_equal") is True)
            else:
                result["semantic_check_passed"] = value.get("status") == "calibration_checks_passed"
    except Exception as error:
        result.update(error_type=type(error).__name__, error=str(error))
    finally:
        if child is not None and child.poll() is None:
            os.killpg(child.pid, signal.SIGKILL)
            child.wait()
        if path.exists():
            result["log"] = record(path)
        result["passed"] = (result["returncode"] == 0 and not result["timed_out"]
                            and result["semantic_check_passed"] and parent_success is True)
        save(directory / (stage["name"] + ".json"), result)
    return result


def run(protocol_path, digest):
    value = protocol(protocol_path, digest)
    if (os.uname().nodename != "bizon" or os.environ.get("SLURM_CPUS_PER_TASK") != "4"
            or os.environ.get("SLURM_MEM_PER_NODE") != "32768" or not os.environ.get("SLURM_JOB_ID")
            or record(sys.executable) != value["interpreter"]):
        raise ValueError("Require the pinned local 4-CPU/32-GiB audit allocation")
    directory = Path(value["output"])
    directory.mkdir(parents=True, exist_ok=False)
    result = dict(status="completion_audits_incomplete", protocol=record(protocol_path), source=record(__file__),
                  job_id=os.environ["SLURM_JOB_ID"], stages=[], production_timing=False, publication_ready=False)
    try:
        argv = ["sacct", "-j", ",".join(map(str, JOBS)), "--parsable2",
                "--format=JobIDRaw,State,ExitCode,Elapsed,AllocCPUS,NodeList,ReqMem"]
        observed = subprocess.run(argv, capture_output=True, text=True, timeout=15, check=True)
        save(directory / "parent_accounting.json", dict(command=argv, stdout=observed.stdout, stderr=observed.stderr))
        for stage in value["stages"]:
            # Both audit attempts remain independent; one failure does not suppress the other.
            stage_result = execute(stage, directory, value["stage_timeout_seconds"],
                                   parent_success=scheduler_success(observed.stdout, stage["parent_job"]))
            result["stages"].append(stage_result)
        protocol(protocol_path, digest)
        result["status"] = "completion_audits_passed" if all(v["passed"] for v in result["stages"]) else "completion_audits_failed"
        result["accounting"] = record(directory / "parent_accounting.json")
    except Exception as error:
        result.update(error_type=type(error).__name__, error=str(error))
        raise
    finally:
        save(directory / "complete.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--protocol-sha256", required=True)
    args = parser.parse_args()
    result = run(args.protocol.resolve(), args.protocol_sha256)
    print(result["status"])
    raise SystemExit(0 if result["status"] == "completion_audits_passed" else 1)
