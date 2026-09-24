"""Bind recovered high-CPM candidate input to completed independent admission."""

import csv
import io
import json
from pathlib import Path
import subprocess

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

JOB = "22155"
COMMIT = "fd91e7432e988b7e462ab9936143fc7252ed10a4"
SOURCE_SHA = "76a49ce9cc1b65a3839fffdb5a9658220e8bf4f07756f7649573173b4d3d6a02"


def completed(accounting):
    rows = [r for r in csv.DictReader(io.StringIO(accounting), delimiter="|") if r["JobID"] == JOB]
    if len(rows) != 1 or tuple(rows[0][k] for k in ("State", "ExitCode", "AllocCPUS", "ReqMem", "NodeList")) != (
            "COMPLETED", "0:0", "2", "64G", "bizon"):
        raise ValueError("Require successful independent recovery admission 22155")
    return rows[0]


def validate_report(report, source, root):
    if (report["status"] != "cpm_checkpoint_recovered_seed_admitted_unscored"
            or report["seed_admitted"] is not True or report["source"] != source
            or report["accuracy_evaluated"] is not False or report["publication_ready"] is not False
            or report["scheduler"]["JobID"] != "22154" or report["scheduler"]["State"] != "COMPLETED"
            or report["scheduler"]["ExitCode"] != "0:0" or report["scheduler"]["AllocCPUS"] != "1"
            or report["scheduler"]["ReqMem"] != "64G" or report["scheduler"]["NodeList"] != "bizon"
            or report["comparison"]["partition_equal"] is not True or report["comparison"]["genes"] != 984137
            or report["original_failure"]["JobID"] != "22081_1" or report["original_failure"]["State"] != "FAILED"):
        raise ValueError("Wrong recovered seed admission")
    directory = root / "benchmarks/results/qfo_cpm_checkpoint_recovery_v1"
    original = root / "benchmarks/results/qfo_parameter_cpm_replay_v3/cpm_high/replay"
    expected = [dict(label=label, origin=origin, output=record(path)) for label, origin, path in (
        ("multipass", "reused", original / "orthogroups_multipass.txt"),
        ("multipass_refined", "reused", original / "orthogroups_multipass_refined.txt"),
        ("strict_profiles", "recovered", directory / "orthogroups_profiles.txt"),
        ("strict_profiles_refined", "recovered", directory / "orthogroups_profiles_refined.txt"))]
    parent = record(directory / "status.json")
    if report["stages"] != expected or report["seed_partition"] != expected[-1]["output"] or report["source_report"] != parent:
        raise ValueError("Recovered candidate seed is not the admitted refined output")
    if any(item not in report["checked_records"] for item in (source, parent, *[row["output"] for row in expected])):
        raise ValueError("Incomplete recovered seed provenance")
    return expected[-1]["output"]


def evidence(root):
    accounting = subprocess.check_output(["sacct", "-j", JOB, "--parsable2",
        "--format=JobID,State,ExitCode,Elapsed,AllocCPUS,ReqMem,NodeList"], text=True)
    scheduler = completed(accounting)
    executor = root / "benchmarks/work/cpm_checkpoint_admission_v1_20260923"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != COMMIT:
        raise ValueError("Recovery admission executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    source = record(executor / "benchmark_tools/admit_cpm_checkpoint_recovery.py")
    if source["sha256"] != SOURCE_SHA:
        raise ValueError("Recovery admission source changed")
    report_record = record(root / "benchmarks/results/qfo_cpm_checkpoint_recovery_admission_v1/status.json")
    report = json.loads(Path(report_record["path"]).read_text())
    seed = validate_report(report, source, root)
    records = [record(__file__), source, report_record, seed, *report["checked_records"]]
    for item in records:
        check(item)
    return dict(scheduler=scheduler, accounting=accounting, executor=str(executor), report_record=report_record,
                seed_partition=seed, report=report, checked_records=records)
