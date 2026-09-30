"""Explicit seed handoff for the separately admitted helper-runtime recovery."""

import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

READBACK = "benchmark_tools/results/qfo_cpm_helper_recovery_admission_readback_20260930.json"
READBACK_SHA = "5fb31ea6cee69501d5ad1826789e9a178f230bb909012537e123d18d21d65850"
ADMISSION = "benchmarks/work/qfo_cpm_helper_recovery_admission_20260930/status.json"
ADMISSION_SHA = "35febb4c1810892988beb79bc8e6ac2c6a2e1ffbf1b111492ad1649f39e04997"
SOURCE_SHA = "6d01dec48b972f7c090ca2c72b08e90da315c5cde3aa9d6b56b92d00ed9485d4"
REVISION = "c2a0137fa197449e7bc17e061bb54488012163b0"
PROTOCOL = "benchmark_tools/results/QFO_CPM_HELPER_CANDIDATE_PROTOCOL_20260930.md"


def scheduler_rows(accounting):
    rows = list(csv.DictReader(io.StringIO(accounting), delimiter="|"))
    selected = {}
    for job, expected in (
        ("22154", ("COMPLETED", "0:0", "1", "64G", "bizon")),
        ("22081_1", ("FAILED", "1:0", "32", "192G", "bizon")),
        ("22155", ("FAILED", "1:0", "2", "64G", "bizon")),
    ):
        matches = [row for row in rows if row["JobID"] == job]
        if len(matches) != 1 or tuple(matches[0][k] for k in
                ("State", "ExitCode", "AllocCPUS", "ReqMem", "NodeList")) != expected:
            raise ValueError("Changed recovery success/failure identity: " + job)
        selected[job] = matches[0]
    cancelled = [row for row in rows if row["JobID"] == "22156"]
    if (len(cancelled) != 1 or cancelled[0]["State"].split()[0] != "CANCELLED"
            or cancelled[0]["ExitCode"] != "0:0" or cancelled[0]["AllocCPUS"] != "0"):
        raise ValueError("Original candidate 22156 must remain cancelled without execution")
    selected["22156"] = cancelled[0]
    return selected


def validate_report(report, readback, root, scheduler):
    if (report["status"] != "cpm_helper_runtime_recovered_seed_admitted_unscored"
            or readback["status"] != "helper_runtime_recovered_seed_admission_independently_read_back_unscored"
            or report["source"] != record(root / "benchmark_tools/admit_helper_cpm_recovery.py")
            or report["source"]["sha256"] != SOURCE_SHA or readback["source_revision"] != REVISION):
        raise ValueError("Wrong explicit helper-runtime seed admission")
    for item in (report, readback):
        if (item["seed_admitted"] is not True or item["native_attempts"] != 0
                or any(item[k] is not False for k in
                       ("accuracy_evaluated", "downstream_admitted", "publication_ready"))):
            raise ValueError("Seed-only admission flags differ")
        for field, job in (("scheduler", "22154"), ("original_failure", "22081_1"),
                           ("historical_failed_admission", "22155")):
            if any(scheduler[job][k] != value for k, value in item[field].items()):
                raise ValueError("Historical scheduler evidence differs")
    amendment = report["refinement_runtime_amendment"]
    if (amendment["scope"] != "independent refinement only"
            or amendment["original_optimizer_runtime_unchanged"] is not True
            or amendment["downstream_runtime_authorized"] is not False):
        raise ValueError("Refinement amendment cannot authorize a candidate runtime change")
    directory = root / "benchmarks/results/qfo_cpm_checkpoint_recovery_v1"
    original = root / "benchmarks/results/qfo_parameter_cpm_replay_v3/cpm_high/replay"
    stages = [dict(label=label, origin=origin, output=record(path)) for label, origin, path in (
        ("multipass", "reused", original / "orthogroups_multipass.txt"),
        ("multipass_refined", "reused", original / "orthogroups_multipass_refined.txt"),
        ("strict_profiles", "recovered", directory / "orthogroups_profiles.txt"),
        ("strict_profiles_refined", "recovered", directory / "orthogroups_profiles_refined.txt"))]
    coverage = [dict(label=stage["label"], genes=984137, groups=groups) for stage, groups in
                zip(stages, (316603, 393142, 314274, 390845))]
    if (report["stages"] != stages or readback["coverage"] != coverage
            or report["coverage"] != [dict(row, output=stage["output"]) for row, stage in zip(coverage, stages)]
            or report["seed_partition"] != stages[-1]["output"]
            or readback["seed_partition"] != stages[-1]["output"]
            or report["source_report"] != record(directory / "status.json")
            or readback["saved_graph"] != report["saved_graph"]
            or readback["constructor_bytes_sha256"] != report["constructor_bytes_sha256"]
            or readback["checked_records_reverified"] != len(report["checked_records"])):
        raise ValueError("Independent graph/coverage/seed lineage differs")
    required = [report["source"], report["source_report"], *[stage["output"] for stage in stages]]
    if any(item not in report["checked_records"] for item in required):
        raise ValueError("Incomplete admitted seed provenance")
    return stages[-1]["output"]


def evidence(root, readback_path, readback_sha, protocol_sha):
    root = Path(root).resolve()
    readback_record = record(readback_path)
    if (Path(readback_record["path"]) != root / READBACK
            or readback_record["sha256"] != readback_sha or readback_sha != READBACK_SHA):
        raise ValueError("Require explicitly selected independent seed readback")
    protocol = record(root / PROTOCOL)
    if protocol["sha256"] != protocol_sha:
        raise ValueError("Unreviewed candidate handoff protocol")
    readback = json.loads(Path(readback_record["path"]).read_bytes())
    report_record = record(root / ADMISSION)
    if report_record != readback["admission"] or report_record["sha256"] != ADMISSION_SHA:
        raise ValueError("Amended recovery admission changed")
    report = json.loads(Path(report_record["path"]).read_bytes())
    accounting = subprocess.check_output(["sacct", "-j", "22081,22154,22155,22156", "-X", "--parsable2",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,AllocCPUS,ReqMem,NodeList"], text=True)
    scheduler = scheduler_rows(accounting)
    seed = validate_report(report, readback, root, scheduler)
    bindings = readback["source_git_bindings"]
    expected_paths = {root / relative for relative in (
        "benchmark_tools/admit_helper_cpm_recovery.py", "tests/unit/test_admit_helper_cpm_recovery.py",
        "benchmark_tools/results/QFO_CPM_HELPER_RECOVERY_ADMISSION_PROTOCOL_20260930.md",
        "benchmark_tools/probe_cpm_partition_parser.py", "benchmark_tools/probe_cpm_private_runtime.py")}
    if len(bindings) != 5 or {Path(item["path"]) for item in bindings} != expected_paths:
        raise ValueError("Missing fixed admission source bindings")
    for item in bindings:
        if item["git_revision"] != REVISION:
            raise ValueError("Wrong admission source revision")
        content = subprocess.check_output(["git", "-C", str(root), "show",
            REVISION + ":" + str(Path(item["path"]).relative_to(root))])
        if len(content) != item["bytes"] or hashlib.sha256(content).hexdigest() != item["sha256"]:
            raise ValueError("Admission source Git blob differs")
    records = [record(__file__), protocol, readback_record, report_record, seed,
               *report["checked_records"], *[{k: item[k] for k in ("path", "bytes", "sha256")} for item in bindings]]
    for item in records:
        check(item)
    return dict(handoff="explicit_helper_runtime_seed_amendment", accounting=accounting,
        scheduler=scheduler, report=report, report_record=report_record, readback=readback_record,
        protocol=protocol, seed_partition=seed, checked_records=records,
        candidate_runtime_policy="unchanged historical runtime, verified against other parameter arms",
        accuracy_evaluated=False, publication_ready=False)
