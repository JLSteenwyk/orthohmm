"""Independently validate job 22153 before checkpoint recovery can use it."""

import csv
import io
import json
from pathlib import Path
import subprocess
import sys

from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

JOB = "22153"
EXECUTOR = "b2d6a877d2fc26dd8b3b857687b74da4a4b1f005"
SOURCE_SHA = "6935f5f51c04188006d76465ede5647b3c2314b892c3f855e3e10404fba38719"


def completed(accounting):
    rows = [r for r in csv.DictReader(io.StringIO(accounting), delimiter="|") if r["JobID"] == JOB]
    if len(rows) != 1 or tuple(rows[0][k] for k in ("State", "ExitCode", "NodeList", "AllocCPUS")) != (
            "COMPLETED", "0:0", "bizon", "1"):
        raise ValueError("Require completed single-CPU refinement reconstruction")
    return rows[0]


def compare_partitions(expected, observed, names, expected_groups):
    if len(names) != len(set(names)):
        raise ValueError("Duplicate reconstruction gene universe")
    universe = set(names)
    left = {frozenset(group) for group in read_partition(expected, universe)}
    right = {frozenset(group) for group in read_partition(observed, universe)}
    if left != right or type(expected_groups) is not int or len(right) != expected_groups:
        raise ValueError("Independent refinement membership/count mismatch")
    return {"genes": len(names), "groups": len(right), "partition_equal": True}


def check_report(parent, child, root, executor, directory, runtime, numeric):
    source = record(executor / "benchmark_tools/reconstruct_cpm_refinement.py")
    expected = record(root / "benchmarks/results/qfo_parameter_cpm_replay_v3/cpm_high/replay/orthogroups_multipass_refined.txt")
    output = record(directory / "reconstructed.txt")
    if source["sha256"] != SOURCE_SHA:
        raise ValueError("Changed frozen reconstruction source")
    if (parent["status"] != "original_cpm_refinement_reproduced_unscored"
            or parent["job_id"] != JOB or type(parent["returncode"]) is not int or parent["returncode"] != 0
            or parent["source"] != source or parent["expected"] != expected
            or any(parent[k] is not False for k in ("recovery_authorized", "accuracy_evaluated", "publication_ready"))
            or parent["runtime_before"] != runtime or parent["runtime_after"] != runtime
            or parent["worker_report"] != record(directory / "worker.json")
            or parent["worker_log"] != record(directory / "worker.log") or parent["result"] != child):
        raise ValueError("Invalid reconstruction parent evidence")
    command = [sys.executable, "-B", source["path"], "--root", str(root),
               "--output", str(directory), "--worker"]
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    modules = [record(launcher / name) for name in (
        "benchmark_tools/replay_high_sensitivity.py", "orthohmm/accuracy.py", "orthohmm/refinement.py",
        "benchmark_tools/audit_historical_profile_ablation.py", "benchmark_tools/audit_accuracy_checkpoint.py")]
    if (parent["command"] != command or child["modules"] != modules or child["output"] != output
            or child["accuracy_evaluated"] is not False or child["partition_equal"] is not True
            or type(child["genes"]) is not int or child["genes"] != 984137
            or type(child["refinement_directed_hits"]) is not int or child["refinement_directed_hits"] != 0
            or child["numeric_checkpoint"] != {**numeric, "auditor": modules[-1]}):
        raise ValueError("Invalid reconstruction scientific input/output evidence")
    return [source, expected, output, *modules, parent["worker_report"], parent["worker_log"]]


def validate(root):
    from benchmark_tools.checked_replay_payload_worker import corrected_evidence
    from benchmark_tools.cpm_replay_context import REPLAY_SHA
    from benchmark_tools.reconstruct_cpm_refinement import PREFIX_SHA
    from benchmark_tools.run_simulation_methods import read_frozen
    from benchmark_tools.verify_qfo_replay_launcher import verify

    accounting = subprocess.check_output(["sacct", "-j", JOB, "--parsable2",
        "--format=JobID,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = completed(accounting)
    executor = root / "benchmarks/work/cpm_refinement_check_v1_20260923"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR:
        raise ValueError("Changed reconstruction executor")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--",
                    "benchmark_tools", "orthohmm"], check=True)
    plan, plan_record, admission_record, names_record = corrected_evidence(
        root / "benchmarks/work/qfo_corrected_replay_commands_20260918.json", REPLAY_SHA)
    original = read_frozen(Path(admission_record["path"]), admission_record["sha256"])
    numeric = original["content"]["numeric_checkpoint"]
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    core = root / "benchmarks/work/publication_method_native_v2"
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    runtime = verify(core, launcher, runtime_path)
    if runtime != plan["runtime"]:
        raise ValueError("Reconstruction runtime differs from frozen plan")
    prefix_path = root / "benchmark_tools/results/qfo_cpm_recovery_predecessors_20260923.json"
    prefix = read_frozen(prefix_path, PREFIX_SHA)
    directory = root / f"benchmarks/work/qfo_cpm_refinement_check_{JOB}"
    status_record = record(directory / "status.json")
    worker_record = record(directory / "worker.json")
    parent, child = [read_frozen(Path(item["path"]), item["sha256"]) for item in (status_record, worker_record)]
    bound = check_report(parent, child, root, executor, directory, runtime, numeric)
    required = [record(prefix_path), *prefix["checked_records"], plan_record, admission_record, names_record]
    if any(item not in parent["checked_records"] for item in required):
        raise ValueError("Missing reconstruction predecessor/input binding")
    records = [record(__file__), status_record, worker_record, *bound, *required, *parent["checked_records"]]
    unique = {}
    for item in records:
        if item["path"] in unique and unique[item["path"]] != item:
            raise ValueError("Conflicting reconstruction provenance")
        unique[item["path"]] = item
        check(item)
    names = Path(names_record["path"]).read_text().splitlines()
    if len(names) != 984137:
        raise ValueError("Wrong reconstruction gene count")
    comparison = compare_partitions(Path(parent["expected"]["path"]), Path(child["output"]["path"]),
                                    names, child["groups"])
    for item in unique.values():
        check(item)
    if verify(core, launcher, runtime_path) != runtime:
        raise ValueError("Runtime changed during reconstruction validation")
    return {"status": "original_cpm_refinement_independently_verified", "scheduler": scheduler,
            "accounting": accounting, "source": record(__file__), "comparison": comparison,
            "checked_records": list(unique.values()), "recovery_authorized": False,
            "accuracy_evaluated": False, "publication_ready": False,
            "limitations": ["Existing multipass refinement only; failed stage and resumed optimization require separate gates.",
                            "Reconstruction is not controlled timing or a new accuracy evaluation."]}
