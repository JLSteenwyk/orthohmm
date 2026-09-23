"""Validate CPM-derived candidate families before inferred phylogeny."""

import argparse
import csv
import io
import json
import math
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_accuracy_checkpoint import audit as audit_numeric
from benchmark_tools.audit_candidate_arm import audit as audit_arm
from benchmark_tools.checked_replay_payload_worker import corrected_evidence
from benchmark_tools.cpm_replay_context import evidence, REPLAY_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_qfo_cpm_candidates import replay_evidence, BASELINE_SHA
from benchmark_tools.run_qfo_cpm_variant import ARMS
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.trace_ob_families import partition, validate_merge_reconstruction
from benchmark_tools.verify_qfo_replay_launcher import verify

JOB = "22062"
EXECUTOR = "b32ced4635a59abd71b25ae2d34efa4f4b29bd51"
SOURCE_SHA = "192d8e2b96001093d06bd15b82a64cc75a9a4e9cb4eda877d5424776de8e2b9c"


def completed_task(accounting, index):
    if type(index) is not int or index not in range(2):
        raise ValueError("Unknown CPM candidate index")
    rows = [r for r in csv.DictReader(io.StringIO(accounting), delimiter="|") if r["JobID"] == f"{JOB}_{index}"]
    if len(rows) != 1 or tuple(rows[0][k] for k in ("State", "ExitCode", "NodeList", "AllocCPUS")) != (
            "COMPLETED", "0:0", "bizon", "2"):
        raise ValueError("Require completed 2-CPU CPM candidate task")
    return rows[0]


def validate_report(report, index, scheduler, context, parameters, source, inputs):
    if (report["status"] != "cpm_candidates_prepared_pending_admission" or report["arm"] != ARMS[index]
            or type(report["index"]) is not int or report["index"] != index
            or report["array_task_id"] != str(index) or report["job_id"] != scheduler["JobIDRaw"]
            or report["executor_commit"] != EXECUTOR or report["source"] != source
            or report["context"] != context or report["inputs"] != inputs
            or report["accuracy_evaluated"] is not False or report["publication_ready"] is not False):
        raise ValueError("Candidate preparation identity/provenance differs")
    control = report["candidate_parameter_control"]
    if (control["label"] != "control" or control["delta"] != {} or control["applied_parameters"] != parameters
            or type(control["engine_calls"]) is not int or control["engine_calls"] != 1
            or control["engine_fixed_profile_report"]["parameters"] != parameters
            or report["candidate_arm"]["expansion"] != control["engine_fixed_profile_report"]):
        raise ValueError("CPM candidates changed fixed expansion parameters")
    seconds = report["incremental_seconds"]
    if type(seconds) not in (int, float) or not math.isfinite(seconds) or seconds < 0:
        raise ValueError("Invalid candidate incremental time")


def admit(root, index, destination):
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    accounting = subprocess.check_output(["sacct", "-j", JOB, "--parsable2",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = completed_task(accounting, index)
    own_helpers = [record(p) for p in sorted(Path(__file__).parent.glob("*.py"))]
    decoder_source = record(sys.modules["orthohmm.accuracy"].__file__)
    plan_path = root / "benchmarks/work/qfo_corrected_replay_commands_20260918.json"
    plan, plan_record, _, names_record = corrected_evidence(plan_path, REPLAY_SHA)
    context = evidence(root, plan, plan_record, ARMS[index])
    replay = replay_evidence(root, index, context)
    executor = root / "benchmarks/work/publication_qfo_cpm_candidates_v3"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR:
        raise ValueError("Candidate executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    source = record(executor / "benchmark_tools/prepare_qfo_cpm_candidates.py")
    if source["sha256"] != SOURCE_SHA:
        raise ValueError("Candidate source changed")
    directory = root / "benchmarks/results/qfo_cpm_candidates_v1" / ARMS[index]
    preparation = record(directory / "manifest.json")
    report = read_frozen(Path(preparation["path"]), preparation["sha256"])
    baseline_path = root / "benchmarks/results/qfo_corrected_factorial_v1/manifest.json"
    baseline = read_frozen(baseline_path, BASELINE_SHA)
    parameters = baseline["candidate_arms"]["p1_c1"]["expansion"]["parameters"]
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    helpers = [record(p) for p in sorted((executor / "benchmark_tools").glob("*.py"))]
    inputs = [plan_record, record(baseline_path), record(runtime_path), *context["checked_records"],
              *replay["checked_records"], *helpers]
    validate_report(report, index, scheduler, context, parameters, source, inputs)
    fresh_path = directory / "fresh_replay_admission.json"
    fresh_record = record(fresh_path)
    command = [sys.executable, "-B", str(Path(replay["executor"]) / "benchmark_tools/admit_qfo_cpm_variant.py"),
               "--root", str(root), "--index", str(index), "--output", str(fresh_path)]
    if (report["replay_admission"] != replay or report["fresh_replay_admission"] != fresh_record
            or report["admission_command"] != command
            or read_frozen(fresh_path, fresh_record["sha256"]) != replay["report"]):
        raise ValueError("Candidate replay authorization differs")
    core, launcher = root / "benchmarks/work/publication_method_native_v2", Path(context["cwd"])
    runtime = verify(core, launcher, runtime_path)
    if not runtime == plan["runtime"] == baseline["runtime_before"] == baseline["runtime_after"] == report["runtime_before"] == report["runtime_after"]:
        raise ValueError("Candidate runtime differs")
    checkpoint = plan["checkpoint_manifest"]
    numeric = audit_numeric(Path(checkpoint["path"]).parent, checkpoint["sha256"])
    prior = report["numeric_checkpoint"]
    if (any(numeric[k] != prior[k] for k in ("status", "manifest", "summary", "accuracy_evaluated"))
            or prior["auditor"] != record(executor / "benchmark_tools/audit_accuracy_checkpoint.py")
            or numeric["summary"] != baseline["numeric_checkpoint"]["summary"]):
        raise ValueError("Candidate numeric checkpoint differs")
    names = Path(names_record["path"]).read_text().splitlines()
    universe = set(names)
    if len(names) != 984137 or len(universe) != len(names):
        raise ValueError("Wrong candidate gene universe")
    arm, seed = report["candidate_arm"], replay["seed_partition"]
    content = audit_arm(arm, seed, directory / "candidate", universe, True)
    if content != report["content_audit"]:
        raise ValueError("Independent candidate content differs")
    seeds, _ = partition(Path(seed["path"]), "plain", universe)
    groups, _ = partition(Path(arm["candidate_partition"]["path"]), "plain", universe)
    events = json.loads(Path(arm["membership_constraints"]["path"]).read_text())
    validate_merge_reconstruction(events, seeds, groups, set())
    records = [preparation, fresh_record, names_record, source, decoder_source, *inputs, *own_helpers,
        *arm["output_files"], *content["checked_records"], record(directory / "admission.log"),
        record(root / f"benchmarks/work/qfo_cpm_candidates_{JOB}_{index}.time.txt")]
    unique = {}
    for item in records:
        if item["path"] in unique and unique[item["path"]] != item:
            raise ValueError("Conflicting candidate provenance")
        unique[item["path"]] = item
        check(item)
    corrected_evidence(plan_path, REPLAY_SHA)
    if verify(core, launcher, runtime_path) != runtime:
        raise ValueError("Runtime changed during candidate admission")
    result = {"status": "cpm_candidates_admitted_unscored", "arm": ARMS[index], "index": index,
        "source": record(__file__), "scheduler": scheduler, "accounting": accounting,
        "preparation": preparation, "candidate_arm": arm, "context": context, "numeric_recheck": numeric,
        "replay_admission": replay["report_record"], "content_audit": content,
        "verification": {"genes": len(universe), "seed_groups": len(seeds), "candidate_groups": len(groups),
                         "reconstructed_merges": len(events)}, "checked_records": list(unique.values()),
        "accuracy_evaluated": False, "publication_ready": False,
        "limitations": ["Content/merge consistency, not independent rescoring of candidate search support.",
                        "Inferred phylogeny, native pair validation and scoring remain required.",
                        "Shared-host preparation times are not controlled efficiency evidence."]}
    with destination.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--index", type=int, choices=range(2), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.index, args.output.resolve())
