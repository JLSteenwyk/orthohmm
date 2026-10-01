"""Independently admit source-corrected recovered CPM candidates, without scores."""

import argparse
import csv
import importlib
import importlib.util
import io
import json
import math
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

JOB = "22385"
COMMIT = "1a6dd9ef97a9b636842b1dbf579b9b8e30e90e47"
EXECUTOR = "benchmarks/work/qfo_cpm_helper_candidates_import_fixed_executor_20260930"
SUBMISSION = "benchmark_tools/results/qfo_cpm_helper_candidates_submission_22385.json"
SUBMISSION_SHA = "a5d5df529977646b8c3b2a9d51aba49422e6663d34f1806ad52cc8fbcf1cf219"
PREPARATION = "benchmarks/results/qfo_cpm_helper_recovered_candidates_v1/manifest.json"
SOURCE_SHA = "194680edee971dcf0b811e0191a04cac82717b5eb9ce444f42b34cfda0d4949a"
GATE_SHA = "01164124f97d06dc5752499c7616b05a94b6a0964c756e8cf6a09e84705c00c0"
PROTOCOL = "benchmark_tools/results/QFO_CPM_HELPER_CANDIDATE_ADMISSION_PROTOCOL_20260930.md"
GENES = 984137


def completed(accounting):
    rows = [row for row in csv.DictReader(io.StringIO(accounting), delimiter="|") if row["JobID"] == JOB]
    if len(rows) != 1 or tuple(rows[0][key] for key in
            ("JobIDRaw", "State", "ExitCode", "AllocCPUS", "ReqMem", "NodeList")) != (
                JOB, "COMPLETED", "0:0", "2", "64G", "bizon"):
        raise ValueError("Require completed source-corrected candidate job " + JOB)
    return rows[0]


def pinned_submission(root):
    submission_record = record(root / SUBMISSION)
    if submission_record["sha256"] != SUBMISSION_SHA:
        raise ValueError("Candidate submission identity changed")
    submission = json.loads(Path(submission_record["path"]).read_bytes())
    executor = root / EXECUTOR
    if (submission["status"] != "source_corrected_explicit_helper_seed_candidate_job_submitted"
            or submission["job_id"] != JOB or submission["executor"] != str(executor)
            or submission["executor_commit"] != COMMIT
            or submission["seed_handoff"] != "explicit_helper_runtime_seed_amendment"
            or any(submission[key] is not False for key in
                   ("candidate_admitted", "accuracy_evaluated", "publication_ready"))):
        raise ValueError("Wrong corrected candidate submission")
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != COMMIT:
        raise ValueError("Corrected candidate executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"],
                   check=True, capture_output=True)
    records = [submission_record, submission["previous_failed_attempt"], submission["seed_readback"],
               *submission["source_records"]]
    for item in records:
        check(item)
    source = record(executor / "benchmark_tools/prepare_recovered_cpm_candidates.py")
    gate = record(executor / "benchmark_tools/helper_recovered_cpm_seed_evidence.py")
    if source["sha256"] != SOURCE_SHA or gate["sha256"] != GATE_SHA or source not in records or gate not in records:
        raise ValueError("Wrong candidate source/gate bindings")
    return submission, submission_record, source, gate, executor, records


def replay_seed_evidence(root, executor, gate, report):
    if Path(gate["path"]) != executor / "benchmark_tools/helper_recovered_cpm_seed_evidence.py":
        raise ValueError("Seed evidence must be replayed from the verified candidate executor")
    check(gate)
    spec = importlib.util.spec_from_file_location("verified_candidate_seed_handoff", gate["path"])
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # Match __file__ provenance to the verified builder executor, not this new executor.
    selected = report["recovery_admission"]
    result = module.evidence(root, root / module.READBACK, module.READBACK_SHA, selected["protocol"]["sha256"])
    if result != selected:
        raise ValueError("Fresh explicit seed handoff differs from candidate preparation")
    return result


def validate_report(report, scheduler, context, parameters, source, builder, inputs):
    if (report["status"] != "recovered_cpm_candidates_prepared_pending_admission"
            or report["arm"] != "cpm_high" or report["job_id"] != scheduler["JobIDRaw"]
            or report["executor_commit"] != COMMIT or report["source"] != source or report["builder"] != builder
            or report["context"] != context or report["inputs"] != inputs
            or report["seed_handoff"] != "explicit_helper_runtime_seed_amendment"
            or report["accuracy_evaluated"] is not False or report["publication_ready"] is not False):
        raise ValueError("Corrected recovered candidate preparation differs")
    control = report["candidate_parameter_control"]
    if (control["label"] != "control" or control["delta"] != {} or control["applied_parameters"] != parameters
            or type(control["engine_calls"]) is not int or control["engine_calls"] != 1
            or control["engine_fixed_profile_report"]["parameters"] != parameters
            or report["candidate_arm"]["expansion"] != control["engine_fixed_profile_report"]):
        raise ValueError("Recovered candidates changed fixed expansion parameters")
    seconds = report["incremental_seconds"]
    if type(seconds) not in (int, float) or not math.isfinite(seconds) or seconds < 0:
        raise ValueError("Invalid candidate incremental time")


def scientific_imports(launcher):
    sys.path.insert(0, str(launcher))
    try:
        engine = importlib.import_module("orthohmm.orthohmm")
        accuracy = importlib.import_module("orthohmm.accuracy")
    finally:
        sys.path.remove(str(launcher))
    if any(Path(module.__file__).resolve().parent != launcher / "orthohmm" for module in (engine, accuracy)):
        raise ValueError("Wrong frozen candidate-admission scientific import")
    return accuracy


def admit(root, preparation_sha, protocol_sha, destination):
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    accounting = subprocess.check_output(["sacct", "-j", JOB, "-X", "--parsable2",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,AllocCPUS,ReqMem,NodeList"], text=True)
    scheduler = completed(accounting)
    protocol = record(root / PROTOCOL)
    if protocol["sha256"] != protocol_sha:
        raise ValueError("Unreviewed candidate admission protocol")
    submission, submission_record, source, gate, executor, submission_records = pinned_submission(root)
    preparation = record(root / PREPARATION)
    if preparation["sha256"] != preparation_sha:
        raise ValueError("Candidate preparation hash differs")
    report = json.loads(Path(preparation["path"]).read_bytes())
    admitted = replay_seed_evidence(root, executor, gate, report)
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    decoder = scientific_imports(launcher)
    # Scientific imports must be selected before any transitive analysis-helper imports.
    from benchmark_tools.audit_accuracy_checkpoint import audit as audit_numeric
    from benchmark_tools.audit_candidate_arm import audit as audit_arm
    from benchmark_tools.checked_replay_payload_worker import corrected_evidence
    from benchmark_tools.cpm_replay_context import evidence, REPLAY_SHA
    from benchmark_tools.prepare_qfo_cpm_candidates import BASELINE_SHA
    from benchmark_tools.trace_ob_families import partition, validate_merge_reconstruction
    from benchmark_tools.verify_qfo_replay_launcher import verify
    from benchmark_tools.run_simulation_methods import read_frozen

    plan, plan_record, _, names_record = corrected_evidence(
        root / "benchmarks/work/qfo_corrected_replay_commands_20260918.json", REPLAY_SHA)
    context = evidence(root, plan, plan_record, "cpm_high")
    baseline_path = root / "benchmarks/results/qfo_corrected_factorial_v1/manifest.json"
    baseline = read_frozen(baseline_path, BASELINE_SHA)
    parameters = baseline["candidate_arms"]["p1_c1"]["expansion"]["parameters"]
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    builder = record(executor / "benchmark_tools/prepare_qfo_cpm_candidates.py")
    helpers = [record(path) for path in sorted((executor / "benchmark_tools").glob("*.py"))]
    inputs = [plan_record, record(baseline_path), record(runtime_path), builder,
              *context["checked_records"], *admitted["checked_records"], *helpers]
    validate_report(report, scheduler, context, parameters, source, builder, inputs)
    core = root / "benchmarks/work/publication_method_native_v2"
    runtime = verify(core, launcher, runtime_path)
    if not runtime == plan["runtime"] == baseline["runtime_before"] == baseline["runtime_after"] == report["runtime_before"] == report["runtime_after"]:
        raise ValueError("Recovered candidate admission runtime differs")
    checkpoint = plan["checkpoint_manifest"]
    numeric = audit_numeric(Path(checkpoint["path"]).parent, checkpoint["sha256"])
    prior = report["numeric_checkpoint"]
    if (any(numeric[key] != prior[key] for key in ("status", "manifest", "summary", "accuracy_evaluated"))
            or prior["auditor"] != record(executor / "benchmark_tools/audit_accuracy_checkpoint.py")
            or numeric["summary"] != baseline["numeric_checkpoint"]["summary"]):
        raise ValueError("Recovered candidate numeric checkpoint differs")
    names = Path(names_record["path"]).read_text().splitlines()
    universe = set(names)
    if len(names) != GENES or len(universe) != GENES:
        raise ValueError("Wrong recovered candidate universe")
    directory = (root / PREPARATION).parent
    arm, seed = report["candidate_arm"], admitted["seed_partition"]
    content = audit_arm(arm, seed, directory / "candidate", universe, True)
    if content != report["content_audit"]:
        raise ValueError("Independent recovered candidate content differs")
    seeds, _ = partition(Path(seed["path"]), "plain", universe)
    groups, _ = partition(Path(arm["candidate_partition"]["path"]), "plain", universe)
    events = json.loads(Path(arm["membership_constraints"]["path"]).read_bytes())
    validate_merge_reconstruction(events, seeds, groups, set())
    records = [record(__file__), protocol, preparation, names_record, source, gate, record(decoder.__file__),
        *submission_records, *inputs, *[record(path) for path in sorted(Path(__file__).parent.glob("*.py"))],
        *arm["output_files"], *content["checked_records"],
        record(root / f"benchmarks/work/qfo_cpm_helper_candidates_{JOB}.log"),
        record(root / f"benchmarks/work/qfo_cpm_helper_candidates_{JOB}.time.txt")]
    unique = {}
    for item in records:
        if item["path"] in unique and unique[item["path"]] != item:
            raise ValueError("Conflicting recovered candidate provenance")
        unique[item["path"]] = item
        check(item)
    if replay_seed_evidence(root, executor, gate, report) != admitted:
        raise ValueError("Seed evidence changed during admission")
    corrected_evidence(root / "benchmarks/work/qfo_corrected_replay_commands_20260918.json", REPLAY_SHA)
    if verify(core, launcher, runtime_path) != runtime or completed(subprocess.check_output(
            ["sacct", "-j", JOB, "-X", "--parsable2",
             "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,AllocCPUS,ReqMem,NodeList"], text=True)) != scheduler:
        raise ValueError("Candidate admission runtime or completion state changed")
    result = dict(status="cpm_helper_recovered_candidates_admitted_unscored", candidate_admitted=True,
        arm="cpm_high", index=1, seed_handoff="explicit_helper_runtime_seed_amendment",
        source=record(__file__), protocol=protocol, scheduler=scheduler, accounting=accounting,
        submission=submission_record, preparation=preparation, candidate_arm=arm, context=context,
        numeric_recheck=numeric, seed_admission=admitted["report_record"], seed_readback=admitted["readback"],
        verification=dict(genes=len(universe), seed_groups=len(seeds), candidate_groups=len(groups),
                          reconstructed_merges=len(events)), content_audit=content,
        checked_records=list(unique.values()), accuracy_evaluated=False, downstream_admitted=False,
        publication_ready=False, limitations=[
            "Independent complete membership/merge reconstruction, not rescoring candidate search support.",
            "Original failures and explicit independent-refinement runtime amendment remain preserved.",
            "Inferred phylogeny, pair conversion/validation and scoring require separate admissions.",
            "Shared-host incremental construction/admission is not controlled comparative timing."])
    with destination.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--preparation-sha256", required=True)
    parser.add_argument("--protocol-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.preparation_sha256, args.protocol_sha256, args.output.resolve())
