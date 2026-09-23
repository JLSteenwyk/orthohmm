"""Independently admit prespecified CPM replay artifacts, not their accuracy."""

import argparse
import csv
import io
import json
import math
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_corrected_replay import STAGES, FILENAMES
from benchmark_tools.audit_corrected_replay_stages import audit
from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.checked_replay_payload_worker import corrected_evidence
from benchmark_tools.cpm_replay_context import evidence, REPLAY_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_qfo_corrected_replay import validate_completion
from benchmark_tools.run_qfo_cpm_variant import ARMS, control_evidence
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_qfo_replay_launcher import verify

JOB = "22059"
EXECUTOR = "b4acb996e8accdcb604a66c95db2f5d8b68ab599"


def completed_task(accounting, index):
    if type(index) is not int or index not in range(2):
        raise ValueError("Unknown CPM variant index")
    rows = [r for r in csv.DictReader(io.StringIO(accounting), delimiter="|")
            if r["JobID"] == f"{JOB}_{index}"]
    if len(rows) != 1 or tuple(rows[0][k] for k in ("State", "ExitCode", "NodeList", "AllocCPUS")) != (
            "COMPLETED", "0:0", "bizon", "32"):
        raise ValueError("Require completed 32-CPU CPM task")
    return rows[0]


def validate_parent(parent, worker, replay, context, scheduler, root, executor, index, helpers):
    if (parent["job_id"] != scheduler["JobIDRaw"] or parent["executor_commit"] != EXECUTOR
            or parent["arm"] != ARMS[index] or type(parent["index"]) is not int or parent["index"] != index
            or parent["status"] != "cpm_variant_replay_complete_pending_admission"
            or type(parent["exit_code"]) is not int or parent["exit_code"] != 0
            or parent["accuracy_evaluated"] is not False or parent["publication_ready"] is not False
            or worker["accuracy_evaluated"] is not False):
        raise ValueError("Variant completion/identity differs")
    if (parent["context"] != context or worker["context"] != context
            or parent["source"] != record(executor / "benchmark_tools/run_qfo_cpm_variant.py")
            or worker["source"] != record(executor / "benchmark_tools/run_qfo_cpm_control.py")
            or parent["helpers"] != helpers):
        raise ValueError("Variant context/source binding differs")
    times = [parent["started_epoch"], parent["finished_epoch"]]
    if (any(type(t) not in (int, float) or not math.isfinite(t) for t in times)
            or not 0 < times[0] <= times[1]):
        raise ValueError("Invalid variant timestamps")
    expected = [sys.executable, "-B", str(executor / "benchmark_tools/run_qfo_cpm_variant.py"),
                "--root", str(root), "--index", str(index), "--worker"]
    if (parent["worker_command"] != expected or replay["command"] != context["native_command"]
            or replay["cwd"] != context["cwd"] or replay["source"] != worker["replay_source"]):
        raise ValueError("Variant command/source/cwd differs")
    parameters = {"accuracy_profile": "high_sensitivity", "cpm_resolution": context["resolution"],
        "profile_expansion": True, "profile_iterations": 1, "jackknife_profile_thresholds": False,
        "jackknife_single_copy_profiles": False, "profile_min_species": 1, "matrix": "BLOSUM62", "leiden_seed": 4}
    if replay["parameters"] != parameters or context["expected_stages"] != STAGES:
        raise ValueError("Variant scientific settings differ")
    validate_completion(worker, replay, STAGES)
    if type(replay["counts"]["species"]) is not int or replay["counts"]["species"] != 78:
        raise ValueError("Wrong variant species universe")
    for stage, filename in zip(replay["stages"], FILENAMES):
        if (stage["output"]["path"] != str(Path(context["output_root"]) / "replay" / filename)
                or "official_orthobench" in stage):
            raise ValueError("Wrong variant output path or unexpected scoring")


def partitions(replay, control, universe):
    if [r["label"] for r in control] != STAGES or [r["label"] for r in replay["stages"]] != STAGES:
        raise ValueError("Wrong comparison stages")
    coverage, comparisons = [], []
    for stage, prior in zip(replay["stages"], control):
        check(stage["output"])
        check(prior["output"])
        groups = read_partition(Path(stage["output"]["path"]), universe)
        baseline = read_partition(Path(prior["output"]["path"]), universe)
        if type(stage["clusters"]) is not int or stage["clusters"] != len(groups):
            raise ValueError("Variant group count differs")
        left, right = {frozenset(g) for g in groups}, {frozenset(g) for g in baseline}
        coverage.append({"label": stage["label"], "groups": len(groups), "output": stage["output"]})
        comparisons.append({"label": stage["label"], "control": prior["output"],
            "partition_equal": left == right, "variant_only_groups": len(left - right),
            "control_only_groups": len(right - left)})
    return coverage, comparisons


def admit(root, index, destination):
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    accounting = subprocess.check_output(["sacct", "-j", JOB, "--parsable2",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = completed_task(accounting, index)
    local_helpers = [record(p) for p in sorted(Path(__file__).parent.glob("*.py"))]
    plan_path = root / "benchmarks/work/qfo_corrected_replay_commands_20260918.json"
    plan, plan_record, admission_record, names_record = corrected_evidence(plan_path, REPLAY_SHA)
    context = evidence(root, plan, plan_record, ARMS[index])
    authorization = control_evidence(root, evidence(root, plan, plan_record, "control"))
    directory = Path(context["output_root"])
    executor = root / "benchmarks/work/publication_qfo_cpm_variant_v2"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR:
        raise ValueError("Variant executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    items = [record(directory / name) for name in ("results.json", "checked_worker.json", "replay.json", "fresh_control_admission.json")]
    parent, worker, replay, fresh = [read_frozen(Path(r["path"]), r["sha256"]) for r in items]
    if (parent["worker"] != items[1] or parent["replay"] != items[2]
            or parent["fresh_control_admission"] != items[3] or parent["authorization"] != authorization
            or fresh != authorization["report"]):
        raise ValueError("Variant control authorization/output binding differs")
    expected_admission = [sys.executable, "-B", str(Path(authorization["executor"]) / "benchmark_tools/admit_qfo_cpm_control.py"),
                          "--root", str(root), "--output", str(directory / "fresh_control_admission.json")]
    if parent["admission_command"] != expected_admission:
        raise ValueError("Independent control validation command differs")
    helpers = [record(p) for p in sorted((executor / "benchmark_tools").glob("*.py"))]
    validate_parent(parent, worker, replay, context, scheduler, root, executor, index, helpers)
    launcher = Path(context["cwd"])
    scientific_source = record(launcher / "benchmark_tools/replay_high_sensitivity.py")
    numeric = read_frozen(Path(admission_record["path"]), admission_record["sha256"])["content"]["numeric_checkpoint"]
    if (replay["source"] != scientific_source
            or any(replay["input"][k] != numeric[k] for k in ("status", "summary", "manifest"))
            or replay["input"]["accuracy_evaluated"] is not False
            or replay["input"]["auditor"] != record(launcher / "benchmark_tools/audit_accuracy_checkpoint.py")):
        raise ValueError("Variant scientific input differs")
    core = root / "benchmarks/work/publication_method_native_v2"
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    runtime = verify(core, launcher, runtime_path)
    if not runtime == plan["runtime"] == parent["runtime_before"] == parent["runtime_after"]:
        raise ValueError("Variant runtime differs")
    stages = audit(root, executor, directory, worker, replay, plan_record, cpm_arm=ARMS[index])
    names = Path(names_record["path"]).read_text().splitlines()
    if len(names) != 984137 or len(set(names)) != len(names):
        raise ValueError("Wrong variant gene universe")
    coverage, comparisons = partitions(replay, authorization["report"]["stage_comparisons"], set(names))
    if coverage != parent["coverage"]:
        raise ValueError("Variant coverage differs from independent reconstruction")
    records = [*items, plan_record, admission_record, names_record, scientific_source, replay["input"]["auditor"],
        *helpers, *local_helpers, *authorization["checked_records"], *context["checked_records"],
        *stages["checked_records"], *[r["output"] for r in coverage],
        *[record(directory / name) for name in ("control_validation.log", "replay.log", "time.txt")]]
    unique = {}
    for item in records:
        if item["path"] in unique and unique[item["path"]] != item:
            raise ValueError("Conflicting variant provenance")
        unique[item["path"]] = item
        check(item)
    corrected_evidence(plan_path, REPLAY_SHA)
    evidence(root, plan, plan_record, ARMS[index])
    if verify(core, launcher, runtime_path) != runtime:
        raise ValueError("Runtime changed during variant admission")
    result = {"status": "cpm_variant_replay_admitted_unscored", "arm": ARMS[index], "index": index,
        "scheduler": scheduler, "accounting": accounting, "context": context, "source": record(__file__),
        "source_report": items[0], "coverage": coverage, "control_comparison": comparisons,
        "clustering": stages["clustering"], "control_admission": authorization["report_record"],
        "checked_records": list(unique.values()), "accuracy_evaluated": False, "publication_ready": False,
        "limitations": ["Replay artifacts only; candidate construction, phylogeny and scoring remain required.",
                        "Different partitions are expected under changed CPM and are not accuracy evidence.",
                        "Cached shared-host replay is not controlled end-to-end timing."]}
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
