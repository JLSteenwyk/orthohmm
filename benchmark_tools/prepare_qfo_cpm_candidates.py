"""Build fixed-profile candidates from each independently admitted CPM seed."""

import argparse
import csv
import importlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, controlled_expansion, record
from benchmark_tools.prepare_qfo_candidate_neighborhood import require_environment
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_qfo_cpm_variant import ARMS

ADMISSION_JOB = "21962"
ADMISSION_COMMIT = "cbd59182a59b151933faebd9c41d7803a38f0726"
ADMISSION_SHA = "96f4129d7ee948fc077ee67e6e4673c61acbb8f56e17cb17777239daa19593d0"
BASELINE_SHA = "d8385c50426e690afd6d32f3c5302e678de6977c9841451d013201b0f75b564a"


def replay_evidence(root, index, context):
    if type(index) is not int or index not in range(2):
        raise ValueError("Unknown CPM index")
    accounting = subprocess.check_output(["sacct", "-j", ADMISSION_JOB, "--parsable2",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    rows = [r for r in csv.DictReader(io.StringIO(accounting), delimiter="|")
            if r["JobID"] == f"{ADMISSION_JOB}_{index}"]
    if len(rows) != 1 or tuple(rows[0][k] for k in ("State", "ExitCode", "NodeList", "AllocCPUS")) != (
            "COMPLETED", "0:0", "bizon", "2"):
        raise ValueError("Require completed CPM replay admission")
    executor = root / "benchmarks/work/publication_qfo_cpm_variant_admission_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != ADMISSION_COMMIT:
        raise ValueError("Changed replay admission executor")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    source = record(executor / "benchmark_tools/admit_qfo_cpm_variant.py")
    if source["sha256"] != ADMISSION_SHA:
        raise ValueError("Changed replay admission source")
    item = record(root / f"benchmarks/work/qfo_cpm_variant_admission_{ADMISSION_JOB}_{index}.json")
    report = read_frozen(Path(item["path"]), item["sha256"])
    if (report["status"] != "cpm_variant_replay_admitted_unscored" or type(report["index"]) is not int or report["index"] != index
            or report["arm"] != ARMS[index] or report["context"] != context or report["source"] != source
            or report["accuracy_evaluated"] is not False or report["publication_ready"] is not False
            or [r["label"] for r in report["coverage"]] != context["expected_stages"]):
        raise ValueError("Wrong CPM replay admission")
    seed = report["coverage"][-1]["output"]
    if seed["path"] != str(Path(context["output_root"]) / "replay/orthogroups_profiles_refined.txt"):
        raise ValueError("CPM candidate seed is not the admitted refined variant")
    records = [source, item, seed, *report["checked_records"]]
    for entry in records:
        check(entry)
    return {"scheduler": rows[0], "accounting": accounting, "executor": str(executor),
            "report": report, "report_record": item, "seed_partition": seed, "checked_records": records}


def build(engine, parameters, seed, names, species, hits, directory, auditor):
    working = directory / "orthohmm_working_res"
    working.mkdir(parents=True, exist_ok=False)
    partition = working / "orthohmm_edges_clustered.txt"
    check(seed)
    shutil.copyfile(seed["path"], partition)
    started = time.monotonic()
    expansion = controlled_expansion(engine, parameters, "control", (str(directory), names, species, hits))
    elapsed = time.monotonic() - started
    arm = {"seed_partition": seed, "candidate_expansion": True, "candidate_partition": record(partition),
        "membership_constraints": record(working / "phylogeny_candidate_merges.json"),
        "expansion": expansion["engine_fixed_profile_report"],
        "output_files": [record(p) for p in sorted(directory.rglob("*")) if p.is_file()]}
    content = auditor(arm, seed, directory, set(names), True)
    check(seed)
    return {"candidate_arm": arm, "content_audit": content, "candidate_parameter_control": expansion,
            "incremental_seconds": elapsed}


def prepare(root, index):
    require_environment(os.environ)
    if type(index) is not int or index not in range(2) or os.environ.get("SLURM_ARRAY_TASK_ID") != str(index):
        raise ValueError("Require matching CPM preparation array index")
    from benchmark_tools.checked_replay_payload_worker import corrected_evidence
    from benchmark_tools.cpm_replay_context import evidence, REPLAY_SHA
    from benchmark_tools.verify_qfo_replay_launcher import verify
    from benchmark_tools.audit_accuracy_checkpoint import audit
    from benchmark_tools.audit_candidate_arm import audit as audit_arm
    plan_path = root / "benchmarks/work/qfo_corrected_replay_commands_20260918.json"
    plan, plan_record, _, _ = corrected_evidence(plan_path, REPLAY_SHA)
    context = evidence(root, plan, plan_record, ARMS[index])
    admitted = replay_evidence(root, index, context)
    baseline_path = root / "benchmarks/results/qfo_corrected_factorial_v1/manifest.json"
    baseline = read_frozen(baseline_path, BASELINE_SHA)
    parameters = baseline["candidate_arms"]["p1_c1"]["expansion"]["parameters"]
    output = root / "benchmarks/results/qfo_cpm_candidates_v1" / ARMS[index]
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    core = root / "benchmarks/work/publication_method_native_v2"
    launcher = Path(context["cwd"])
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    runtime = verify(core, launcher, runtime_path)
    if not runtime == plan["runtime"] == baseline["runtime_before"] == baseline["runtime_after"]:
        raise ValueError("Candidate runtime differs")
    sys.path.insert(0, str(launcher))
    try:
        engine = importlib.import_module("orthohmm.orthohmm")
        accuracy = importlib.import_module("orthohmm.accuracy")
    finally:
        sys.path.remove(str(launcher))
    for module in (engine, accuracy):
        if Path(module.__file__).resolve().parent != launcher / "orthohmm":
            raise ValueError("Wrong frozen scientific import")
    checkpoint_record = plan["checkpoint_manifest"]
    checkpoint = Path(checkpoint_record["path"]).parent
    numeric = audit(checkpoint, checkpoint_record["sha256"])
    if numeric["summary"] != baseline["numeric_checkpoint"]["summary"] or numeric["summary"]["species"] != 78:
        raise ValueError("Changed candidate hit checkpoint")
    names, species, queries, targets, scores = accuracy.load_accuracy_checkpoint(checkpoint, verify=False)
    if len(names) != 984137 or len(set(names)) != len(names):
        raise ValueError("Wrong candidate gene universe")
    helpers = [record(p) for p in sorted(Path(__file__).parent.glob("*.py"))]
    inputs = [plan_record, record(baseline_path), record(runtime_path), *context["checked_records"],
              *admitted["checked_records"], *helpers]
    for item in inputs:
        check(item)
    output.mkdir(parents=True, exist_ok=False)
    report = {"status": "preparing_unscored", "arm": ARMS[index], "index": index, "context": context,
        "source": record(__file__), "inputs": inputs, "replay_admission": admitted,
        "job_id": os.environ["SLURM_JOB_ID"], "array_task_id": str(index),
        "executor_commit": subprocess.check_output(["git", "-C", str(Path(__file__).resolve().parent.parent),
            "rev-parse", "HEAD"], text=True).strip(), "runtime_before": runtime, "numeric_checkpoint": numeric,
        "accuracy_evaluated": False, "publication_ready": False}
    try:
        fresh = output / "fresh_replay_admission.json"
        command = [sys.executable, "-B", str(Path(admitted["executor"]) / "benchmark_tools/admit_qfo_cpm_variant.py"),
                   "--root", str(root), "--index", str(index), "--output", str(fresh)]
        report["admission_command"] = command
        with (output / "admission.log").open("x") as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
        fresh_record = record(fresh)
        if read_frozen(fresh, fresh_record["sha256"]) != admitted["report"]:
            raise ValueError("Fresh replay admission disagrees")
        report["fresh_replay_admission"] = fresh_record
        report.update(build(engine, parameters, admitted["seed_partition"], names, species,
                            (queries, targets, scores), output / "candidate", audit_arm))
        if audit(checkpoint, checkpoint_record["sha256"]) != numeric:
            raise ValueError("Candidate checkpoint changed during preparation")
        report["runtime_after"] = verify(core, launcher, runtime_path)
        if report["runtime_after"] != runtime:
            raise ValueError("Candidate runtime changed during preparation")
        for item in [*inputs, fresh_record, *report["candidate_arm"]["output_files"]]:
            check(item)
        report.update(status="cpm_candidates_prepared_pending_admission", limitations=[
            "Candidate thresholds are unchanged; only the CPM-derived seed partition differs.",
            "Independent candidate admission, phylogeny and scoring remain required.",
            "Shared-host incremental preparation is not controlled end-to-end timing."])
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        (output / "manifest.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--index", type=int, choices=range(2), required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.index)
