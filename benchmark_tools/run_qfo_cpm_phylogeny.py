"""Run admitted CPM-specific candidates through independently inferred phylogeny."""

import argparse
import csv
import io
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.checked_replay_payload_worker import corrected_evidence
from benchmark_tools.cpm_replay_context import evidence, REPLAY_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_ob_candidate_neighborhood import variant_cell
from benchmark_tools.run_qfo_cpm_variant import ARMS
from benchmark_tools.run_qfo_parameter_phylogeny import verify_baseline
from benchmark_tools.run_qfo_factorial_cell import native_command
from benchmark_tools.run_simulation_methods import read_frozen, execution_environment, execute

ADMISSION_JOB = "22086"
ADMISSION_COMMIT = "497c54f1de8010482c0aaabd92c2c4e68dff43da"
ADMISSION_SHA = "afc76384bf31e720f9491659978a403fe38d084fc552f4087a59034e5bd3bde8"


def select_arm(report, index, context):
    if (type(index) is not int or index not in range(2)
            or report["status"] != "cpm_candidates_admitted_unscored"
            or type(report["index"]) is not int or report["index"] != index or report["arm"] != ARMS[index]
            or report["context"] != context or report["accuracy_evaluated"] is not False
            or report["publication_ready"] is not False):
        raise ValueError("Require matching unscored CPM candidate admission")
    arm = report["candidate_arm"]
    expected = Path(context["output_root"]) / "replay/orthogroups_profiles_refined.txt"
    if arm["seed_partition"]["path"] != str(expected) or arm["candidate_expansion"] is not True:
        raise ValueError("Wrong CPM-specific seed or expansion mode")
    return {"label": ARMS[index], "partition": arm["candidate_partition"],
            "constraints": arm["membership_constraints"], "seed_partition": arm["seed_partition"]}


def verify_sources(root, index):
    if type(index) is not int or index not in range(2):
        raise ValueError("Unknown CPM phylogeny index")
    accounting = subprocess.check_output(["sacct", "-j", ADMISSION_JOB, "--parsable2",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    rows = [r for r in csv.DictReader(io.StringIO(accounting), delimiter="|") if r["JobID"] == f"{ADMISSION_JOB}_{index}"]
    if len(rows) != 1 or tuple(rows[0][k] for k in ("State", "ExitCode", "NodeList", "AllocCPUS")) != (
            "COMPLETED", "0:0", "bizon", "2"):
        raise ValueError("Require completed CPM candidate admission")
    executor = root / "benchmarks/work/publication_qfo_cpm_candidates_admission_v3"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != ADMISSION_COMMIT:
        raise ValueError("Candidate admission executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    source = record(executor / "benchmark_tools/admit_qfo_cpm_candidates.py")
    if source["sha256"] != ADMISSION_SHA:
        raise ValueError("Candidate admission source changed")
    path = root / f"benchmarks/work/qfo_cpm_candidates_admission_{ADMISSION_JOB}_{index}.json"
    item = record(path)
    report = read_frozen(path, item["sha256"])
    plan, plan_record, _, _ = corrected_evidence(root / "benchmarks/work/qfo_corrected_replay_commands_20260918.json", REPLAY_SHA)
    context = evidence(root, plan, plan_record, ARMS[index])
    arm = select_arm(report, index, context)
    if report["source"] != source:
        raise ValueError("Wrong candidate admission source binding")
    records = [source, item, *report["checked_records"], *context["checked_records"],
               arm["partition"], arm["constraints"], arm["seed_partition"]]
    for entry in records:
        check(entry)
    manifest, original, launcher, prepared, environment, baseline_records = verify_baseline(root)
    return {"arm": arm, "context": context, "admission": report, "admission_record": item,
        "admission_executor": str(executor), "scheduler": rows[0], "accounting": accounting,
        "manifest": manifest, "original": original, "launcher": str(launcher), "prepared": str(prepared),
        "environment": environment, "checked_records": [*records, *baseline_records]}


def run(root, index):
    if (type(index) is not int or index not in range(2) or not os.environ.get("SLURM_JOB_ID")
            or os.environ.get("SLURM_CPUS_PER_TASK") != "32" or os.environ.get("SLURM_JOB_NODELIST") != "bizon"
            or os.environ.get("SLURM_ARRAY_TASK_ID") != str(index)):
        raise ValueError("Require scheduled 32-CPU CPM phylogeny task")
    output = root / "benchmarks/results/qfo_cpm_phylogeny_v1" / ARMS[index]
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    verified = verify_sources(root, index)
    launcher, prepared = Path(verified["launcher"]), Path(verified["prepared"])
    cell = variant_cell(verified["original"], verified["arm"], output)
    argv, equivalence = native_command(cell, launcher, prepared)
    env, resolved = execution_environment(verified["environment"])
    env.update(verified["manifest"]["environment_overrides"], PYTHONPATH=str(launcher))
    helpers = [record(p) for p in sorted(Path(__file__).parent.glob("*.py"))]
    provenance = {"source": record(__file__), "helpers": helpers, "verified": verified, "cell": cell,
        "executed_argv": argv, "launcher_source_equivalence": equivalence, "resolved_tools": resolved,
        "cwd": str(launcher), "job_id": os.environ["SLURM_JOB_ID"], "array_task_id": str(index),
        "scope": "CPM-specific seed/candidate families; independently inferred phylogeny; no accuracy scoring"}
    output.mkdir(parents=True, exist_ok=False)
    (output / "preflight.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    postflight = {"status": "running", "accuracy_evaluated": False, "native_outputs_validated": False}
    cwd = Path.cwd()
    try:
        fresh = output / "fresh_candidate_admission.json"
        command = [sys.executable, "-B", str(Path(verified["admission_executor"]) / "benchmark_tools/admit_qfo_cpm_candidates.py"),
                   "--root", str(root), "--index", str(index), "--output", str(fresh)]
        postflight["admission_command"] = command
        with (output / "admission.log").open("x") as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
        fresh_record = record(fresh)
        if read_frozen(fresh, fresh_record["sha256"]) != verified["admission"]:
            raise ValueError("Fresh candidate admission disagrees")
        postflight["fresh_candidate_admission"] = fresh_record
        method = {"argv": argv, "output": str(output / "output"), "metrics": str(output / "metrics.json")}
        inputs = {"status": "ready", "inputs": [{**r, "absolute_path": r["path"]} for r in verified["manifest"]["input_fastas"]]}
        os.chdir(launcher)
        execution = execute({"label": cell["label"], "methods": {cell["label"]: method}}, [cell["label"]],
                            env, output / "execution", inputs, provenance)
        os.chdir(cwd)
        if execution["failed_methods"]:
            raise RuntimeError("CPM phylogeny failed; preserve outputs without retry")
        if verify_sources(root, index) != verified:
            raise ValueError("CPM phylogeny inputs changed during execution")
        for item in [provenance["source"], fresh_record, *helpers, *verified["checked_records"]]:
            check(item)
        for pair in equivalence:
            check(pair["prepared"])
            check(pair["executed"])
        postflight.update(status="complete_pending_native_validation", cell=cell)
    except BaseException as error:
        postflight.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        os.chdir(cwd)
        (output / "postflight.json").write_text(json.dumps(postflight, indent=2, sort_keys=True) + "\n")
    return postflight


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--index", type=int, choices=range(2), required=True)
    args = parser.parse_args()
    run(args.root.resolve(), args.index)
