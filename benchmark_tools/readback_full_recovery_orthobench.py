"""Admit and score only the frozen full canonical recovery OrthoBench run."""

import argparse
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

from benchmark_tools.audit_installed_orthobench import baseline_inputs, read_root_hogs, verify_input_inventory
from benchmark_tools.audit_publication_pipeline import audit as scientific_audit
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.readback_canonical_ob_phylogeny import comparison, SCORER_SHA
from benchmark_tools.run_full_recovery_orthobench import validate_plan, INSTALL_SHA
from benchmark_tools.run_publication_pipeline import arguments, save
from benchmark_tools.score_orthobench_partition import score_partition
from benchmark_tools.verify_ygob_validation import require_completed_job

PLAN_SHA = "57878c0712b3216376ea1d431b1d8e9e8e3ae8eb58e3bd3069e15cf7dddb7005"
PROTOCOL_SHA = "d03a695d670a42d15c265a84eae46881c4908317609bd2fc5b5fa3d3583de8ac"
JOB = 22326


def verify_execution(plan, plan_record, submission, execution_started, execution, native_started,
                     native_complete, native_started_record, native_complete_record, job):
    directory = Path(plan["directory"])
    if (submission["job_id"] != str(job) or submission["plan"] != plan_record
            or execution_started != dict(plan=plan_record, job_id=str(job))
            or execution["status"] != "native_complete_pending_independent_readback"
            or execution["returncode"] != 0 or execution["job_id"] != str(job)
            or execution["plan"] != plan_record or execution["native_complete"] != native_complete_record):
        raise ValueError("Native execution is not bound to the requested completed attempt")
    expected_command = ["/usr/bin/time", "-v", "-o", str(directory / "time.txt"), *plan["command"]]
    if execution["command"] != expected_command:
        raise ValueError("Changed native command")
    if (native_complete["status"] != "native_complete_pending_scientific_readback"
            or native_complete["started"] != native_started_record
            or native_complete["accuracy_evaluated"] is not False
            or native_started["command"] != plan["command"][2:]
            or native_started["executable"] != plan["command"][0]
            or native_started["attempts"] != 1 or native_started["checkpoint_reuse"] is not False
            or native_started["production_default_changed"] is not False):
        raise ValueError("Wrong native interpreter, command or scope")
    command = plan["command"]
    value = lambda key: command[command.index(key) + 1]
    expected_args = arguments(SimpleNamespace(input=Path(value("--input")),
        output=Path(value("--output")), cpu=int(value("--cpu")),
        aligner=Path(value("--aligner")), tree_builder=Path(value("--tree-builder"))))
    if native_started["arguments"] != expected_args:
        raise ValueError("Frozen scientific arguments differ")
    keys = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "PYTHONHASHSEED", "MAFFT_BINARIES", "PATH")
    if native_started["environment"] != {k: plan["environment"][k] for k in keys}:
        raise ValueError("Native environment differs")


def admit(repo, directory, job):
    if job != JOB:
        raise ValueError("This readback is bound to the prespecified full recovery job")
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, job)
    if scheduler["AllocCPUS"] != "32":
        raise ValueError("Wrong CPU allocation")
    path = directory / "plan.json"
    plan = validate_plan(path, PLAN_SHA)
    if Path(plan["directory"]) != directory or Path(plan["repo"]) != repo:
        raise ValueError("Wrong run directory or repository")
    protocol_path = repo / "benchmark_tools/results/full_recovery_orthobench_protocol_20260926.json"
    if record(protocol_path)["sha256"] != PROTOCOL_SHA:
        raise ValueError("Changed prespecified protocol")
    protocol = json.loads(protocol_path.read_text())
    if protocol["plan"] != record(path):
        raise ValueError("Protocol plan differs")
    files = ("submission.json", "execution_started.json", "execution.json",
             "native/started.json", "native/complete.json")
    submission, started, executed, native, complete = [json.loads((directory / n).read_text()) for n in files]
    if submission["protocol"] != record(protocol_path):
        raise ValueError("Submission protocol differs")
    if (directory / "native/failure.json").exists():
        raise ValueError("Native failure retained")
    verify_execution(plan, record(path), submission, started, executed, native, complete,
                     record(directory / "native/started.json"), record(directory / "native/complete.json"), job)
    checked = [record(directory / n) for n in files] + [record(protocol_path), *protocol["references"]]
    pinned = {r["path"]: r for r in plan["checked_records"]}
    for item in [native["source"], native["policy"], *native["inputs"],
                 *native["scientific_sources"], *native["tools"]]:
        if pinned.get(item["path"]) != item:
            raise ValueError("Native source/input/tool was not pinned")
    for name in ("install_preparation.json", "install_before.json", "install_after.json"):
        item = record(directory / name)
        if item["sha256"] != INSTALL_SHA:
            raise ValueError("Recovery installation changed")
        checked.append(item)
    if executed["logs"] != [record(directory / n) for n in ("native.log", "time.txt")]:
        raise ValueError("Missing execution logs")
    checked.extend(executed["logs"])
    for item in checked:
        check(item)
    return dict(status="full_recovery_native_bound_pending_scientific_readback", plan=record(path),
        job_id=job, accounting=accounting, scheduler=scheduler, checked_records=checked,
        source=record(__file__), scientific_scores_admitted=False)


def readback(repo, directory, job, output):
    if output.exists():
        raise FileExistsError(output)
    admitted = admit(repo, directory, job)
    plan = json.loads((directory / "plan.json").read_text())
    scorer = record(repo / "benchmark_tools/score_orthobench_partition.py")
    if scorer["sha256"] != SCORER_SHA:
        raise ValueError("Changed frozen scorer")
    baseline, refs, uncertain, baseline_records = baseline_inputs(repo)
    input_dir = Path(plan["command"][plan["command"].index("--input") + 1])
    universe, _ = verify_input_inventory(input_dir, plan)
    output.mkdir(parents=True)
    save(output / "admission.json", admitted)
    scientific = scientific_audit(directory / "native", output / "scientific")
    current_path = directory / "native/inference/orthohmm_phylogeny/orthohmm_root_hogs.tsv"
    paths = dict(historical=Path(plan["baseline_partition"]["path"]), current=current_path)
    partitions = {key: read_root_hogs(path, universe) for key, path in paths.items()}
    scores = {key: score_partition(partition, refs, uncertain) for key, partition in partitions.items()}
    if scores["historical"] != baseline["scores"]["p1_c1_r1"]:
        raise ValueError("Historical score did not reproduce exactly")
    for item in baseline_records:
        check(item)
    if admit(repo, directory, job) != admitted:
        raise ValueError("Admission changed during scientific readback")
    result = dict(status="full_recovery_orthobench_scientific_readback_complete", job_id=job,
        plan=admitted["plan"], genes=len(universe), scores=scores,
        comparison=comparison(partitions["historical"], partitions["current"], scores["historical"], scores["current"]),
        partitions={key: record(path) for key, path in paths.items()},
        readbacks=[record(output / "admission.json"), record(output / "scientific/result.json")],
        summary=scientific["summary"], scorer=scorer, source=record(__file__),
        historical_scores_replaced=False, publication_ready=False,
        limitations=["Development-exposed same-host reproducibility experiment, not independent generalization.",
                     "Shared-host resources do not establish controlled comparative speed."])
    save(output / "result.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "directory", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--job", type=int, required=True)
    args = parser.parse_args()
    readback(args.repo.resolve(), args.directory.resolve(), args.job, args.output.resolve())
