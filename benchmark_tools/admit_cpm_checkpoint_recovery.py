"""Independently admit a recovered high-CPM seed, not downstream accuracy."""

import argparse
import csv
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_blast_recovery_batch import save_status
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.prepare_cpm_checkpoint_recovery import prepare
from benchmark_tools.validate_cpm_refinement_reconstruction import compare_partitions
from benchmark_tools.verify_blast_recovery_panel import unique_records

RUNNER_SHA = "a277ab01e7fbcbfaa15f7a63a092c78183baaabbbcc25cd23640f6a750eabbe2"


def completed(accounting, job):
    if not job or not job.isascii() or not job.isdigit():
        raise ValueError("Require explicit single recovery job ID")
    rows = [row for row in csv.DictReader(io.StringIO(accounting), delimiter="|") if row["JobID"] == job]
    if len(rows) != 1 or tuple(rows[0][k] for k in ("State", "ExitCode", "AllocCPUS", "ReqMem", "NodeList")) != (
            "COMPLETED", "0:0", "1", "64G", "bizon"):
        raise ValueError("Require completed one-CPU 64-GiB recovery")
    return rows[0]


def parent_contract(parent, job, commit, source, root, directory):
    if (parent["status"] != "cpm_checkpoint_recovered_pending_independent_admission"
            or parent["job_id"] != job or parent["executor_commit"] != commit or parent["source"] != source
            or parent["optimizer_attempted"] is not True
            or any(parent[key] is not False for key in ("accuracy_evaluated", "downstream_admitted", "publication_ready"))
            or parent["missing_original_statistics"] != ["profile_counters", "successful_stage_timings", "full_run_time"]):
        raise ValueError("Invalid recovered parent identity or status")
    phases = parent["phases"]
    if [row["mode"] for row in phases] != ["optimize", "refine", "repeat-refinement"]:
        raise ValueError("Repeated or missing recovery phase")
    for row in phases:
        expected = [sys.executable, "-B", source["path"], "--root", str(root),
                    "--output", str(directory), "--mode", row["mode"]]
        if (row["command"] != expected or row["status"] != "completed"
                or type(row["returncode"]) is not int or row["returncode"] != 0
                or type(row["wall_s"]) not in (int, float) or not math.isfinite(row["wall_s"]) or row["wall_s"] < 0
                or row["log"] != record(directory / f"{row['mode']}.log")):
            raise ValueError("Invalid completed recovery phase")
    original = root / "benchmarks/results/qfo_parameter_cpm_replay_v3/cpm_high/replay"
    expected_stages = [dict(label=label, origin=origin, output=record(path)) for label, origin, path in (
        ("multipass", "reused", original / "orthogroups_multipass.txt"),
        ("multipass_refined", "reused", original / "orthogroups_multipass_refined.txt"),
        ("strict_profiles", "recovered", directory / "orthogroups_profiles.txt"),
        ("strict_profiles_refined", "recovered", directory / "orthogroups_profiles_refined.txt"))]
    if parent["stages"] != expected_stages:
        raise ValueError("Reused/recovered stage identity differs")
    return expected_stages


def admit(root, job, commit, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    accounting = subprocess.check_output(["sacct", "-j", job, "--parsable2",
        "--format=JobID,State,ExitCode,Elapsed,AllocCPUS,ReqMem,NodeList"], text=True)
    scheduler = completed(accounting, job)
    executor = root / "benchmarks/work/cpm_checkpoint_recovery_v1_20260923"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != commit:
        raise ValueError("Changed recovery executor")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    source = record(executor / "benchmark_tools/run_cpm_checkpoint_recovery.py")
    if source["sha256"] != RUNNER_SHA:
        raise ValueError("Unreviewed recovery runner")
    directory = root / "benchmarks/results/qfo_cpm_checkpoint_recovery_v1"
    parent_record = record(directory / "status.json")
    parent = read_frozen(Path(parent_record["path"]), parent_record["sha256"])
    stages = parent_contract(parent, job, commit, source, root, directory)
    preflight_record = record(root / "benchmarks/results/qfo_cpm_checkpoint_recovery_preflight_v1/status.json")
    if parent["preflight"] != preflight_record:
        raise ValueError("Recovery preflight report differs")
    original_preflight = read_frozen(Path(preflight_record["path"]), preflight_record["sha256"])
    if original_preflight["status"] != "cpm_checkpoint_preflight_verified_unscored" or original_preflight["preflight_passed"] is not True:
        raise ValueError("Original recovery preflight incomplete")
    records = unique_records([record(__file__), source, parent_record, preflight_record,
        *parent["checked_records"], *original_preflight["checked_records"], *[r["output"] for r in stages],
        *[r["log"] for r in parent["phases"]], *parent["refinement_reports"]])
    for item in records:
        check(item)
    output.mkdir()
    report = dict(status="recovery_admission_running", source=record(__file__), scheduler=scheduler,
        accounting=accounting, source_report=parent_record, checked_records=records,
        seed_admitted=False, accuracy_evaluated=False, publication_ready=False)
    save_status(output / "status.json", report)
    try:
        fresh = prepare(root, output / "preflight")
        if any(fresh[key] != original_preflight[key] for key in ("context", "runtime", "saved_graph", "saved_payload")):
            raise ValueError("Fresh recovery preflight disagrees")
        if parent["runtime_after"] != fresh["runtime"]:
            raise ValueError("Recovered runtime differs from fresh preflight")
        # Load only the reviewed runner's validation function, never its CLI.
        spec = importlib.util.spec_from_file_location("reviewed_cpm_recovery", source["path"])
        runner = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(runner)
        launcher = Path(fresh["context"]["cwd"])
        overrides = dict(PYTHONPATH=str(launcher), PYTHONHASHSEED="0", OMP_NUM_THREADS="1",
                         OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONNOUSERSITE="1")
        optimizer = runner.optimizer_evidence(root, directory, fresh, overrides)
        if optimizer != parent["optimizer"]:
            raise ValueError("Independent optimizer audit disagrees")
        if any(optimizer["partition"][key] != stages[2]["output"][key] for key in ("bytes", "sha256")):
            raise ValueError("Profile partition is not the recovered optimizer output")
        expected_reports = [record(directory / name) for name in ("refinement.json", "refinement_repeat.json")]
        if parent["refinement_reports"] != expected_reports:
            raise ValueError("Recovered refinement report inventory differs")
        children = [read_frozen(Path(item["path"]), item["sha256"]) for item in expected_reports]
        reference_record = record(root / "benchmarks/work/qfo_cpm_refinement_check_22153/worker.json")
        if reference_record not in fresh["checked_records"]:
            raise ValueError("Fresh preflight lacks original refinement reference")
        reference = read_frozen(Path(reference_record["path"]), reference_record["sha256"])
        for child, filename in zip(children, ("orthogroups_profiles_refined.txt", "refinement_repeat.txt")):
            if (child["numeric_checkpoint"] != reference["numeric_checkpoint"] or child["modules"] != reference["modules"]
                    or child["output"] != record(directory / filename) or child["accuracy_evaluated"] is not False
                    or child["genes"] != 984137 or child["refinement_directed_hits"] != 0):
                raise ValueError("Recovered refinement evidence differs")
        names = (directory / "payload/gene_names.txt").read_text().splitlines()
        comparison = compare_partitions(directory / "orthogroups_profiles_refined.txt", directory / "refinement_repeat.txt",
                                        names, children[0]["groups"])
        if children[0]["groups"] != children[1]["groups"] or comparison != parent["refinement_comparison"]:
            raise ValueError("Recorded refinement comparison differs")
        repeated = output / "refinement"
        repeated.mkdir()
        (repeated / "payload").symlink_to(directory / "payload", target_is_directory=True)
        (repeated / "orthogroups_profiles.txt").symlink_to(directory / "orthogroups_profiles.txt")
        env = os.environ.copy()
        env.update(overrides)
        for key in ("PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH"):
            env.pop(key, None)
        command = [sys.executable, "-B", source["path"], "--root", str(root), "--output", str(repeated),
                   "--mode", "repeat-refinement"]
        report["refinement_command"] = command
        with (output / "refinement.log").open("x") as log:
            subprocess.run(command, cwd=launcher, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        repeated_record = record(repeated / "refinement_repeat.json")
        result = read_frozen(Path(repeated_record["path"]), repeated_record["sha256"])
        if ({k: result[k] for k in result if k != "output"} != {k: children[0][k] for k in children[0] if k != "output"}
                or result["output"] != record(repeated / "refinement_repeat.txt")):
            raise ValueError("Independent refinement execution differs")
        compare_partitions(directory / "orthogroups_profiles_refined.txt", repeated / "refinement_repeat.txt",
                           names, result["groups"])
        records = unique_records([*records, *fresh["checked_records"], *optimizer["checked_records"],
            record(output / "preflight/status.json"), repeated_record, result["output"], record(output / "refinement.log"),
            reference_record, *[child["output"] for child in children]])
        for item in records:
            check(item)
        from benchmark_tools.verify_qfo_replay_launcher import verify
        if verify(root / "benchmarks/work/publication_method_native_v2", launcher,
                  root / "benchmark_tools/results/publication_native_runtime_20260916.json") != fresh["runtime"]:
            raise ValueError("Runtime changed during recovery admission")
        report.update(status="cpm_checkpoint_recovered_seed_admitted_unscored", seed_admitted=True,
            checked_records=records, stages=stages, seed_partition=stages[-1]["output"], comparison=comparison,
            original_failure=fresh["original_scheduler"],
            limitations=["Recovered seed only: candidate preparation, phylogeny, pair conversion and scoring remain required.",
                         "Original failure retained; no controlled timing or independent generalization claim."])
    except BaseException as error:
        report.update(status="recovery_admission_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        save_status(output / "status.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--job", required=True)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.job, args.commit, args.output.absolute())
