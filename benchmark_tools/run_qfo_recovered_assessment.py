"""Assess one frozen recovered QfO stage in a fresh, non-destructive namespace."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_qfo_recovered_pairs import record, STAGES
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_ygob_validation import require_completed_job


def select_stage(pairs, index):
    if isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < 4:
        raise ValueError("Require stage index0-3")
    rows = pairs["stages"]
    if (pairs["status"] != "four_stage_pairs_prepared_unscored" or pairs["accuracy_evaluated"] is not False
            or len(rows) != 4 or [r["stage"] for r in rows] != list(STAGES)):
        raise ValueError("Require complete unscored four-stage pair inventory")
    for i, row in enumerate(rows):
        if (row["index"] != i or row["participant"] != f"ohmm_checked_v2_{i}"
                or not 0 < row["retained_pairs"] <= row["total_pairs"]
                or row["removed_mapping_pairs"] != row["total_pairs"] - row["retained_pairs"]):
            raise ValueError("Pair stage identity or counts differ")
    return rows[index]


def checked_record(item):
    if record(item["path"]) != item:
        raise ValueError("Scoring input/runtime checksum changed: " + item["path"])


def environment_records(manifest):
    if manifest["status"] != "local_qfo_assessment_environment_frozen" or manifest["accuracy_evaluated"] is not False:
        raise ValueError("Require frozen local scorer manifest")
    return [manifest["source"], manifest["execution_config"], manifest["singularity_config"],
            *[r for key in ("pipeline_files", "reference_files", "java_files", "images", "executables", "singularity_support")
              for r in manifest[key]]]


def command_for(root, stage, manifest, work, results):
    if len(str(work)) + len("/xx/xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx/predictions.db") > 159:
        raise ValueError("Darwin work path exceeds validated length limit")
    pipeline = Path(manifest["pipeline"])
    return ["/home/bizon/bin/nextflow", "-c", manifest["execution_config"]["path"], "run", str(pipeline / "main.nf"),
            "-profile", "singularity", "-ansi-log", "false", "--input", stage["filtered_pairs"]["path"],
            "--participant_id", stage["participant"], "--event_year", "2020",
            "--challenges_ids", "GO EC VGNC SwissTrees TreeFam-A FAS", "--goldstandard_dir", str(pipeline / "reference_data/2020"),
            "--assess_dir", str(pipeline / "reference_data/data"), "--results_dir", str(results), "-work-dir", str(work)]


def run(root, index, environment_sha, pairs_sha):
    manifest_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    pairs_path = root / "benchmarks/results/qfo_recovered_stage_pairs_v1/results.json"
    manifest, pairs = read_frozen(manifest_path, environment_sha), read_frozen(pairs_path, pairs_sha)
    stage = select_stage(pairs, index)
    accounting = subprocess.check_output(["sacct", "-j", "21522", "--parsable2", "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 21522)
    if pairs["job_id"] != "21522":
        raise ValueError("Pair preparation belongs to another job")
    output = root / "benchmarks/results/qfo_recovered_assessment_v1" / f"stage_{index}"
    work = root / "qfo_benchmark/w" / f"qrv2_{index}"
    results = root / "qfo_benchmark/scoring" / f"checked_v2_{index}"
    for path in (output, work, results):
        if path.exists():
            raise FileExistsError(path)
    records = [*environment_records(manifest), *pairs["sources"], *pairs["input_fastas"], pairs["mapping"], pairs["admission"],
               stage["partition"], stage["pairs"], stage["filtered_pairs"], stage["conversion_log"]]
    for item in records:
        checked_record(item)
    command = command_for(root, stage, manifest, work, results)
    env = {**os.environ, **manifest["environment_overrides"],
           "NXF_SINGULARITY_CACHEDIR": str(root / "qfo_benchmark/scoring/container_cache")}
    output.mkdir(parents=True)
    report = {"status": "running", "accuracy_admitted": False, "source": record(__file__), "stage": stage,
              "environment_manifest": record(manifest_path), "pairs_manifest": record(pairs_path),
              "pair_preparation_scheduler": scheduler, "command": command, "cwd": str(output),
              "job_id": os.environ.get("SLURM_JOB_ID"), "array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"),
              "array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"), "work": str(work), "results": str(results),
              "verified_records": records}
    (output / "preflight.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    try:
        with (output / "scoring.log").open("x") as log:
            process = subprocess.run(command, cwd=output, env=env, stdout=log, stderr=subprocess.STDOUT)
        report["exit_code"] = process.returncode
        for item in records:
            checked_record(item)
        read_frozen(manifest_path, environment_sha)
        read_frozen(pairs_path, pairs_sha)
        checked_record(report["source"])
        report["outputs"] = [record(p) for p in sorted(results.rglob("*")) if p.is_file()]
        report["status"] = "process_succeeded_pending_independent_admission" if process.returncode == 0 else "scoring_failed"
    except Exception as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["log"] = record(output / "scoring.log")
        (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--environment-sha256", required=True)
    parser.add_argument("--pairs-sha256", required=True)
    args = parser.parse_args()
    result = run(args.root.resolve(), args.index, args.environment_sha256, args.pairs_sha256)
    raise SystemExit(0 if result["status"] == "process_succeeded_pending_independent_admission" else 1)
