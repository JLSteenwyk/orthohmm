"""Run the six frozen QfO endpoints on corrected native comparator pairs."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_qfo_recovered_assessment import command_for, environment_records
from benchmark_tools.prepare_qfo_corrected_comparator_pairs import ENV_SHA
from benchmark_tools.prepare_qfo_corrected_orthofinder_pairs import SEMANTICS as OF_SEMANTICS
from benchmark_tools.verify_ygob_validation import require_completed_job

WORK_NAMES = {"proteinortho": "qc_p", "sonic": "qc_s",
              "orthofinder_full": "qc_of", "orthofinder_sequence_only": "qc_om"}
METHODS = tuple(WORK_NAMES)
OF_CONVERTER = "aa8da7800c4801684726151ad249da7c82a3b88d"


def converter_source(root, method):
    if method not in METHODS:
        raise ValueError("Unknown corrected comparator")
    if method not in OF_SEMANTICS:
        return None
    executor = root / "benchmarks/work/publication_qfo_corrected_orthofinder_pairs_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != OF_CONVERTER:
        raise ValueError("OrthoFinder converter executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    return record(executor / "benchmark_tools/prepare_qfo_corrected_orthofinder_pairs.py")


def validate_stage(stage, method, scheduler):
    if method not in METHODS:
        raise ValueError("Unknown corrected comparator")
    status = ("corrected_orthofinder_pairs_prepared_unscored" if method in OF_SEMANTICS
              else "corrected_comparator_pairs_prepared_unscored")
    if (stage["status"] != status or stage["accuracy_evaluated"] is not False
            or stage["method"] != method or stage["participant"] != "qfo_corrected_" + method):
        raise ValueError("Wrong corrected conversion identity/status")
    if method in OF_SEMANTICS and (stage["semantics"] != OF_SEMANTICS[method]
                                   or stage["publication_ready"] is not False):
        raise ValueError("Wrong OrthoFinder prediction semantics")
    if (scheduler["State"] != "COMPLETED" or scheduler["ExitCode"] != "0:0"
            or scheduler["NodeList"] != "bizon" or scheduler["AllocCPUS"] != "2"
            or stage["job_id"] != scheduler["JobIDRaw"]):
        raise ValueError("Wrong conversion scheduler identity/allocation")
    if (any(type(stage[k]) is not int for k in ("total_pairs", "retained_pairs", "removed_mapping_pairs"))
            or stage["total_pairs"] <= 0
            or stage["retained_pairs"] != stage["total_pairs"] or stage["removed_mapping_pairs"] != 0):
        raise ValueError("Invalid corrected pair counts")
    if any(stage["pairs"][k] != stage["filtered_pairs"][k] for k in ("bytes", "sha256")):
        raise ValueError("Unexpected corrected reference-filter change")


def prepare(root, method, pairs_sha, conversion_job):
    if method not in METHODS:
        raise ValueError("Unknown corrected comparator")
    accounting = subprocess.check_output(["sacct", "-j", str(conversion_job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, conversion_job)
    pairs_path = root / "benchmarks/results/qfo_corrected_comparator_pairs_v1" / method / "results.json"
    stage = read_frozen(pairs_path, pairs_sha)
    validate_stage(stage, method, scheduler)
    converter = converter_source(root, method)
    if converter is not None and (stage["source"] != converter or converter not in stage["checked_records"]):
        raise ValueError("Wrong frozen OrthoFinder conversion source")
    env_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    manifest = read_frozen(env_path, ENV_SHA)
    if method in OF_SEMANTICS and [r for r in manifest["reference_files"]
            if Path(r["path"]).name == "mapping.json.gz"] != [stage["mapping"]]:
        raise ValueError("Conversion and assessment reference mappings differ")
    records = [record(pairs_path), record(env_path), *environment_records(manifest),
               stage["source"], *stage["checked_records"], stage["pairs"], stage["filtered_pairs"],
               record(Path(__file__).with_name("run_qfo_recovered_assessment.py")),
               record(Path(__file__).with_name("prepare_qfo_corrected_comparator_pairs.py"))]
    if converter is not None:
        records.append(converter)
    for item in records:
        check(item)
    output = root / "benchmarks/results/qfo_corrected_assessment_v1" / method
    work = root / "qfo_benchmark/w" / WORK_NAMES[method]
    results = root / "qfo_benchmark/scoring" / ("corrected_" + method)
    for path in (output, work, results):
        if path.exists():
            raise FileExistsError(path)
    return {"status": "prepared_unrun", "method": method, "stage": stage,
            "source": record(__file__), "pairs_manifest": record(pairs_path), "environment_manifest": record(env_path),
            "conversion_scheduler": scheduler, "conversion_accounting": accounting,
            "command": command_for(root, stage, manifest, work, results), "cwd": str(output),
            "work": str(work), "results": str(results), "verified_records": records,
            "environment_overrides": manifest["environment_overrides"], "accuracy_admitted": False}


def run(root, method, pairs_sha, conversion_job, check_only=False):
    report = prepare(root, method, pairs_sha, conversion_job)
    if check_only:
        return report
    if os.environ.get("SLURM_CPUS_PER_TASK") != "8" or not os.environ.get("SLURM_JOB_ID"):
        raise ValueError("Require scheduled eight-CPU assessment")
    output = Path(report["cwd"])
    output.mkdir(parents=True, exist_ok=False)
    report.update(status="running", job_id=os.environ["SLURM_JOB_ID"])
    with (output / "preflight.json").open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    env = {**os.environ, **report["environment_overrides"],
           "NXF_SINGULARITY_CACHEDIR": str(root / "qfo_benchmark/scoring/container_cache")}
    try:
        with (output / "scoring.log").open("x") as log:
            done = subprocess.run(report["command"], cwd=output, env=env, stdout=log, stderr=subprocess.STDOUT)
        report["exit_code"] = done.returncode
        for item in [report["source"], *report["verified_records"]]:
            check(item)
        report["outputs"] = [record(p) for p in sorted(Path(report["results"]).rglob("*")) if p.is_file()]
        if done.returncode:
            raise RuntimeError(f"Native scoring failed: {done.returncode}")
        report["status"] = "process_succeeded_pending_independent_admission"
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        log_path = output / "scoring.log"
        if log_path.exists():
            report["log"] = record(log_path)
        with (output / "results.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--method", choices=METHODS, required=True)
    parser.add_argument("--pairs-sha256", required=True)
    parser.add_argument("--conversion-job", type=int, required=True)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    report = run(args.root.resolve(), args.method, args.pairs_sha256, args.conversion_job, args.check_only)
    print(json.dumps({"status": report["status"], "command": report["command"]}))
