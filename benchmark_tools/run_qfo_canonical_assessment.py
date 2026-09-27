"""Run the prespecified six QfO endpoints for validated canonical native pairs."""

import argparse
import json
import os
from pathlib import Path
import subprocess

from benchmark_tools.run_qfo_order_replay import record, save
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_qfo_recovered_assessment import command_for, environment_records
from benchmark_tools.prepare_qfo_corrected_comparator_pairs import ENV_SHA
from benchmark_tools.prepare_qfo_canonical_pairs import check_records, require_counts
from benchmark_tools.verify_ygob_validation import require_completed_job

VERIFIED_SHA = "6d6e6699ddbb822a81458fa0f30e79d76aa28c6cb743865c078c223f347ff3e2"


def verify_stage(stage, verified, execution):
    if (stage["status"] != "canonical_native_pairs_prepared_unscored"
            or stage["job_id"] != "22335"
            or stage["participant"] != "ohmm_qfo_canonical_20260927"
            or stage["semantics"] != "native phylogenetically inferred pairs"
            or stage["accuracy_evaluated"] is not False
            or verified["status"] != "canonical_conversion_independently_verified"
            or verified["removed_mapping_pairs"] != 0 or stage["removed_mapping_pairs"] != 0):
        raise ValueError("Wrong canonical conversion scope")
    require_counts(verified["rows_compared"], stage["total_pairs"], stage["retained_pairs"])
    for key in ("pairs", "filtered_pairs", "mapping", "native_input"):
        if stage[key] != verified[key]:
            raise ValueError("Conversion verification identity differs")
    if execution != dict(status="conversion_complete", job_id="22335",
                         plan=verified["plan"], result=verified["result"]):
        raise ValueError("Wrong conversion completion binding")


def prepare(repo, directory):
    if directory.exists():
        raise FileExistsError(directory)
    verified_path = repo / "benchmark_tools/results/qfo_canonical_conversion_verified_22335.json"
    verified = read_frozen(verified_path, VERIFIED_SHA)
    stage = read_frozen(Path(verified["result"]["path"]), verified["result"]["sha256"])
    execution = json.loads(Path(verified["execution"]["path"]).read_text())
    verify_stage(stage, verified, execution)
    accounting = subprocess.check_output(["sacct", "-j", "22335", "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, 22335)
    if scheduler["AllocCPUS"] != "2":
        raise ValueError("Wrong conversion allocation")
    env_path = repo / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment = read_frozen(env_path, ENV_SHA)
    if [x for x in environment["reference_files"] if Path(x["path"]).name == "mapping.json.gz"] != [stage["mapping"]]:
        raise ValueError("Scoring and conversion mappings differ")
    paths = sorted(set((repo / "benchmark_tools").glob("*.py")) | set((repo / "orthohmm").rglob("*.py"))
                   | set((repo / "qfo_benchmark").glob("*.py")))
    records = [record(p) for p in paths] + [record(verified_path), record(env_path),
        *environment_records(environment), *stage["checked_records"],
        *[verified[k] for k in ("plan", "execution", "result", "submission", "pairs", "filtered_pairs")]]
    check_records(records)
    output = repo / "benchmarks/results/qfo_canonical_assessment_20260927"
    work, results = repo / "qfo_benchmark/w/qcan27", repo / "qfo_benchmark/scoring/canonical_20260927"
    if any(p.exists() for p in (output, work, results)):
        raise FileExistsError("Existing canonical scoring namespace")
    plan = dict(repo=str(repo), directory=str(directory), output=str(output), work=str(work), results=str(results),
        stage=stage, conversion_verification=record(verified_path), conversion_scheduler=scheduler,
        environment_manifest=record(env_path), environment_overrides=environment["environment_overrides"],
        command=command_for(repo, stage, environment, work, results), checked_records=records,
        resources=dict(cpus=8, memory_gib=128, hours=12), attempts=1, accuracy_admitted=False)
    directory.mkdir(parents=True)
    save(directory / "plan.json", plan)
    return record(directory / "plan.json")


def run(path, sha):
    if not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "8":
        raise ValueError("Require scheduled eight-CPU scoring")
    plan = read_frozen(path, sha)
    check_records(plan["checked_records"])
    output = Path(plan["output"])
    if any(Path(plan[k]).exists() for k in ("output", "work", "results")):
        raise FileExistsError("Existing scoring namespace")
    output.mkdir(parents=True)
    report = dict(plan=record(path), status="running", job_id=os.environ["SLURM_JOB_ID"],
                  source=record(__file__), accuracy_admitted=False, command=plan["command"])
    save(output / "preflight.json", report)
    env = {**os.environ, **plan["environment_overrides"],
           "NXF_SINGULARITY_CACHEDIR": str(Path(plan["repo"]) / "qfo_benchmark/scoring/container_cache")}
    try:
        with (output / "scoring.log").open("x") as log:
            done = subprocess.run(plan["command"], cwd=output, env=env, stdout=log, stderr=subprocess.STDOUT)
        report["exit_code"] = done.returncode
        check_records(plan["checked_records"])
        if record(path) != report["plan"]:
            raise ValueError("Scoring plan changed")
        if done.returncode:
            raise RuntimeError(f"Native scoring failed: {done.returncode}")
        report["outputs"] = [record(p) for p in sorted(Path(plan["results"]).rglob("*")) if p.is_file()]
        report["status"] = "process_succeeded_pending_independent_admission"
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        if (output / "scoring.log").exists():
            report["log"] = record(output / "scoring.log")
        save(output / "results.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    args = parser.parse_args()
    run(args.plan.resolve(), args.sha256)
