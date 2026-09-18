"""Run frozen QfO endpoints for one independently converted corrected cell."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.bootstrap_qfo_factorial import CELLS
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_qfo_corrected_comparator_pairs import ENV_SHA
from benchmark_tools.run_qfo_recovered_assessment import command_for, environment_records
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_ygob_validation import require_completed_job

CONVERTERS = (
    ("group", "cb6489d6b3d06fb5c140f3c2ac297016abca0a1b", "cross-species group-derived clique pairs"),
    ("native", "1e67abba484ae22ac28639e90cd84df83e969b9f", "native phylogenetically inferred pairs"),
)


def validate_stage(stage, index, scheduler):
    if type(index) is not int or not 0 <= index < 8:
        raise ValueError("Unknown corrected factorial cell")
    kind, _, semantics = CONVERTERS[index % 2]
    if (stage["status"] != f"corrected_factorial_{kind}_pairs_prepared_unscored"
            or stage["cell"] != CELLS[index] or stage["index"] != index
            or stage["participant"] != "ohmm_qfo_corrected_factorial_" + CELLS[index]
            or stage["semantics"] != semantics or stage["accuracy_evaluated"] is not False
            or stage["publication_ready"] is not False):
        raise ValueError("Wrong corrected conversion identity or semantics")
    if (scheduler["State"], scheduler["ExitCode"], scheduler["NodeList"], scheduler["AllocCPUS"]) != (
            "COMPLETED", "0:0", "bizon", "2") or stage["job_id"] != scheduler["JobIDRaw"]:
        raise ValueError("Wrong conversion completion/resources")
    counts = [stage[k] for k in ("total_pairs", "retained_pairs", "removed_mapping_pairs")]
    if any(type(v) is not int for v in counts) or not 0 < counts[0] == counts[1] or counts[2] != 0:
        raise ValueError("Invalid corrected pair counts")
    if any(stage["pairs"][k] != stage["filtered_pairs"][k] for k in ("bytes", "sha256")):
        raise ValueError("Corrected reference filtering changed predictions")
    if not index % 2 and stage["expected_pairs"] != stage["total_pairs"]:
        raise ValueError("Group-derived pair count differs")
    if index % 2 and stage["native_admission_recheck"] not in stage["checked_records"]:
        raise ValueError("Missing fresh independent native admission")


def prepare(root, index, pairs_sha, conversion_job):
    if type(index) is not int or not 0 <= index < 8:
        raise ValueError("Unknown corrected factorial cell")
    accounting = subprocess.check_output(["sacct", "-j", str(conversion_job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, conversion_job)
    pairs_path = root / "benchmarks/results/qfo_corrected_factorial_pairs_v1" / CELLS[index] / "results.json"
    stage = read_frozen(pairs_path, pairs_sha)
    validate_stage(stage, index, scheduler)
    kind, commit, _ = CONVERTERS[index % 2]
    executor = root / f"benchmarks/work/publication_qfo_corrected_{kind}_pairs_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != commit:
        raise ValueError("Conversion executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--",
                    "benchmark_tools", "orthohmm", "qfo_benchmark/og_to_pairwise.py"], check=True)
    source = record(executor / f"benchmark_tools/prepare_qfo_corrected_{kind}_pairs.py")
    if source not in stage["checked_records"]:
        raise ValueError("Wrong frozen conversion source")
    env_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    manifest = read_frozen(env_path, ENV_SHA)
    mappings = [r for r in manifest["reference_files"] if Path(r["path"]).name == "mapping.json.gz"]
    if mappings != [stage["mapping"]]:
        raise ValueError("Conversion and scoring reference mappings differ")
    records = [record(pairs_path), record(env_path), *environment_records(manifest), *stage["checked_records"],
               stage["pairs"], stage["filtered_pairs"], record(Path(__file__).with_name("run_qfo_recovered_assessment.py"))]
    for item in records:
        check(item)
    output = root / "benchmarks/results/qfo_corrected_factorial_assessment_v1" / CELLS[index]
    work = root / "qfo_benchmark/w" / f"qcf{index}"
    results = root / "qfo_benchmark/scoring" / f"corrected_factorial_{index}"
    for path in (output, work, results):
        if path.exists():
            raise FileExistsError(path)
    return {"status": "prepared_unrun", "index": index, "cell": CELLS[index], "stage": stage,
            "source": record(__file__), "pairs_manifest": record(pairs_path), "environment_manifest": record(env_path),
            "conversion_scheduler": scheduler, "conversion_accounting": accounting,
            "converter_commit": commit, "command": command_for(root, stage, manifest, work, results),
            "cwd": str(output), "work": str(work), "results": str(results), "verified_records": records,
            "environment_overrides": manifest["environment_overrides"], "accuracy_admitted": False}


def run(root, index, pairs_sha, conversion_job, check_only=False):
    if not check_only and (os.environ.get("SLURM_CPUS_PER_TASK") != "8" or not os.environ.get("SLURM_JOB_ID")):
        raise ValueError("Require scheduled eight-CPU assessment")
    report = prepare(root, index, pairs_sha, conversion_job)
    if check_only:
        return report
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
    parser.add_argument("--index", type=int, choices=range(8), required=True)
    parser.add_argument("--pairs-sha256", required=True)
    parser.add_argument("--conversion-job", required=True)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    result = run(args.root.resolve(), args.index, args.pairs_sha256, args.conversion_job, args.check_only)
    print(json.dumps({"status": result["status"], "command": result["command"]}))
