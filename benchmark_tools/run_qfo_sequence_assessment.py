"""Run frozen QfO endpoints on admitted sequence-control clique pairs."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_sequence_numeric import completed
from benchmark_tools.compare_qfo_search_coverage import frozen
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_qfo_corrected_comparator_pairs import ENV_SHA
from benchmark_tools.run_qfo_recovered_assessment import command_for, environment_records
from benchmark_tools.run_simulation_methods import read_frozen

CONVERTER = "4f0c30e5cdf287a35c9600886aec0a41bcc0b720"
WORK_NAMES = {"all_hits": "qcs_a", "top100": "qcs_t"}


def validate_stage(stage, label, scheduler):
    if (label not in WORK_NAMES or stage["variant"] != label
            or stage["status"] != "corrected_sequence_group_pairs_prepared_unscored"
            or stage["participant"] != "ohmm_qfo_corrected_sequence_" + label
            or stage["semantics"] != "cross-species group-derived clique pairs"
            or stage["accuracy_evaluated"] is not False or stage["publication_ready"] is not False):
        raise ValueError("Wrong sequence pair identity or semantics")
    if (scheduler["State"], scheduler["ExitCode"], scheduler["NodeList"], scheduler["AllocCPUS"]) != (
            "COMPLETED", "0:0", "bizon", "2") or stage["job_id"] != scheduler["JobIDRaw"]:
        raise ValueError("Wrong conversion completion or allocation")
    values = [stage[key] for key in ("total_pairs", "retained_pairs", "expected_pairs", "removed_mapping_pairs")]
    if any(type(value) is not int or value < 0 for value in values) or not values[0] == values[1] == values[2] or values[3]:
        raise ValueError("Pair accounting differs")
    if any(stage["pairs"][key] != stage["filtered_pairs"][key] for key in ("bytes", "sha256")):
        raise ValueError("Reference filtering changed predictions")
    if stage["graph_admission"] not in stage["checked_records"] or stage["prediction"] not in stage["checked_records"]:
        raise ValueError("Missing graph-admission provenance")


def prepare(root, label, pairs_sha, conversion_job, require_fresh=True):
    if label not in WORK_NAMES:
        raise ValueError("Unknown sequence variant")
    scheduler = completed(conversion_job, 2, "64G")
    path = root / "benchmarks/results/qfo_sequence_pairs_v1" / label / "results.json"
    stage = read_frozen(path, pairs_sha)
    validate_stage(stage, label, scheduler)
    converter = frozen(root, "publication_qfo_sequence_pairs_v1", CONVERTER)
    source = record(converter / "benchmark_tools/prepare_qfo_sequence_pairs.py")
    if stage["source"] != source or source not in stage["checked_records"]:
        raise ValueError("Wrong frozen converter source")
    environment_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    manifest = read_frozen(environment_path, ENV_SHA)
    if [item for item in manifest["reference_files"] if Path(item["path"]).name == "mapping.json.gz"] != [stage["mapping"]]:
        raise ValueError("Conversion and assessment mappings differ")
    directory = root / "benchmarks/results/qfo_sequence_assessment_v1" / label
    work = root / "qfo_benchmark/w" / WORK_NAMES[label]
    results = root / "qfo_benchmark/scoring" / ("corrected_sequence_" + label)
    records = [record(path), record(environment_path), *environment_records(manifest), *stage["checked_records"],
               stage["source"], stage["pairs"], stage["filtered_pairs"], stage["conversion_log"],
               record(Path(__file__).with_name("run_qfo_recovered_assessment.py"))]
    for item in records:
        check(item)
    if require_fresh:
        for candidate in (directory, work, results):
            if candidate.exists():
                raise FileExistsError(candidate)
    return dict(status="prepared_unrun", variant=label, source=record(__file__), stage=stage,
        pairs_manifest=record(path), environment_manifest=record(environment_path), conversion_scheduler=scheduler,
        converter_commit=CONVERTER, command=command_for(root, stage, manifest, work, results),
        cwd=str(directory), work=str(work), results=str(results), verified_records=records,
        environment_overrides=manifest["environment_overrides"], accuracy_admitted=False)


def run(root, label, pairs_sha, conversion_job, check_only=False):
    if not check_only and (os.environ.get("SLURM_CPUS_PER_TASK") != "8" or not os.environ.get("SLURM_JOB_ID")):
        raise ValueError("Require scheduled eight-CPU assessment")
    report = prepare(root, label, pairs_sha, conversion_job)
    if check_only:
        return report
    directory = Path(report["cwd"])
    directory.mkdir(parents=True, exist_ok=False)
    report.update(status="running", job_id=os.environ["SLURM_JOB_ID"])
    (directory / "preflight.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    environment = {**os.environ, **report["environment_overrides"],
                   "NXF_SINGULARITY_CACHEDIR": str(root / "qfo_benchmark/scoring/container_cache")}
    try:
        with (directory / "scoring.log").open("x") as log:
            done = subprocess.run(report["command"], cwd=directory, env=environment, stdout=log, stderr=subprocess.STDOUT)
        report["exit_code"] = done.returncode
        for item in [report["source"], *report["verified_records"]]:
            check(item)
        report["outputs"] = [record(path) for path in sorted(Path(report["results"]).rglob("*")) if path.is_file()]
        if done.returncode:
            raise RuntimeError("Native scoring failed; no implicit retry or fabricated score")
        report["status"] = "process_succeeded_pending_independent_admission"
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        if (directory / "scoring.log").exists():
            report["log"] = record(directory / "scoring.log")
        (directory / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--variant", choices=WORK_NAMES, required=True)
    parser.add_argument("--pairs-sha256", required=True)
    parser.add_argument("--conversion-job", required=True)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    run(args.root.resolve(), args.variant, args.pairs_sha256, args.conversion_job, args.check_only)
