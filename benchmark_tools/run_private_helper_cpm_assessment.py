"""Run unchanged QfO endpoints for losslessly converted recovered private CPM pairs."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import prepare_private_helper_cpm_pairs as converter
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_qfo_recovered_assessment import command_for, environment_records
from benchmark_tools.run_simulation_methods import read_frozen

CONVERTER_COMMIT = "b502bc6965d37916224c87595a46ef8ab5c9c4d1"
CONVERTER_SHA = "a5f9fd7de5ba163c7cdfcfc487ece50ecea7e4edb4923f7d7ffe07f330fd918f"
CONVERTER_PROTOCOL_SHA = "e8bda78ed6d4c2b62652a88791771d51c309e7a755517d9f0b1d32bdf2964456"
CONVERTER_EXECUTOR = "benchmarks/work/qfo_private_cpm_pair_conversion_executor_20261001"
PROTOCOL = "benchmark_tools/results/QFO_PRIVATE_CPM_ASSESSMENT_PROTOCOL_20261001.md"
OUTPUT = "benchmarks/results/qfo_cpm_private_assessment_v1/cpm_high"
WORK = "qfo_benchmark/w/qcpp1"
RESULTS = "qfo_benchmark/scoring/private_cpm_high_v1"


def bound_record(records, path):
    matches = [ref for ref in records if ref["path"] == str(path)]
    if not matches or any(ref != matches[0] for ref in matches):
        raise ValueError("Missing or conflicting converted evidence binding")
    return matches[0]


def validate_stage(stage, scheduler, source):
    if (stage["status"] != "private_recovered_cpm_native_pairs_prepared_unscored"
            or stage["arm"] != "cpm_high" or type(stage["index"]) is not int or stage["index"] != 1
            or stage["participant"] != "ohmm_qfo_parameter_cpm_high"
            or stage["semantics"] != "native phylogenetically inferred pairs"
            or stage["job_id"] != scheduler["JobIDRaw"] or stage["source"] != source
            or any(stage[key] is not False for key in
                   ("accuracy_evaluated", "scoring_admitted", "controlled_timing", "publication_ready"))):
        raise ValueError("Wrong recovered private conversion identity or semantics")
    counts = [stage[key] for key in ("written_pairs", "total_pairs", "retained_pairs", "removed_mapping_pairs")]
    if any(type(value) is not int for value in counts) or not 0 < counts[0] == counts[1] == counts[2] or counts[3] != 0:
        raise ValueError("Invalid recovered private conversion/mapping counts")
    if any(stage["pairs"][key] != stage["filtered_pairs"][key] for key in ("bytes", "sha256")):
        raise ValueError("Recovered private reference filtering changed predictions")
    for key in ("source", "protocol", "pairs", "filtered_pairs", "native_admission", "native_admission_recheck", "conversion_counts"):
        if stage[key] not in stage["checked_records"]:
            raise ValueError("Missing checked recovered private conversion evidence")


def verify_conversion(root, job, conversion_sha, submission_sha):
    observed = converter.accounting(job)
    scheduler = converter.completed_admission(observed, job)
    path = root / converter.OUTPUT / "results.json"
    stage = read_frozen(path, conversion_sha)
    executor = root / CONVERTER_EXECUTOR
    source = record(executor / "benchmark_tools/prepare_private_helper_cpm_pairs.py")
    if source["sha256"] != CONVERTER_SHA or record(converter.__file__)["sha256"] != CONVERTER_SHA:
        raise ValueError("Frozen or reused private converter source changed")
    validate_stage(stage, scheduler, source)
    for key, name in (("pairs", "pairs.tsv"), ("filtered_pairs", "pairs.qfo.tsv"),
                      ("conversion_counts", "conversion_counts.json")):
        if stage[key] != record(root / converter.OUTPUT / name):
            raise ValueError("Recovered private converted output path or bytes differ")
    preflight_ref = bound_record(stage["checked_records"], root / converter.OUTPUT / "preflight.json")
    preflight = read_frozen(Path(preflight_ref["path"]), preflight_ref["sha256"])
    if (preflight["status"] != "preparing"
            or any(stage.get(key) != value for key, value in preflight.items() if key not in ("status", "checked_records"))
            or stage["checked_records"][:len(preflight["checked_records"])] != preflight["checked_records"]):
        raise ValueError("Recovered private conversion preflight differs")
    protocol = record(root / converter.PROTOCOL)
    if stage["protocol"] != protocol or protocol["sha256"] != CONVERTER_PROTOCOL_SHA:
        raise ValueError("Recovered private conversion protocol changed")
    submission_path = root / f"benchmark_tools/results/qfo_private_cpm_pair_conversion_submission_{job}.json"
    submission = read_frozen(submission_path, submission_sha)
    native_submission = bound_record(stage["checked_records"], root /
        f"benchmark_tools/results/qfo_private_cpm_native_admission_submission_{stage['admission_job']}.json")
    script = "benchmark_tools/results/qfo_private_cpm_pair_conversion_20261001.sh"
    if (submission["status"] != "private_recovered_qfo_high_cpm_pair_conversion_submitted"
            or submission["job_id"] != job or submission["executor"] != str(executor)
            or submission["executor_commit"] != CONVERTER_COMMIT or submission["executor_clean"] is not True
            or any(submission[key] is not False for key in
                   ("native_inference", "accuracy_evaluated", "controlled_timing", "publication_ready"))
            or submission["native_admission"] != stage["native_admission"]
            or submission["parent_submission"] != native_submission
            or submission["submission_argv"] != ["sbatch", "--parsable", str(executor / script), str(executor),
                CONVERTER_COMMIT, CONVERTER_PROTOCOL_SHA, stage["admission_job"],
                stage["native_admission"]["sha256"], native_submission["sha256"]]):
        raise ValueError("Recovered private conversion submission differs")
    names = ("benchmark_tools/prepare_private_helper_cpm_pairs.py", converter.PROTOCOL,
             script, "tests/unit/test_prepare_private_helper_cpm_pairs.py")
    submitted_sources = [record(executor / name) for name in names]
    if submission["source_records"] != submitted_sources:
        raise ValueError("Recovered private conversion submitted source inventory differs")
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != CONVERTER_COMMIT:
        raise ValueError("Recovered private converter revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm",
                    "tests/unit/test_prepare_private_helper_cpm_pairs.py"], check=True, capture_output=True)
    native_accounting = converter.accounting(stage["admission_job"])
    native_scheduler = converter.completed_admission(native_accounting, stage["admission_job"])
    checker, native_sources = converter.admission_submission(root, stage["admission_job"], native_submission["sha256"])
    if (stage["admission_scheduler"] != native_scheduler or stage["native_admission"] != record(root / converter.ADMISSION)
            or stage["native_admission_recheck"] != record(root / converter.OUTPUT / "native_admission_recheck.json")):
        raise ValueError("Recovered private native admission binding differs")
    native = read_frozen(Path(stage["native_admission"]["path"]), stage["native_admission"]["sha256"])
    fresh = read_frozen(Path(stage["native_admission_recheck"]["path"]), stage["native_admission_recheck"]["sha256"])
    if fresh != native:
        raise ValueError("Converted fresh native admission differs")
    producer, verified = converter.native_checker.producer(root, root / converter.native_checker.EXECUTOR)
    cell, _, _, _ = producer.native_command(verified, root / producer.OUTPUT)
    converter.validate_native(native, checker, verified, cell, root)
    if (stage["native_input"] != native["native_pairs"] or stage["written_pairs"] != native["native_pair_count"]
            or stage["input_fastas"] != verified["baseline"]["manifest"]["input_fastas"]
            or stage["candidate_admission"] != verified["candidates"]["admission_record"]
            or stage["private_native_admission"] != verified["private_control"]["admission"]):
        raise ValueError("Recovered private conversion input context differs")
    command = [converter.PYTHON, "-B", checker["path"], "--root", str(root),
        "--submission-sha256", converter.NATIVE_SUBMISSION_SHA, "--protocol-sha256", converter.NATIVE_PROTOCOL_SHA,
        "--output", stage["native_admission_recheck"]["path"]]
    if stage["native_recheck_command"] != command:
        raise ValueError("Recovered private conversion recheck command differs")
    counts = read_frozen(Path(stage["conversion_counts"]["path"]), stage["conversion_counts"]["sha256"])
    if counts != {key: stage[key] for key in ("written_pairs", "total_pairs", "retained_pairs", "removed_mapping_pairs")}:
        raise ValueError("Recovered private conversion count sidecar differs")
    checked = [record(path), record(submission_path), *submitted_sources, *native_sources,
               *stage["checked_records"], *verified["checked_records"]]
    for ref in checked:
        check(ref)
    return dict(stage=stage, scheduler=scheduler, accounting=observed, pairs_manifest=record(path),
        submission=record(submission_path), checked_records=checked,
        context=dict(candidate_admission=stage["candidate_admission"], private_native_admission=stage["private_native_admission"],
                     recovered_arm=verified["candidates"]["arm"]))


def prepare(root, job, conversion_sha, submission_sha, protocol_sha):
    output, work, results = (root / relative for relative in (OUTPUT, WORK, RESULTS))
    for namespace in (output, work, results):
        if namespace.resolve() != namespace or namespace.exists() or namespace.is_symlink():
            raise FileExistsError("Recovered private assessment namespaces must be fresh and direct")
    verified = verify_conversion(root, job, conversion_sha, submission_sha)
    stage = verified["stage"]
    protocol = record(root / PROTOCOL)
    if protocol["sha256"] != protocol_sha:
        raise ValueError("Recovered private assessment protocol changed")
    environment_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment = read_frozen(environment_path, converter.ENV_SHA)
    if [ref for ref in environment["reference_files"] if Path(ref["path"]).name == "mapping.json.gz"] != [stage["mapping"]]:
        raise ValueError("Recovered private conversion and scoring mapping differ")
    helpers = [record(module.__file__) for name, module in sorted(sys.modules.items())
        if name.startswith("benchmark_tools.") and getattr(module, "__file__", None)]
    base = [protocol, verified["pairs_manifest"], verified["submission"], record(environment_path),
            *environment_records(environment), *verified["checked_records"]]
    for ref in [record(__file__), *base, *helpers]:
        check(ref)
    return dict(status="prepared_unrun", arm="cpm_high", index=1, context=verified["context"], stage=stage,
        source=record(__file__), protocol=protocol, pairs_manifest=verified["pairs_manifest"],
        conversion_submission=verified["submission"], environment_manifest=record(environment_path),
        conversion_scheduler=verified["scheduler"], conversion_accounting=verified["accounting"],
        conversion_job=job, conversion_sha256=conversion_sha, conversion_submission_sha256=submission_sha,
        converter_commit=CONVERTER_COMMIT, command=command_for(root, stage, environment, work, results),
        cwd=str(output), work=str(work), results=str(results), verified_records=[*base, *helpers], helpers=helpers,
        environment_overrides=environment["environment_overrides"], accuracy_admitted=False,
        controlled_timing=False, publication_ready=False), verified


def run(root, job, conversion_sha, submission_sha, protocol_sha):
    if (not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "8"
            or os.environ.get("SLURM_JOB_NODELIST") != "bizon" or os.environ.get("SLURM_MEM_PER_NODE") != "65536"
            or os.environ.get("SLURM_ARRAY_TASK_ID") or sys.executable != converter.PYTHON):
        raise ValueError("Require standalone eight-CPU/64-GiB bizon assessment controller")
    if Path.cwd().resolve() != root:
        raise ValueError("Run from original repository verification directory")
    report, verified = prepare(root, job, conversion_sha, submission_sha, protocol_sha)
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
            process = subprocess.run(report["command"], cwd=output, env=env, stdout=log, stderr=subprocess.STDOUT)
        report["exit_code"] = process.returncode
        if process.returncode:
            raise RuntimeError("Recovered private QfO scoring failed; preserve evidence without retry")
        for ref in [report["source"], *report["verified_records"]]:
            check(ref)
        if verify_conversion(root, job, conversion_sha, submission_sha) != verified:
            raise ValueError("Recovered private conversion evidence/runtime changed during scoring")
        for ref in [report["source"], *report["verified_records"]]:
            check(ref)
        report["outputs"] = [record(path) for path in sorted(Path(report["results"]).rglob("*")) if path.is_file()]
        report["status"] = "process_succeeded_pending_independent_admission"
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        results = Path(report["results"])
        if report["status"] == "failed" and results.exists():
            try:
                report["outputs"] = [record(path) for path in sorted(results.rglob("*")) if path.is_file()]
            except Exception as error:
                report["output_capture_error"] = {"type": type(error).__name__, "message": str(error)}
        if (output / "scoring.log").exists():
            report["log"] = record(output / "scoring.log")
        with (output / "results.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--conversion-job", required=True)
    parser.add_argument("--conversion-sha256", required=True)
    parser.add_argument("--conversion-submission-sha256", required=True)
    parser.add_argument("--protocol-sha256", required=True)
    args = parser.parse_args()
    run(args.root.resolve(), args.conversion_job, args.conversion_sha256,
        args.conversion_submission_sha256, args.protocol_sha256)
