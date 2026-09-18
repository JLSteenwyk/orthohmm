"""Independently validate terminal corrected-comparator QfO assessments."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_qfo_corrected_comparator_assessment import (
    validate_stage, ENV_SHA, METHODS, WORK_NAMES, OF_SEMANTICS, converter_source,
)
from benchmark_tools.run_qfo_recovered_assessment import command_for, environment_records
from benchmark_tools.admit_qfo_recovered_assessment import validate_trace
from benchmark_tools.validate_qfo_native_assessment import validate_directory
from benchmark_tools.verify_ygob_validation import require_completed_job

EXECUTOR = "74afad5376b7ee11fdabfba386851fd8d3c02857"
OF_EXECUTOR = "5c34f8baad47a9895659b43696e6887aa29dbaaf"


def executor_identity(root, method):
    if method not in METHODS:
        raise ValueError("Unknown corrected comparator")
    if method in OF_SEMANTICS:
        return root / "benchmarks/work/publication_qfo_corrected_of_assessment_v1", OF_EXECUTOR
    return root / "benchmarks/work/publication_qfo_corrected_assessment_v1", EXECUTOR


def validate_completion(report, preflight, scheduler):
    if (scheduler["State"] != "COMPLETED" or scheduler["ExitCode"] != "0:0"
            or scheduler["NodeList"] != "bizon" or scheduler["AllocCPUS"] != "8"
            or report["job_id"] != scheduler["JobIDRaw"]):
        raise ValueError("Assessment scheduler mismatch")
    if (report["status"] != "process_succeeded_pending_independent_admission"
            or report.get("exit_code") != 0 or report.get("accuracy_admitted") is not False):
        raise ValueError("Require successful unadmitted assessment")
    if preflight["status"] != "running" or any(report.get(k) != v for k, v in preflight.items() if k != "status"):
        raise ValueError("Changed assessment preflight")


def accounting(job):
    text = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    return text, require_completed_job(text, job)


def admit(root, method, job, conversion_job, pairs_sha, output):
    executor, executor_commit = executor_identity(root, method)
    if output.exists():
        raise FileExistsError(output)
    text, scheduler = accounting(job)
    directory = root / "benchmarks/results/qfo_corrected_assessment_v1" / method
    report_path, preflight_path = directory / "results.json", directory / "preflight.json"
    checked = [record(report_path), record(preflight_path)]
    report, preflight = (json.loads(p.read_text()) for p in (report_path, preflight_path))
    validate_completion(report, preflight, scheduler)
    pairs_path = root / "benchmarks/results/qfo_corrected_comparator_pairs_v1" / method / "results.json"
    stage = read_frozen(pairs_path, pairs_sha)
    conversion_text, conversion_scheduler = accounting(conversion_job)
    validate_stage(stage, method, conversion_scheduler)
    converter = converter_source(root, method)
    if converter is not None and (stage["source"] != converter or converter not in stage["checked_records"]):
        raise ValueError("Wrong frozen OrthoFinder conversion source")
    env_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    manifest = read_frozen(env_path, ENV_SHA)
    if method in OF_SEMANTICS and [r for r in manifest["reference_files"]
            if Path(r["path"]).name == "mapping.json.gz"] != [stage["mapping"]]:
        raise ValueError("Conversion and assessment reference mappings differ")
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != executor_commit:
        raise ValueError("Assessment executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    work = root / "qfo_benchmark/w" / WORK_NAMES[method]
    results = root / "qfo_benchmark/scoring" / ("corrected_" + method)
    records = [record(pairs_path), record(env_path), *environment_records(manifest), stage["source"],
               *stage["checked_records"], stage["pairs"], stage["filtered_pairs"],
               record(executor / "benchmark_tools/run_qfo_recovered_assessment.py"),
               record(executor / "benchmark_tools/prepare_qfo_corrected_comparator_pairs.py")]
    if converter is not None:
        records.append(converter)
    expected = {"method": method, "stage": stage, "source": record(executor / "benchmark_tools/run_qfo_corrected_comparator_assessment.py"),
                "pairs_manifest": record(pairs_path), "environment_manifest": record(env_path),
                "conversion_scheduler": conversion_scheduler, "conversion_accounting": conversion_text,
                "command": command_for(root, stage, manifest, work, results), "cwd": str(directory),
                "work": str(work), "results": str(results), "verified_records": records,
                "environment_overrides": manifest["environment_overrides"]}
    for key, value in expected.items():
        if report[key] != value:
            raise ValueError("Changed execution provenance: " + key)
    if report["log"]["path"] != str(directory / "scoring.log"):
        raise ValueError("Unexpected log path")
    checked.extend([report["source"], report["log"], *records])
    for item in checked:
        check(item)
    observed = [record(p) for p in sorted(results.rglob("*")) if p.is_file()]
    if observed != report["outputs"]:
        raise ValueError("Changed output set/content")
    traces = list((results / "stats").glob("trace_*.txt"))
    if len(traces) != 1:
        raise ValueError("Require exactly one native scoring trace")
    tasks = validate_trace(traces[0].read_text())
    assessment, paths = validate_directory(results, stage["participant"], Path(manifest["pipeline"]) / "reference_data")
    metric_records = [record(p) for p in paths]
    for item in [*checked, *observed, *metric_records]:
        check(item)
    result = {"status": "corrected_comparator_assessment_admitted", "method": method,
              "source": record(__file__), "scheduler": scheduler, "accounting": text,
              "conversion_scheduler": conversion_scheduler, "pairs_manifest": record(pairs_path),
              "execution_report": record(report_path), "preflight": record(preflight_path),
              "assessment": assessment, "metric_files": metric_records, "native_tasks": tasks,
              "native_trace": record(traces[0]), "checked_records": checked,
              "validator_sources": [record(Path(__file__).with_name(name)) for name in (
                  "run_qfo_corrected_comparator_assessment.py", "run_qfo_recovered_assessment.py",
                  "admit_qfo_recovered_assessment.py", "validate_qfo_native_assessment.py")],
              "accuracy_admitted": True, "publication_ready": False,
              "limitations": ["Native score/provenance validation, not independent biological validation or paired uncertainty.",
                              "The six-endpoint mean is a project-defined secondary summary, not official QfO F1.",
                              "Native error fields retain endpoint-specific semantics; original-release scores are not reused."]}
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--method", choices=METHODS, required=True)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--conversion-job", type=int, required=True)
    parser.add_argument("--pairs-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.method, args.job, args.conversion_job, args.pairs_sha256, args.output.resolve())
