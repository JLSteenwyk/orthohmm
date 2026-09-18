"""Independently admit a fresh terminal six-challenge QfO factorial assessment."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_qfo_factorial_assessment import cell_label, verify_conversion_identity
from benchmark_tools.run_qfo_recovered_assessment import command_for, environment_records
from benchmark_tools.admit_qfo_recovered_assessment import ENV_SHA, validate_trace
from benchmark_tools.validate_qfo_native_assessment import validate_directory
from benchmark_tools.verify_ygob_validation import require_completed_job

EXECUTOR_COMMIT = "25f328d994765369cfae0382a21c3e7fdb3b7dab"


def check_completion(report, preflight, scheduler, index):
    if index in (0, 4):
        raise ValueError("Reused assessments require their separate revalidation evidence")
    if (report["cell"] != cell_label(index) or report["index"] != index or
            report["job_id"] != scheduler["JobIDRaw"] or scheduler["State"] != "COMPLETED" or
            scheduler["ExitCode"] != "0:0" or scheduler["NodeList"] != "bizon" or scheduler["AllocCPUS"] != "8"):
        raise ValueError("Wrong assessment scheduler identity or allocation")
    if (report["status"] != "process_succeeded_pending_independent_admission" or report.get("exit_code") != 0
            or report.get("accuracy_admitted") is not False):
        raise ValueError("Fresh successful assessment not established")
    if preflight["status"] != "running" or any(report.get(k) != v for k, v in preflight.items() if k != "status"):
        raise ValueError("Assessment preflight changed")


def admit(root, index, job, pairs_sha):
    if Path.cwd().resolve() != root:
        raise ValueError("Run from repository verification directory")
    label = cell_label(index)
    directory = root / "benchmarks/results/qfo_factorial_assessment_v1" / label
    report_path, preflight_path = directory / "results.json", directory / "preflight.json"
    checked = [record(report_path), record(preflight_path)]
    report = json.loads(report_path.read_text())
    preflight = json.loads(preflight_path.read_text())
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, job)
    check_completion(report, preflight, scheduler, index)
    pairs_path = root / "benchmarks/results/qfo_factorial_pairs_v1" / label / "results.json"
    stage = read_frozen(pairs_path, pairs_sha)
    array = 21675 if index % 2 else 21674
    conversion_accounting = subprocess.check_output(["sacct", "-j", str(array), "--parsable2",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed"], text=True)
    conversion_scheduler = verify_conversion_identity(stage, index, conversion_accounting)
    env_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment = read_frozen(env_path, ENV_SHA)
    executor = root / "benchmarks/work/publication_qfo_factorial_assessment_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR_COMMIT:
        raise ValueError("Assessment executor revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    work = root / "qfo_benchmark/w" / f"qfx_{index}"
    results = root / "qfo_benchmark/scoring" / f"factorial_v1_{index}"
    expected = {"stage": stage, "source": record(executor / "benchmark_tools/run_qfo_factorial_assessment.py"),
                "pairs_manifest": record(pairs_path), "environment_manifest": record(env_path),
                "conversion_scheduler": conversion_scheduler, "command": command_for(root, stage, environment, work, results),
                "cwd": str(directory), "work": str(work), "results": str(results),
                "helper_sources": [record(executor / "benchmark_tools" / name) for name in
                    ("admit_qfo_recovered_assessment.py", "run_qfo_recovered_assessment.py", "validate_qfo_native_assessment.py")]}
    for key, value in expected.items():
        if report[key] != value:
            raise ValueError(f"Changed assessment execution provenance: {key}")
    records = [*environment_records(environment), *stage["checked_records"], stage["pairs"], stage["filtered_pairs"], record(pairs_path)]
    if report["verified_records"] != records:
        raise ValueError("Assessment verification inventory differs")
    checked.extend([*records, report["source"], *report["helper_sources"], report["log"]])
    for item in checked:
        check(item)
    outputs = [record(p) for p in sorted(results.rglob("*")) if p.is_file()]
    if outputs != report["outputs"]:
        raise ValueError("Assessment output set or content changed")
    traces = list((results / "stats").glob("trace_*.txt"))
    if len(traces) != 1:
        raise ValueError("Require exactly one native scoring trace")
    trace = validate_trace(traces[0].read_text())
    assessment, paths = validate_directory(results, stage["participant"], Path(environment["pipeline"]) / "reference_data")
    inventory = {item["path"] for item in outputs}
    if any(str(p.resolve()) not in inventory for p in paths if p.is_relative_to(results)):
        raise ValueError("Native metric missing from execution inventory")
    for item in checked:
        check(item)
    return {"status": "fresh_factorial_assessment_admitted", "cell": label, "index": index,
            "accuracy_admitted": True, "publication_ready": False, "source": record(__file__),
            "scheduler": scheduler, "accounting": accounting, "conversion_accounting": conversion_accounting,
            "execution_report": checked[0], "preflight": checked[1], "pairs_manifest": record(pairs_path),
            "environment_manifest": record(env_path), "conversion": stage, "assessment": assessment,
            "metric_files": [record(p) for p in paths], "native_trace": record(traces[0]), "native_tasks": trace,
            "helper_sources": [record(Path(__file__).with_name(n)) for n in
                ("run_qfo_factorial_assessment.py", "run_qfo_recovered_assessment.py", "admit_qfo_recovered_assessment.py", "validate_qfo_native_assessment.py")],
            "limitations": ["Native score/provenance admission; not a paired uncertainty analysis or independent biological validation.",
                "Six-metric mean is project-defined secondary summary, not an official QfO F1.",
                "FAS retains native sampling and endpoint-specific stderr semantics; no new uncertainty interpretation."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--index", type=int, choices=range(8), required=True)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--pairs-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = admit(args.root.resolve(), args.index, args.job, args.pairs_sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
