"""Admit all four terminal QfO stage assessments with native metric and provenance checks."""

import argparse
import csv
import io
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.assemble_simulation_results import terminal_tasks
from benchmark_tools.prepare_qfo_recovered_pairs import record
from benchmark_tools.run_qfo_recovered_assessment import select_stage, command_for, checked_record, environment_records
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.validate_qfo_native_assessment import validate_directory

ENV_SHA = "e86545fd04cb644ed642cf2a31b4993fb2225ca4c05743cd03e661db129e82bc"
PAIRS_SHA = "ce6f19cd005b886a91dc3aee63cb7b41e7f448d8cad5e65ff7cb10d92a636100"
EXECUTOR = "dca732540f1ac3a56629f9cc6cdf0bbf58fc6d25"


def validate_trace(text):
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    expected = {"validate_input_file", "convertPredictions", "consolidate", "vgnc_benchmark (1)",
                "ec_benchmark (1)", "go_benchmark (1)", "fas_benchmark (1)",
                "reference_genetrees_benchmark (SwissTrees)", "reference_genetrees_benchmark (TreeFam-A)"}
    expected |= {f"scheduleMetrics ({i})" for i in range(1, 7)}
    if len(rows) != 15 or {r["name"] for r in rows} != expected or len({r["task_id"] for r in rows}) != 15:
        raise ValueError("Incomplete or duplicated Nextflow task inventory")
    for row in rows:
        expected_exit = "-" if row["name"].startswith("scheduleMetrics (") else "0"
        if row["status"] != "COMPLETED" or row["exit"] != expected_exit:
            raise ValueError("Native scoring task did not complete freshly and successfully")
    return rows


def validate_stage(root, index, scheduler, manifest, pairs, executor, env_path, pairs_path):
    stage = select_stage(pairs, index)
    directory = root / "benchmarks/results/qfo_recovered_assessment_v1" / f"stage_{index}"
    path = directory / "results.json"
    common = {"index": index, "stage": stage["stage"], "participant": stage["participant"], "scheduler": scheduler}
    if not path.exists():
        if scheduler["State"] == "COMPLETED":
            raise ValueError("Successful scoring task lacks report")
        return {**common, "status": "failed", "reason": "Unsuccessful terminal task without report"}
    report = json.loads(path.read_text())
    common["execution_report"] = record(path)
    work = root / "qfo_benchmark/w" / f"qrv2_{index}"
    results = root / "qfo_benchmark/scoring" / f"checked_v2_{index}"
    if (report["job_id"] != scheduler["JobIDRaw"] or report["array_job_id"] != "21548"
            or report["array_task_id"] != str(index) or report["stage"] != stage
            or report["source"] != record(executor / "benchmark_tools/run_qfo_recovered_assessment.py")
            or report["command"] != command_for(root, stage, manifest, work, results)
            or report["cwd"] != str(directory) or report["work"] != str(work) or report["results"] != str(results)
            or report["environment_manifest"] != record(env_path) or report["pairs_manifest"] != record(pairs_path)
            or report["accuracy_admitted"] is not False):
        raise ValueError("Scoring execution identity or command changed")
    preflight = json.loads((directory / "preflight.json").read_text())
    if preflight["status"] != "running" or any(report[k] != v for k, v in preflight.items() if k != "status"):
        raise ValueError("Scoring preflight changed")
    expected = [*environment_records(manifest), *pairs["sources"], *pairs["input_fastas"], pairs["mapping"], pairs["admission"],
                stage["partition"], stage["pairs"], stage["filtered_pairs"], stage["conversion_log"]]
    if report["verified_records"] != expected:
        raise ValueError("Scoring provenance inventory differs")
    if scheduler["State"] != "COMPLETED":
        return {**common, "status": "failed", "reason": report.get("error", report["status"])}
    if report["status"] != "process_succeeded_pending_independent_admission" or report["exit_code"] != 0:
        raise ValueError("Scheduler success contradicts scoring report")
    for item in expected:
        checked_record(item)
    actual_outputs = [record(p) for p in sorted(results.rglob("*")) if p.is_file()]
    if actual_outputs != report["outputs"]:
        raise ValueError("Scoring output inventory changed")
    checked_record(report["log"])
    traces = list((results / "stats").glob("trace_*.txt"))
    if len(traces) != 1:
        raise ValueError("Expected exactly one fresh scoring trace")
    trace_rows = validate_trace(traces[0].read_text())
    metrics, paths = validate_directory(results, stage["participant"], Path(manifest["pipeline"]) / "reference_data")
    inventoried = {r["path"] for r in actual_outputs}
    if any(str(p.resolve()) not in inventoried for p in paths if p.is_relative_to(results)):
        raise ValueError("Native metric output not inventoried")
    return {**common, "status": "admitted", "assessment": metrics,
            "metric_files": [record(p) for p in paths], "conversion": stage,
            "native_trace": record(traces[0]), "native_tasks": trace_rows}


def admit(root, output):
    if output.exists():
        raise FileExistsError(output)
    accounting = subprocess.check_output(["sacct", "-j", "21548", "--parsable2",
                                         "--format=JobID,JobIDRaw,State,ExitCode,Elapsed"], text=True)
    tasks = terminal_tasks(accounting, 21548, 4)
    env_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    pairs_path = root / "benchmarks/results/qfo_recovered_stage_pairs_v1/results.json"
    manifest, pairs = read_frozen(env_path, ENV_SHA), read_frozen(pairs_path, PAIRS_SHA)
    executor = root / "benchmarks/work/publication_qfo_recovered_assessment_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR:
        raise ValueError("Frozen scorer executor revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    for item in environment_records(manifest):
        checked_record(item)
    output.mkdir(parents=True)
    report = {"status": "validating", "publication_ready": False, "source": record(__file__), "records": [],
              "environment_manifest": record(env_path), "pairs_manifest": record(pairs_path), "accounting": accounting,
              "helper_sources": [record(Path(__file__).with_name(n)) for n in
                ("validate_qfo_native_assessment.py", "run_qfo_recovered_assessment.py", "prepare_qfo_recovered_pairs.py",
                 "assemble_simulation_results.py", "qfo_summarize_scores.py")]}
    try:
        for index, scheduler in enumerate(tasks):
            report["records"].append(validate_stage(root, index, scheduler, manifest, pairs, executor, env_path, pairs_path))
        for item in [report["source"], *report["helper_sources"], *environment_records(manifest)]:
            checked_record(item)
        read_frozen(env_path, ENV_SHA)
        read_frozen(pairs_path, PAIRS_SHA)
        report.update(status="four_stage_assessments_checked", limitations=[
            "Failed tasks are retained, not assigned zero or historical scores.",
            "Native standard errors are not paired stage-difference confidence intervals.",
            "Development-exposed stage analysis, not full HMM/phylogeny factorial or independent validation."])
    except Exception as error:
        report.update(status="validation_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.output.resolve())
