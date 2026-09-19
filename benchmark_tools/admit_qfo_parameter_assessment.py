"""Independently validate terminal parameter-variant QfO assessments."""

import argparse
import csv
import io
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_corrected_comparator_assessment import validate_completion
from benchmark_tools.admit_qfo_corrected_factorial_assessment import compare_execution, verify_checkout
from benchmark_tools.admit_qfo_recovered_assessment import validate_trace
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_qfo_parameter_assessment import (
    VARIANTS, ENV_SHA, CONVERSION_JOB, CONVERTER_COMMIT, completed_conversion, validate_stage,
)
from benchmark_tools.run_qfo_recovered_assessment import command_for, environment_records
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.validate_qfo_native_assessment import validate_directory

ASSESSMENT_JOB = "21944"
EXECUTOR_COMMIT = "f5a42ab3fc6ee5076e9dfb2b239b63b131aa0d73"


def completed_assessment(accounting, index):
    if type(index) is not int or index not in range(4):
        raise ValueError("Unknown parameter variant")
    rows = [r for r in csv.DictReader(io.StringIO(accounting), delimiter="|")
            if r["JobID"] == f"{ASSESSMENT_JOB}_{index}"]
    if len(rows) != 1 or tuple(rows[0][k] for k in ("State", "ExitCode", "NodeList", "AllocCPUS")) != (
            "COMPLETED", "0:0", "bizon", "8"):
        raise ValueError("Require successfully completed eight-CPU assessment")
    return rows[0]


def bind_records(report, base, executor):
    records = report["verified_records"]
    if records[:len(base)] != base:
        raise ValueError("Changed verified assessment inputs")
    helpers = records[len(base):]
    paths = [Path(r["path"]) for r in helpers]
    required = {"run_qfo_recovered_assessment.py", "run_simulation_methods.py",
                "prepare_qfo_corrected_comparator_pairs.py", "run_qfo_parameter_phylogeny.py",
                "prepare_ob_candidate_neighborhood.py"}
    if (not required.issubset({p.name for p in paths}) or len(paths) != len(set(paths))
            or any(p.parent != executor / "benchmark_tools" for p in paths)):
        raise ValueError("Missing, duplicate or foreign scoring helpers")
    for item in records:
        check(item)
    return records


def check_output_inventory(report, results):
    observed = [record(p) for p in sorted(results.rglob("*")) if p.is_file()]
    if observed != report["outputs"]:
        raise ValueError("Changed native assessment output inventory")
    traces = list((results / "stats").glob("trace_*.txt"))
    if len(traces) != 1:
        raise ValueError("Require exactly one native scoring trace")
    return observed, traces[0]


def admit(root, index, output):
    if type(index) is not int or index not in range(4):
        raise ValueError("Unknown parameter variant")
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    fmt = "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"
    accounting = subprocess.check_output(["sacct", "-j", ASSESSMENT_JOB, "--parsable2", fmt], text=True)
    scheduler = completed_assessment(accounting, index)
    source = record(__file__)
    helpers = [record(m.__file__) for n, m in sorted(sys.modules.items())
               if n.startswith("benchmark_tools.") and getattr(m, "__file__", None)]
    directory = root / "benchmarks/results/qfo_parameter_assessment_v1" / VARIANTS[index]
    report_path, preflight_path = directory / "results.json", directory / "preflight.json"
    checked = [record(report_path), record(preflight_path)]
    report, preflight = (json.loads(p.read_text()) for p in (report_path, preflight_path))
    validate_completion(report, preflight, scheduler)
    conversion_text = subprocess.check_output(["sacct", "-j", CONVERSION_JOB, "--parsable2", fmt], text=True)
    conversion_scheduler = completed_conversion(conversion_text, index)
    pairs_path = root / "benchmarks/results/qfo_parameter_pairs_v1" / VARIANTS[index] / "results.json"
    pair_record = record(pairs_path)
    stage = read_frozen(pairs_path, pair_record["sha256"])
    validate_stage(stage, index, conversion_scheduler)
    converter = root / "benchmarks/work/publication_qfo_parameter_pairs_v1"
    verify_checkout(converter, CONVERTER_COMMIT, ["benchmark_tools", "orthohmm", "qfo_benchmark/og_to_pairwise.py"])
    if record(converter / "benchmark_tools/prepare_qfo_parameter_pairs.py") not in stage["checked_records"]:
        raise ValueError("Wrong converter source")
    executor = root / "benchmarks/work/publication_qfo_parameter_assessment_v1"
    verify_checkout(executor, EXECUTOR_COMMIT, ["benchmark_tools", "orthohmm"])
    env_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment = read_frozen(env_path, ENV_SHA)
    if [r for r in environment["reference_files"] if Path(r["path"]).name == "mapping.json.gz"] != [stage["mapping"]]:
        raise ValueError("Conversion and scoring mapping differ")
    base = [pair_record, record(env_path), *environment_records(environment), *stage["checked_records"],
            stage["pairs"], stage["filtered_pairs"], stage["conversion_counts"]]
    records = bind_records(report, base, executor)
    counts = json.loads(Path(stage["conversion_counts"]["path"]).read_text())
    if counts != {k: stage[k] for k in ("written_pairs", "total_pairs", "retained_pairs", "removed_mapping_pairs")}:
        raise ValueError("Conversion count sidecar differs")
    work = root / "qfo_benchmark/w" / f"qpv{index}"
    results = root / "qfo_benchmark/scoring" / f"parameter_v1_{index}"
    expected = {"index": index, "variant": VARIANTS[index], "stage": stage,
        "source": record(executor / "benchmark_tools/run_qfo_parameter_assessment.py"),
        "pairs_manifest": pair_record, "environment_manifest": record(env_path),
        "conversion_scheduler": conversion_scheduler, "conversion_accounting": conversion_text,
        "converter_commit": CONVERTER_COMMIT, "command": command_for(root, stage, environment, work, results),
        "cwd": str(directory), "work": str(work), "results": str(results), "verified_records": records,
        "environment_overrides": environment["environment_overrides"], "array_task_id": str(index)}
    compare_execution(report, expected)
    if report["log"]["path"] != str(directory / "scoring.log"):
        raise ValueError("Unexpected scoring log")
    checked.extend([report["source"], report["log"], *records])
    for item in checked:
        check(item)
    observed, trace = check_output_inventory(report, results)
    tasks = validate_trace(trace.read_text())
    assessment, paths = validate_directory(results, stage["participant"], Path(environment["pipeline"]) / "reference_data")
    metrics = [record(p) for p in paths]
    for item in [source, *helpers, *checked, *observed, *metrics]:
        check(item)
    check_output_inventory(report, results)
    result = {"status": "corrected_parameter_assessment_admitted", "variant": VARIANTS[index], "index": index,
        "source": source, "validator_sources": helpers, "scheduler": scheduler, "accounting": accounting,
        "conversion_scheduler": conversion_scheduler, "conversion_accounting": conversion_text,
        "pairs_manifest": pair_record, "conversion": stage, "assessment": assessment,
        "execution_report": checked[0], "preflight": checked[1], "environment_manifest": record(env_path),
        "metric_files": metrics, "native_tasks": tasks, "native_trace": record(trace), "checked_records": checked,
        "accuracy_admitted": True, "publication_ready": False,
        "limitations": ["Development-exposed parameter robustness, not independent biological validation.",
            "Six-endpoint mean is a project-defined secondary summary, not official QfO F1.",
            "Paired uncertainty and prespecified multiplicity handling remain separate requirements."]}
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--index", required=True, type=int, choices=range(4))
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    admit(args.root.resolve(), args.index, args.output.resolve())
