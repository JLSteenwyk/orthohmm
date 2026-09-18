"""Independently admit completed corrected factorial QfO endpoint evidence."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_corrected_comparator_assessment import accounting, validate_completion
from benchmark_tools.admit_qfo_recovered_assessment import validate_trace
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_qfo_corrected_factorial_assessment import CELLS, CONVERTERS, ENV_SHA, validate_stage
from benchmark_tools.run_qfo_recovered_assessment import command_for, environment_records
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.validate_qfo_native_assessment import validate_directory

EXECUTOR = "97bc99beafeb93a212aad191c589f8a7746bcbe3"


def verify_checkout(path, commit, paths):
    if subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"], text=True).strip() != commit:
        raise ValueError("Frozen executor revision changed")
    subprocess.run(["git", "-C", str(path), "diff", "--exit-code", "HEAD", "--", *paths], check=True)


def compare_execution(report, expected):
    for key, value in expected.items():
        if report.get(key) != value:
            raise ValueError("Changed execution provenance: " + key)


def admit(root, index, job, conversion_job, pairs_sha, output):
    if type(index) is not int or not 0 <= index < 8:
        raise ValueError("Unknown corrected factorial cell")
    if output.exists():
        raise FileExistsError(output)
    text, scheduler = accounting(job)
    directory = root / "benchmarks/results/qfo_corrected_factorial_assessment_v1" / CELLS[index]
    report_path, preflight_path = directory / "results.json", directory / "preflight.json"
    checked = [record(report_path), record(preflight_path)]
    report, preflight = (json.loads(p.read_text()) for p in (report_path, preflight_path))
    validate_completion(report, preflight, scheduler)
    pairs_path = root / "benchmarks/results/qfo_corrected_factorial_pairs_v1" / CELLS[index] / "results.json"
    stage = read_frozen(pairs_path, pairs_sha)
    conversion_text, conversion_scheduler = accounting(conversion_job)
    validate_stage(stage, index, conversion_scheduler)
    kind, converter_commit, _ = CONVERTERS[index % 2]
    converter = root / f"benchmarks/work/publication_qfo_corrected_{kind}_pairs_v1"
    verify_checkout(converter, converter_commit, ["benchmark_tools", "orthohmm", "qfo_benchmark/og_to_pairwise.py"])
    converter_source = record(converter / f"benchmark_tools/prepare_qfo_corrected_{kind}_pairs.py")
    if converter_source not in stage["checked_records"]:
        raise ValueError("Wrong frozen converter source")
    env_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    manifest = read_frozen(env_path, ENV_SHA)
    mappings = [r for r in manifest["reference_files"] if Path(r["path"]).name == "mapping.json.gz"]
    if mappings != [stage["mapping"]]:
        raise ValueError("Reference mapping differs from converted predictions")
    executor = root / "benchmarks/work/publication_qfo_corrected_factorial_assessment_v1"
    verify_checkout(executor, EXECUTOR, ["benchmark_tools"])
    work = root / "qfo_benchmark/w" / f"qcf{index}"
    results = root / "qfo_benchmark/scoring" / f"corrected_factorial_{index}"
    records = [record(pairs_path), record(env_path), *environment_records(manifest), *stage["checked_records"],
               stage["pairs"], stage["filtered_pairs"], record(executor / "benchmark_tools/run_qfo_recovered_assessment.py")]
    expected = {"index": index, "cell": CELLS[index], "stage": stage,
        "source": record(executor / "benchmark_tools/run_qfo_corrected_factorial_assessment.py"),
        "pairs_manifest": record(pairs_path), "environment_manifest": record(env_path),
        "conversion_scheduler": conversion_scheduler, "conversion_accounting": conversion_text,
        "converter_commit": converter_commit, "command": command_for(root, stage, manifest, work, results),
        "cwd": str(directory), "work": str(work), "results": str(results), "verified_records": records,
        "environment_overrides": manifest["environment_overrides"]}
    compare_execution(report, expected)
    if report["log"]["path"] != str(directory / "scoring.log"):
        raise ValueError("Unexpected assessment log path")
    checked.extend([report["source"], report["log"], *records])
    for item in checked:
        check(item)
    observed = [record(p) for p in sorted(results.rglob("*")) if p.is_file()]
    if observed != report["outputs"]:
        raise ValueError("Changed native output set or content")
    traces = list((results / "stats").glob("trace_*.txt"))
    if len(traces) != 1:
        raise ValueError("Require exactly one native scoring trace")
    tasks = validate_trace(traces[0].read_text())
    assessment, paths = validate_directory(results, stage["participant"], Path(manifest["pipeline"]) / "reference_data")
    metrics = [record(p) for p in paths]
    for item in [*checked, *observed, *metrics]:
        check(item)
    result = {"status": "corrected_factorial_assessment_admitted", "index": index, "cell": CELLS[index],
        "source": record(__file__), "scheduler": scheduler, "accounting": text,
        "conversion_scheduler": conversion_scheduler, "conversion_accounting": conversion_text,
        "pairs_manifest": record(pairs_path), "conversion": stage, "assessment": assessment,
        "execution_report": record(report_path), "preflight": record(preflight_path),
        "environment_manifest": record(env_path), "metric_files": metrics, "native_tasks": tasks,
        "native_trace": record(traces[0]), "checked_records": checked,
        "validator_sources": [record(Path(__file__).with_name(name)) for name in (
            "run_qfo_corrected_factorial_assessment.py", "admit_qfo_corrected_comparator_assessment.py",
            "run_qfo_recovered_assessment.py", "admit_qfo_recovered_assessment.py", "validate_qfo_native_assessment.py")],
        "accuracy_admitted": True, "publication_ready": False,
        "limitations": ["Native score/provenance validation, not independent biological validation or paired uncertainty.",
            "The six-endpoint mean is a project-defined secondary summary, not official QfO F1.",
            "Original-release scores and group-derived versus native-pair semantics are not interchangeable."]}
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--index", type=int, choices=range(8), required=True)
    for name in ("job", "conversion-job", "pairs-sha256"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.index, args.job, args.conversion_job, args.pairs_sha256, args.output.resolve())
