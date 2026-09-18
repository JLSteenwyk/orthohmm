"""Validate terminal sequence-control QfO endpoint evidence before admitting scores."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_corrected_comparator_assessment import validate_completion
from benchmark_tools.admit_qfo_corrected_factorial_assessment import compare_execution
from benchmark_tools.admit_qfo_recovered_assessment import validate_trace
from benchmark_tools.admit_qfo_sequence_numeric import completed
from benchmark_tools.compare_qfo_search_coverage import frozen
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_qfo_sequence_assessment import prepare, WORK_NAMES, ENV_SHA
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.validate_qfo_native_assessment import validate_directory

EXECUTOR = "425f0a7f5d5dc9e1438ab0a1766c45b13596dc14"


def outputs(results, expected):
    observed = [record(path) for path in sorted(results.rglob("*")) if path.is_file()]
    if not observed or observed != expected:
        raise ValueError("Native output inventory/content differs")
    traces = list((results / "stats").glob("trace_*.txt"))
    if len(traces) != 1:
        raise ValueError("Require exactly one native scoring trace")
    return observed, traces[0]


def admit(root, label, job, conversion_job, pairs_sha, output):
    if label not in WORK_NAMES:
        raise ValueError("Unknown sequence variant")
    if output.exists():
        raise FileExistsError(output)
    scheduler = completed(job, 8, "64G")
    executor = frozen(root, "publication_qfo_sequence_assessment_v1", EXECUTOR)
    # Rebuild expected provenance only with helpers identical to the frozen runner.
    for name in ("run_qfo_sequence_assessment.py", "run_qfo_recovered_assessment.py"):
        current, original = record(Path(__file__).with_name(name)), record(executor / "benchmark_tools" / name)
        if any(current[key] != original[key] for key in ("sha256", "bytes")):
            raise ValueError("Assessment reconstruction helper differs from frozen runner")
    expected = prepare(root, label, pairs_sha, conversion_job, require_fresh=False)
    expected["source"] = record(executor / "benchmark_tools/run_qfo_sequence_assessment.py")
    own_helper = record(Path(__file__).with_name("run_qfo_recovered_assessment.py"))
    expected["verified_records"] = [record(executor / "benchmark_tools/run_qfo_recovered_assessment.py")
        if item == own_helper else item for item in expected["verified_records"]]
    directory = Path(expected["cwd"])
    report_path, preflight_path = directory / "results.json", directory / "preflight.json"
    checked = [record(report_path), record(preflight_path)]
    report, preflight = (json.loads(path.read_text()) for path in (report_path, preflight_path))
    validate_completion(report, preflight, scheduler)
    compare_execution(report, {key: value for key, value in expected.items() if key != "status"})
    if report["log"]["path"] != str(directory / "scoring.log"):
        raise ValueError("Unexpected scoring log location")
    checked.extend([report["source"], report["log"], *expected["verified_records"]])
    for item in checked:
        check(item)
    results = Path(report["results"])
    observed, trace = outputs(results, report["outputs"])
    tasks = validate_trace(trace.read_text())
    environment = read_frozen(Path(report["environment_manifest"]["path"]), ENV_SHA)
    assessment, paths = validate_directory(results, report["stage"]["participant"],
                                           Path(environment["pipeline"]) / "reference_data")
    metrics = [record(path) for path in paths]
    for item in [*checked, *observed, *metrics]:
        check(item)
    result = dict(status="corrected_sequence_assessment_admitted", variant=label, source=record(__file__),
        scheduler=scheduler, conversion_scheduler=report["conversion_scheduler"],
        pairs_manifest=report["pairs_manifest"], conversion=report["stage"], assessment=assessment,
        execution_report=record(report_path), preflight=record(preflight_path),
        environment_manifest=report["environment_manifest"], metric_files=metrics, native_tasks=tasks,
        native_trace=record(trace), checked_records=checked, accuracy_admitted=True, publication_ready=False,
        validator_sources=[record(Path(__file__).with_name(name)) for name in (
            "admit_qfo_sequence_assessment.py", "run_qfo_sequence_assessment.py",
            "admit_qfo_corrected_comparator_assessment.py", "admit_qfo_corrected_factorial_assessment.py",
            "admit_qfo_recovered_assessment.py", "validate_qfo_native_assessment.py")],
        limitations=["Native score/provenance validation, not independent biological confirmation or paired uncertainty.",
                     "The six-endpoint mean is a project-defined secondary summary, not official QfO F1.",
                     "Group-derived clique pairs are not native phylogenetic pair predictions.",
                     "Equal search cutoffs do not establish matched sensitivity or computational effort."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("job", "conversion-job", "pairs-sha256"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--variant", choices=WORK_NAMES, required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.variant, args.job, args.conversion_job, args.pairs_sha256, args.output.resolve())
