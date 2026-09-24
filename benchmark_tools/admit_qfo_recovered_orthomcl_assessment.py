"""Independently admit six QfO endpoints from the recovered OrthoMCL workflow."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_corrected_comparator_assessment import validate_completion
from benchmark_tools.admit_qfo_corrected_orthomcl import completed, frozen, unique_records
from benchmark_tools.admit_qfo_recovered_assessment import validate_trace
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_qfo_recovered_orthomcl_assessment import prepare, ENV_SHA
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.validate_qfo_native_assessment import validate_directory

RUNNER_SHA = "05b58b4e4479ecf4550feef395b9abd71f093e0e488886b8266d3126075db392"


def admit(root, job, conversion_job, digest, executor, commit, converter, converter_commit, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if not executor.resolve().is_relative_to(root / "benchmarks/work"):
        raise ValueError("Require retained assessment executor")
    scheduler, accounting = completed(job, 8, "64G")
    frozen(executor, commit)
    source = record(executor / "benchmark_tools/run_qfo_recovered_orthomcl_assessment.py")
    if source["sha256"] != RUNNER_SHA:
        raise ValueError("Unreviewed recovered assessment runner")
    directory = root / "benchmarks/results/qfo_blast_recovery_assessment_v1"
    report_record, preflight_record = record(directory / "results.json"), record(directory / "preflight.json")
    report = read_frozen(Path(report_record["path"]), report_record["sha256"])
    preflight = read_frozen(Path(preflight_record["path"]), preflight_record["sha256"])
    validate_completion(report, preflight, scheduler)
    if type(report["exit_code"]) is not int or report["publication_ready"] is not False:
        raise ValueError("Invalid recovered assessment completion flags")
    expected = prepare(root, digest, conversion_job, converter, converter_commit,
                       require_fresh=False, helpers=executor / "benchmark_tools")
    for key, value in expected.items():
        if key != "status" and report.get(key) != value:
            raise ValueError("Changed recovered assessment provenance: " + key)
    if report["source"] != source or report["log"]["path"] != str(directory / "scoring.log"):
        raise ValueError("Wrong recovered assessment source/log")
    checked = unique_records([report_record, preflight_record, source, report["log"],
                              *expected["verified_records"]])
    for item in checked:
        check(item)
    results = Path(expected["results"])
    observed = [record(p) for p in sorted(results.rglob("*")) if p.is_file()]
    if observed != report["outputs"]:
        raise ValueError("Changed recovered scoring output inventory")
    traces = list((results / "stats").glob("trace_*.txt"))
    if len(traces) != 1:
        raise ValueError("Require one recovered native scoring trace")
    tasks = validate_trace(traces[0].read_text())
    manifest = read_frozen(Path(expected["environment_manifest"]["path"]), ENV_SHA)
    assessment, paths = validate_directory(results, expected["stage"]["participant"],
                                            Path(manifest["pipeline"]) / "reference_data")
    metrics = [record(path) for path in paths]
    for item in [*checked, *observed, *metrics]:
        check(item)
    stage = expected["stage"]
    result = dict(status="recovered_orthomcl_assessment_admitted", method="orthomcl", source=record(__file__),
        scheduler=scheduler, accounting=accounting, conversion_scheduler=expected["conversion_scheduler"],
        pairs_manifest=expected["pairs_manifest"], execution_report=report_record, preflight=preflight_record,
        assessment=assessment, metric_files=metrics, native_tasks=tasks, native_trace=record(traces[0]),
        checked_records=checked, query_coverage=stage["query_coverage"], group_coverage=stage["content"],
        pair_semantics=stage["semantics"], group_audit=stage["group_audit"],
        validator_sources=[record(Path(__file__).with_name(name)) for name in (
            "run_qfo_recovered_orthomcl_assessment.py", "admit_qfo_corrected_comparator_assessment.py",
            "admit_qfo_recovered_assessment.py", "validate_qfo_native_assessment.py")],
        accuracy_admitted=True, publication_ready=False, limitations=[
            "Native score/provenance validation, not independent biological validation or paired uncertainty.",
            "The six-endpoint mean is a project-defined secondary summary, not official QfO F1.",
            "Final-group cliques retain ungrouped inputs and failed queries; no missing search hits are repaired.",
            "Native error fields retain endpoint-specific meanings; shared-host timings are uncontrolled."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("root", "executor", "converter", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    for flag in ("job", "conversion-job"):
        parser.add_argument("--" + flag, type=int, required=True)
    for flag in ("commit", "converter-commit", "pairs-sha256"):
        parser.add_argument("--" + flag, required=True)
    args = parser.parse_args()
    print(admit(args.root.resolve(), args.job, args.conversion_job, args.pairs_sha256,
                args.executor.resolve(), args.commit, args.converter.resolve(), args.converter_commit,
                args.output.absolute())["status"])
