"""Independently validate canonical assessment execution and native endpoints."""

import argparse
import json
from pathlib import Path
import subprocess

from benchmark_tools.run_qfo_order_replay import record, save
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.prepare_qfo_canonical_pairs import check_records
from benchmark_tools.run_qfo_canonical_assessment import verify_stage, VERIFIED_SHA
from benchmark_tools.run_qfo_recovered_assessment import command_for, environment_records
from benchmark_tools.prepare_qfo_corrected_comparator_pairs import ENV_SHA
from benchmark_tools.admit_qfo_recovered_assessment import validate_trace
from benchmark_tools.validate_qfo_native_assessment import validate_directory
from benchmark_tools.verify_ygob_validation import require_completed_job

PLAN_SHA = "41910a6cfc113c028bbcbe0bed2bb9d8ca20c05cf8bb530d1868e31c79ca499b"
SUBMISSION_SHA = "96adc6fd4b6a07dbc16c8d2b31dae27f0927ef78fb11edf76b5800d7613ca9d3"


def verify_completion(report, preflight, plan_record, command, source):
    expected = dict(plan=plan_record, status="running", job_id="22336",
                    source=source, accuracy_admitted=False, command=command)
    if preflight != expected:
        raise ValueError("Wrong assessment preflight")
    if any(report.get(k) != v for k, v in expected.items() if k != "status"):
        raise ValueError("Changed assessment provenance")
    if report.get("status") != "process_succeeded_pending_independent_admission" or report.get("exit_code") != 0:
        raise ValueError("Assessment not successful")


def admit(repo, output):
    if output.exists():
        raise FileExistsError(output)
    accounting = subprocess.check_output(["sacct", "-j", "22336", "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, 22336)
    if scheduler["AllocCPUS"] != "8":
        raise ValueError("Wrong scoring CPU allocation")
    directory = repo / "benchmarks/work/qfo_canonical_assessment_20260927"
    plan_path = directory / "plan.json"
    plan = read_frozen(plan_path, PLAN_SHA)
    submission_path = directory / "submission.json"
    submission = read_frozen(submission_path, SUBMISSION_SHA)
    if submission["job_id"] != "22336" or submission["plan"] != record(plan_path):
        raise ValueError("Wrong scoring submission")
    check_records(plan["checked_records"])
    verified = read_frozen(Path(plan["conversion_verification"]["path"]), VERIFIED_SHA)
    execution = json.loads(Path(verified["execution"]["path"]).read_text())
    verify_stage(plan["stage"], verified, execution)
    environment = read_frozen(Path(plan["environment_manifest"]["path"]), ENV_SHA)
    check_records(environment_records(environment))
    results = repo / "qfo_benchmark/scoring/canonical_20260927"
    work = repo / "qfo_benchmark/w/qcan27"
    execution_dir = repo / "benchmarks/results/qfo_canonical_assessment_20260927"
    if (plan["repo"] != str(repo) or plan["directory"] != str(directory)
            or plan["results"] != str(results) or plan["work"] != str(work)
            or plan["output"] != str(execution_dir)
            or plan["environment_overrides"] != environment["environment_overrides"]):
        raise ValueError("Wrong scoring namespace or environment")
    command = command_for(repo, plan["stage"], environment, work, results)
    if plan["command"] != command:
        raise ValueError("Changed scoring command")
    report_path, preflight_path = execution_dir / "results.json", execution_dir / "preflight.json"
    report, preflight = [json.loads(p.read_text()) for p in (report_path, preflight_path)]
    verify_completion(report, preflight, record(plan_path), command,
                      record(repo / "benchmark_tools/run_qfo_canonical_assessment.py"))
    if report["log"] != record(execution_dir / "scoring.log"):
        raise ValueError("Wrong or changed scoring log")
    observed = [record(p) for p in sorted(results.rglob("*")) if p.is_file()]
    if observed != report["outputs"]:
        raise ValueError("Changed scoring output inventory")
    traces = list((results / "stats").glob("trace_*.txt"))
    if len(traces) != 1:
        raise ValueError("Require one scoring trace")
    tasks = validate_trace(traces[0].read_text())
    assessment, paths = validate_directory(results, plan["stage"]["participant"],
                                           Path(environment["pipeline"]) / "reference_data")
    checked = [*plan["checked_records"], *observed, *[record(p) for p in paths],
               *[record(p) for p in (plan_path, submission_path, report_path, preflight_path, Path(__file__))]]
    check_records(checked)
    result = dict(status="canonical_qfo_assessment_admitted", scheduler=scheduler,
        plan=record(plan_path), submission=record(submission_path), execution_report=record(report_path),
        conversion=plan["conversion_verification"], assessment=assessment, native_tasks=tasks,
        native_trace=record(traces[0]), metric_files=[record(p) for p in paths], checked_records=checked,
        source=record(__file__), accuracy_admitted=True, publication_ready=False,
        limitations=["Development-exposed downstream ordering experiment, not independent generalization",
            "FAS is an unseeded native sample; differences are not solely attributable to ordering",
            "Six-endpoint mean is a project-defined secondary summary, not F1",
            "Native standard errors are not paired difference confidence intervals"])
    save(output, result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = admit(args.repo.resolve(), args.output.resolve())
    print(json.dumps(dict(status=result["status"], output=record(args.output))))
