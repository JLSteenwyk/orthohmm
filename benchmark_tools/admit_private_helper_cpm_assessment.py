"""Independently admit recovered private high-CPM native QfO scores."""

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import run_private_helper_cpm_assessment as runner
from benchmark_tools.admit_qfo_corrected_comparator_assessment import validate_completion
from benchmark_tools.admit_qfo_corrected_factorial_assessment import compare_execution, verify_checkout
from benchmark_tools.admit_qfo_parameter_assessment import check_output_inventory
from benchmark_tools.admit_qfo_recovered_assessment import validate_trace
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_qfo_recovered_assessment import command_for, environment_records
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.validate_qfo_native_assessment import validate_directory

EXECUTOR = "benchmarks/work/qfo_private_cpm_assessment_executor_20261001"
COMMIT = "903ea2d71a3b478b5198dd68dbbc088b84df4578"
RUNNER_SHA = "0a17e39f816331f2004d20aaa0dcfdc8f5033388a2062a6f3f813718e791a97e"
RUNNER_PROTOCOL_SHA = "873b5ae06d47c64b54bed520f82de2155b5af1c7a5e09f6dc83aecae3d54b031"
PROTOCOL = "benchmark_tools/results/QFO_PRIVATE_CPM_SCORE_ADMISSION_PROTOCOL_20261001.md"
REQUIRED_HELPERS = {"prepare_private_helper_cpm_pairs.py", "admit_private_helper_cpm_phylogeny.py",
    "run_qfo_recovered_assessment.py", "run_simulation_methods.py", "prepare_ob_candidate_neighborhood.py",
    "prepare_qfo_corrected_comparator_pairs.py", "prepare_qfo_parameter_pairs.py"}


def assessment_job(job):
    if type(job) is not str or not job.isascii() or not job.isdecimal() or job.startswith("0"):
        raise ValueError("Require explicit standalone assessment job identity")


def completed_assessment(text, job):
    assessment_job(job)
    rows = list(csv.DictReader(io.StringIO(text), delimiter="|"))
    parents = [row for row in rows if row["JobID"] == job]
    if len(parents) != 1 or tuple(parents[0][key] for key in
            ("JobIDRaw", "State", "ExitCode", "NodeList", "AllocCPUS", "ReqMem")) != (
                job, "COMPLETED", "0:0", "bizon", "8", "64G"):
        raise ValueError("Require successfully completed private eight-CPU assessment")
    steps = [row for row in rows if row["JobID"].startswith(job + ".")]
    if (len(rows) != 1 + len(steps) or len({row["JobID"] for row in steps}) != len(steps)
            or not any(row["JobID"] == job + ".batch" for row in steps)
            or any(row["State"] != "COMPLETED" or row["ExitCode"] != "0:0" for row in steps)):
        raise ValueError("Require successful unique terminal assessment steps")
    return parents[0]


def git_record(executor, ref):
    path = Path(ref["path"])
    if path.resolve() != path or not path.is_relative_to(executor):
        raise ValueError("Foreign or indirect submitted assessment source")
    blob = subprocess.check_output(["git", "-C", str(executor), "show", COMMIT + ":" + str(path.relative_to(executor))])
    if len(blob) != ref["bytes"] or hashlib.sha256(blob).hexdigest() != ref["sha256"]:
        raise ValueError("Submitted assessment source differs from frozen Git blob")
    check(ref)


def submitted_assessment(root, job, sha, verified):
    path = root / f"benchmark_tools/results/qfo_private_cpm_assessment_submission_{job}.json"
    receipt = read_frozen(path, sha)
    executor = root / EXECUTOR
    script = "benchmark_tools/results/qfo_private_cpm_assessment_20261001.sh"
    names = ("benchmark_tools/run_private_helper_cpm_assessment.py", runner.PROTOCOL,
             script, "tests/unit/test_run_private_helper_cpm_assessment.py")
    sources = [record(executor / name) for name in names]
    if (receipt["status"] != "private_recovered_qfo_high_cpm_assessment_submitted"
            or receipt["job_id"] != job or receipt["executor"] != str(executor)
            or receipt["executor_commit"] != COMMIT or receipt["executor_clean"] is not True
            or any(receipt[key] is not False for key in
                   ("native_inference", "accuracy_evaluated", "controlled_timing", "publication_ready"))
            or receipt["source_records"] != sources
            or receipt["pair_conversion"] != verified["pairs_manifest"]
            or receipt["parent_submission"] != verified["submission"]
            or receipt["submission_argv"] != ["sbatch", "--parsable", str(executor / script), str(executor), COMMIT,
                RUNNER_PROTOCOL_SHA, verified["scheduler"]["JobIDRaw"], verified["pairs_manifest"]["sha256"],
                verified["submission"]["sha256"]]):
        raise ValueError("Private assessment actual submission binding differs")
    verify_checkout(executor, COMMIT, ["benchmark_tools", "orthohmm", names[-1]])
    if subprocess.check_output(["git", "-C", str(executor), "status", "--porcelain", "--untracked-files=all"], text=True).strip():
        raise ValueError("Private assessment executor is not clean")
    for ref in sources:
        git_record(executor, ref)
    if (sources[0]["sha256"] != RUNNER_SHA or record(runner.__file__)["sha256"] != RUNNER_SHA
            or sources[1]["sha256"] != RUNNER_PROTOCOL_SHA):
        raise ValueError("Private assessment source or reused verifier changed")
    return sources[0], [record(path), *sources]


def bind_records(report, base, executor):
    helpers = report["helpers"]
    paths = [Path(ref["path"]) for ref in helpers]
    if (report["verified_records"] != [*base, *helpers]
            or not REQUIRED_HELPERS.issubset({path.name for path in paths})
            or len(paths) != len(set(paths))
            or any(path.parent != executor / "benchmark_tools" for path in paths)):
        raise ValueError("Missing, duplicate, foreign or changed private assessment helper bindings")
    for ref in helpers:
        git_record(executor, ref)
        local = Path(__file__).with_name(Path(ref["path"]).name)
        if record(local)["sha256"] != ref["sha256"]:
            raise ValueError("Reused scoring helper differs from submitted source")
    for ref in [*base, *helpers]:
        check(ref)
    return [*base, *helpers]


def admit(root, job, report_sha, submission_sha, protocol_sha, output):
    if not output.is_absolute() or output.resolve() != output:
        raise ValueError("Require direct absolute private score-admission destination")
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    assessment_job(job)
    text = runner.converter.accounting(job)
    scheduler = completed_assessment(text, job)
    source = record(__file__)
    protocol = record(root / PROTOCOL)
    if protocol["sha256"] != protocol_sha:
        raise ValueError("Private score-admission protocol changed")
    helpers = [record(module.__file__) for name, module in sorted(sys.modules.items())
        if name.startswith("benchmark_tools.") and getattr(module, "__file__", None)]
    checked = [source, protocol, *helpers]
    result = dict(status="validating", source=source, protocol=protocol, validator_sources=helpers,
        scheduler=scheduler, accounting=text, checked_records=checked, accuracy_admitted=False,
        controlled_timing=False, publication_ready=False)
    try:
        directory = root / runner.OUTPUT
        report_path, preflight_path = directory / "results.json", directory / "preflight.json"
        report = read_frozen(report_path, report_sha)
        checked.extend([record(report_path), record(preflight_path)])
        preflight = read_frozen(preflight_path, checked[-1]["sha256"])
        validate_completion(report, preflight, scheduler)
        if (report["controlled_timing"] is not False or report["publication_ready"] is not False
                or type(report["index"]) is not int or report["index"] != 1 or report["arm"] != "cpm_high"
                or type(report["exit_code"]) is not int):
            raise ValueError("Private assessment arm or claim flags differ")
        verified = runner.verify_conversion(root, report["conversion_job"], report["conversion_sha256"],
                                            report["conversion_submission_sha256"])
        executor_source, submitted = submitted_assessment(root, job, submission_sha, verified)
        stage = verified["stage"]
        env_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
        environment = read_frozen(env_path, runner.converter.ENV_SHA)
        if [ref for ref in environment["reference_files"] if Path(ref["path"]).name == "mapping.json.gz"] != [stage["mapping"]]:
            raise ValueError("Private scoring and conversion mapping differ")
        scoring_protocol = record(root / runner.PROTOCOL)
        if scoring_protocol["sha256"] != RUNNER_PROTOCOL_SHA:
            raise ValueError("Frozen private scoring protocol changed")
        base = [scoring_protocol, verified["pairs_manifest"], verified["submission"], record(env_path),
                *environment_records(environment), *verified["checked_records"]]
        records = bind_records(report, base, root / EXECUTOR)
        work, results = root / runner.WORK, root / runner.RESULTS
        expected = dict(arm="cpm_high", index=1, context=verified["context"], stage=stage,
            source=executor_source, protocol=scoring_protocol, pairs_manifest=verified["pairs_manifest"],
            conversion_submission=verified["submission"], environment_manifest=record(env_path),
            conversion_scheduler=verified["scheduler"], conversion_accounting=verified["accounting"],
            conversion_job=verified["scheduler"]["JobIDRaw"], conversion_sha256=verified["pairs_manifest"]["sha256"],
            conversion_submission_sha256=verified["submission"]["sha256"], converter_commit=runner.CONVERTER_COMMIT,
            command=command_for(root, stage, environment, work, results), cwd=str(directory), work=str(work),
            results=str(results), verified_records=records, environment_overrides=environment["environment_overrides"],
            accuracy_admitted=False, controlled_timing=False, publication_ready=False)
        compare_execution(report, expected)
        if report["log"] != record(directory / "scoring.log"):
            raise ValueError("Unexpected private scoring log")
        checked.extend([*submitted, report["source"], report["log"], *records])
        for ref in checked:
            check(ref)
        observed, trace = check_output_inventory(report, results)
        tasks = validate_trace(trace.read_text())
        assessment, paths = validate_directory(results, stage["participant"], Path(environment["pipeline"]) / "reference_data")
        metrics = [record(path) for path in paths]
        available = {ref["path"] for ref in [*observed, *records]}
        if any(ref["path"] not in available for ref in metrics) or record(trace) not in observed:
            raise ValueError("Uninventoried metric, reference or native trace")
        checked.extend([*observed, *metrics])
        if runner.verify_conversion(root, report["conversion_job"], report["conversion_sha256"],
                                    report["conversion_submission_sha256"]) != verified:
            raise ValueError("Private conversion/native context changed during score admission")
        if completed_assessment(runner.converter.accounting(job), job) != scheduler:
            raise ValueError("Private scoring completion changed during admission")
        if submitted_assessment(root, job, submission_sha, verified)[1] != submitted:
            raise ValueError("Private assessment submission changed during admission")
        for ref in checked:
            check(ref)
        check_output_inventory(report, results)
        result.update(status="private_recovered_cpm_assessment_admitted", arm="cpm_high", index=1,
            context=verified["context"], pairs_manifest=verified["pairs_manifest"], conversion=stage,
            conversion_scheduler=verified["scheduler"], conversion_accounting=verified["accounting"],
            execution_report=record(report_path), preflight=record(preflight_path), submission=submitted[0],
            environment_manifest=record(env_path), assessment=assessment, native_tasks=tasks,
            native_trace=record(trace), metric_files=metrics, accuracy_admitted=True,
            limitations=["Development-exposed recovered parameter arm, not independent biological validation.",
                "Shared-host incremental inference and scoring are not controlled timing evidence.",
                "Six-metric mean is a project-defined secondary summary, not official QfO F1.",
                "Native FAS sampling remains unseeded; sample attrition and dependence require separate analysis.",
                "Paired uncertainty and frozen seven-arm/18-endpoint multiplicity remain separate requirements."])
    except BaseException as error:
        result.update(status="validation_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with output.open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("job", "report-sha256", "submission-sha256", "protocol-sha256"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.job, args.report_sha256, args.submission_sha256, args.protocol_sha256, args.output)
