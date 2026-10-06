"""Diagnose terminal collector cadence failures without repairing or admitting them."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.measure_native_hierarchy_step import interval_point
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_interval_cpu import MAX_GAP_S, interval, validate_point
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.replay_threadripper_scaling import point_inventory
from benchmark_tools.run_native_factorial_cost import ROOT, read, validate_plan, validate_request, verify_terminal
from benchmark_tools.validate_native_factorial_outputs import require


ERROR = "Missing or irregular observation interval"


def cadence(points, job):
    """Stream standard points; retain every rejection from the unchanged interval API."""
    previous = None
    failures, durations = [], []
    count = 0
    for point in points:
        validate_point(point, job)
        if previous is not None:
            duration = (sum(point["native_read_ns"]) - sum(previous["native_read_ns"])) / 2e9
            durations.append(duration)
            try:
                interval(previous, point, job)
            except ValueError as error:
                failures.append(dict(left_index=count - 1, right_index=count,
                    wall_s=duration, error=str(error), error_type=type(error).__name__,
                    left_native_read_ns=previous["native_read_ns"],
                    right_native_read_ns=point["native_read_ns"]))
        previous = point
        count += 1
    require(count >= 2, "Require at least two observation points")
    return dict(points=count, intervals=count - 1,
        minimum_period_s=min(durations), maximum_period_s=max(durations),
        unchanged_cadence_bounds_s=[.5, MAX_GAP_S], failures=failures,
        failure_counts=dict(Counter(item["error"] for item in failures)),
        full_resource_replay=False, scientific_timings_admitted=False,
        native_outputs_validated=False, next_identity_authorized=False)


def raw_points(directory, job, references, digest):
    paths = point_inventory(directory)
    for path in paths:
        require(path.resolve() == path and not path.is_symlink(), "Require direct raw points")
        before = path.stat()
        data = path.read_bytes()
        after = path.stat()
        require((before.st_ino, before.st_size, before.st_mtime_ns) ==
                (after.st_ino, after.st_size, after.st_mtime_ns), "Point changed during census")
        ref = dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
        references.append(ref)
        digest.update(json.dumps(ref, sort_keys=True, separators=(",", ":")).encode() + b"\n")
        point = json.loads(data)
        require(point.get("schema") == "native_lineage_v1", "Wrong raw observation schema")
        yield interval_point(point, job)
    require(point_inventory(directory) == paths, "Raw point inventory changed")


def audit(request_ref, failed_review_ref, destination):
    destination = Path(destination)
    require(destination.is_absolute() and destination.resolve() == destination
        and destination.is_relative_to(ROOT) and destination.parent.is_dir()
        and not destination.exists() and not destination.is_symlink(), "Require fresh direct output")
    request = read(request_ref)
    plan_ref = request["plan"]
    plan = read(plan_ref)
    run = validate_plan(plan)[request["index"]]
    validate_request(request, plan_ref, request["job_id"])
    failed_review = read(failed_review_ref)
    require(failed_review.get("status") == "terminal_factorial_review_failed"
        and failed_review.get("request") == request_ref and failed_review.get("plan") == plan_ref
        and failed_review.get("index") == run["index"]
        and failed_review.get("job_id") == request["job_id"]
        and failed_review.get("terminal_reviewed") is False
        and failed_review.get("next_identity_authorized") is False,
        "Require original failed terminal review")
    terminal = verify_terminal(request["job_id"])
    fields = terminal["verified"].get("fields", terminal["verified"])
    require(fields.get("JobState", fields.get("State")) == "FAILED" and fields["ExitCode"] == "1:0",
        "Do not relabel the enclosing scheduler failure")
    root = Path(run["output_root"])
    directory = root / "measurement"
    require(not destination.is_relative_to(root), "Keep diagnosis separate from original evidence")
    refs = [request_ref, plan_ref, failed_review_ref, record(__file__)]
    refs += [record(path) for path in (
        Path(plan["panel_root"]) / "sessions" / f"run_{run['index']:02d}" / "result.json",
        root / "verification.json", root / "native_execution.json", root / "metrics.json",
        directory / "done.json", directory / "report_finalization.json",
        directory / "process_stream_review.json", directory / "native_completion.json")]
    session, wrapper, native, metrics, done, finalization, environment, completion = [read(ref) for ref in refs[4:]]
    require(session.get("status") == "factorial_attempt_failed_retained"
        and session.get("request") == request_ref and session.get("plan") == plan_ref
        and session.get("index") == run["index"] and session.get("job_id") == request["job_id"]
        and session.get("wrapper") == wrapper, "Session and wrapper binding differ")
    require(wrapper.get("status") == "verified_wrapper_failed" and wrapper.get("error") == ERROR
        and finalization.get("status") == "reporting_failed" and finalization.get("error") == ERROR,
        "Not the retained cadence finalization failure")
    require(native.get("status") == "native_factorial_completed_pending_output_review"
        and native.get("index") == run["index"] and native.get("plan") == plan_ref
        and native.get("cell") == run["cell"] and native.get("counts") == metrics.get("counts")
        and metrics.get("status") == "complete" and type(done.get("exit_code")) is int
        and done["exit_code"] == 0 and done.get("timed_out") is False and not completion.get("errors"),
        "Incomplete native completion records; this audit is not output validation")
    missing = [name for name in ("lineage_report.json", "root_context_report.json")
               if not (directory / name).exists()]
    require("lineage_report.json" in missing, "Do not diagnose a repaired measurement")
    point_refs, digest = [], hashlib.sha256()
    census = cadence(raw_points(directory, request["job_id"], point_refs, digest), request["job_id"])
    require(census["failure_counts"].get(ERROR, 0) > 0, "Cadence failure did not reproduce")
    failed_points = sorted({i for item in census["failures"]
                            for i in (item["left_index"], item["right_index"])})
    for ref in [*refs, *plan["helper_sources"]]:
        check(ref)
    result = dict(schema="native_factorial_cadence_failure_audit_v1",
        status="retained_measurement_cadence_failure_reproduced", job_id=request["job_id"],
        index=run["index"], cell=run["cell"], plan=plan_ref, request=request_ref,
        scheduler=terminal, evidence=refs, census=census,
        raw_inventory=dict(points=len(point_refs), bytes=sum(r["bytes"] for r in point_refs),
            ordered_record_digest_sha256=digest.hexdigest(),
            encoding="SHA256 of sorted-key compact JSON path/bytes/sha256 records plus LF in point order",
            first=point_refs[0], last=point_refs[-1],
            failing_interval_points=[point_refs[i] for i in failed_points]),
        missing_measurement_reports=missing,
        native_completion_records=dict(status=native["status"], metric_status=metrics["status"],
            exit_code=done["exit_code"], timed_out=done["timed_out"]),
        retained_pressure_review=environment["pressure_review"],
        frozen_helper_sources_checked=len(plan["helper_sources"]),
        primary_resources=None, full_resource_replay=False, scientific_timings_admitted=False,
        native_outputs_validated=False, next_identity_authorized=False,
        original_receipts_rewritten=False, inference_reexecuted=False, automatic_retry=False,
        limitations=["Cadence diagnosis only, not full lineage/affinity/resource or scientific-output admission.",
            "Successful native completion records do not establish valid scientific outputs or benchmark scores.",
            "Original failed allocation and review remain failed; absent measurement reports are not synthesized.",
            "Unchanged interval API retains all rejections; no interpolation, threshold relaxation or timing correction.",
            "Shared-host contention effects and the cause of the cadence violations remain unknown."])
    save(destination, result)
    return record(destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("request", "failed-review", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("request-sha256", "failed-review-sha256"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    request_ref, review_ref = record(args.request), record(args.failed_review)
    require(request_ref["sha256"] == args.request_sha256
        and review_ref["sha256"] == args.failed_review_sha256, "Explicit evidence digest differs")
    print(json.dumps(audit(request_ref, review_ref, args.output.resolve()), sort_keys=True))


if __name__ == "__main__":
    main()
