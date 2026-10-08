"""Replay native11 runtime brackets and exact late additions, without admission."""

import argparse
import copy
from datetime import datetime
import json
from pathlib import Path
import time
from zoneinfo import ZoneInfo

from benchmark_tools.classify_postterminal_runtime_additions import classify
from benchmark_tools.native_factorial_allocated_execution import amendment, validate_request
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.review_allocated_native_factorial_attempt import bind_session
from benchmark_tools.review_native_factorial_attempt import runtime_review, scheduler_fields
from benchmark_tools.run_native_factorial_cost import ROOT, read, verify_terminal
from benchmark_tools.snapshot_runtime_trees import inventory
from benchmark_tools.validate_native_factorial_outputs import Evidence, require


SCHEMA = "native11_postterminal_runtime_component_v1"
REQUEST = ROOT / "benchmarks/work/native_factorial_launch_20261004/request_11_allocated_v1.json"
REQUEST_SHA = "7bf63b80bd5932b9edbd1b2c5ff3fb77f5557e4f6c64e077045d6a50c8d366a1"
CLASSIFICATION = ROOT / "benchmark_tools/results/native11_postterminal_additions_classification_20261008_v1.json"
CLASSIFICATION_SHA = "8718376a789f2175305d83496984ad057aaa559b57469638e97f05ec1fd9537b"
DEFAULT_DESTINATION = ROOT / "benchmarks/work/native11_postterminal_runtime_review_20261008_v1"
SOURCES = {
    "benchmark_tools/review_native_factorial_attempt.py":
        "63e7d7fdda52afa7a36eecd89492a2684e8310260febbfb5adf6d23f88fd26de",
    "benchmark_tools/snapshot_runtime_trees.py":
        "ab0008651e75ebadc1be3c87503f09bffd9e5be87f2e6b0c3eec49e552d9bd4c",
    "benchmark_tools/classify_postterminal_runtime_additions.py":
        "782defc4f5103511460658c44ad31fbe87c01d97b0506ecc784580c35b43751e",
}


def compare_inventory(expected, observed):
    def indexed(value):
        rows = value["records"]
        require(isinstance(rows, list) and bool(rows), "Require nonempty inventory records")
        by_path = {row["path"]: row for row in rows}
        require(len(by_path) == len(rows), "Duplicate inventory path")
        return by_path

    old, new = indexed(expected), indexed(observed)
    metadata_old = {k: v for k, v in expected.items() if k != "records"}
    metadata_new = {k: v for k, v in observed.items() if k != "records"}
    return dict(equal=expected == observed, expected_records=len(old), observed_records=len(new),
        added=[new[p] for p in sorted(new.keys() - old.keys())],
        missing=[old[p] for p in sorted(old.keys() - new.keys())],
        changed=[dict(expected=old[p], observed=new[p]) for p in sorted(old.keys() & new.keys())
                 if old[p] != new[p]],
        metadata_differences={k: dict(expected=metadata_old.get(k), observed=metadata_new.get(k))
            for k in sorted(metadata_old.keys() | metadata_new.keys())
            if k not in metadata_old or k not in metadata_new or metadata_old[k] != metadata_new[k]})


def historical_checker(specs, prior):
    require(prior.get("status") == "runtime_brackets_and_lookup_replayed"
        and prior.get("continuous_runtime_integrity_established") is False,
        "Require retained historical runtime report")
    descriptors = prior["fresh_terminal_inventory"]
    require(len(descriptors) == len(specs) and len(specs) == 2,
        "Historical runtime inventory count differs")
    for (path, digest), descriptor in zip(specs, descriptors):
        manifest = read(dict(record(path), sha256=digest))
        require(descriptor == dict(path=path, sha256=digest,
            status="runtime_tree_identity_matches", records=len(manifest["records"]),
            scientific_execution_authorized=False), "Historical runtime descriptor differs")
    return copy.deepcopy(descriptors)


def chronology(classification, failure):
    require(failure.get("schema") == "native11_full_review_runtime_refusal_observation_v1"
        and failure.get("native_job_id") == 23985
        and failure.get("original_failed_review_job_id") == 23986
        and failure.get("failed_review_job_id") == 24033
        and failure.get("original_reviewer_error") == "Runtime inventory changed",
        "Wrong retained failure observation")
    require(classification.get("prior_inventory_time_is_review_terminal_upper_bound") is True
        and classification.get("timezone") == "America/New_York",
        "Require explicit historical-check upper bound")
    rows = [line.split("|") for line in failure["accounting"]["stdout"].splitlines()]
    require(len(rows) == 3 and all(len(row) == 12 for row in rows)
        and {row[0] for row in rows} == {"23985", "23986", "24033"},
        "Retained accounting identity or fields differ")
    jobs = {row[0]: row for row in rows}
    require(jobs["23985"][1:3] == ["COMPLETED", "0:0"]
        and jobs["23986"][1:3] == ["FAILED", "0:11"]
        and jobs["24033"][1:3] == ["FAILED", "1:0"], "Retained terminal outcomes differ")
    line = classification["package_install_matching_line"]
    require(line in failure["package_log_matching_lines"]["/var/log/dpkg.log"]
        and line.endswith(" install lftp:amd64 <none> 4.9.2-2ubuntu1.1"),
        "Package installation evidence differs")
    zone = ZoneInfo("America/New_York")
    times = (datetime.fromisoformat(jobs["23985"][4]).replace(tzinfo=zone),
             datetime.fromisoformat(jobs["23986"][4]).replace(tzinfo=zone),
             datetime.strptime(line[:19], "%Y-%m-%d %H:%M:%S").replace(tzinfo=zone))
    result = classification["result"]
    require([t.isoformat() for t in times] == [result[k] for k in
        ("native_end", "prior_inventory_latest_at", "installed_at")],
        "Chronology differs from bound evidence")
    return times


def replay_runtime(plan, run, session, verification, prior, evidence, additions, times):
    lookup = read(plan["runtime_lookup"])
    binding = read(lookup["binding"])
    specs = binding["runtime_specs"]
    historical_checker(specs, prior)
    historical = runtime_review(plan, run, session, verification, evidence,
        tree_checker=lambda current_specs: historical_checker(current_specs, prior))
    require(historical["phases"] == prior["phases"], "Historical bracket replay differs")
    comparisons = []
    started = time.monotonic()
    for number, (path, digest) in enumerate(specs):
        manifest_ref = dict(record(path), sha256=digest)
        manifest = read(manifest_ref)
        evidence.bind(path)
        comparison = compare_inventory(manifest, inventory(manifest["roots"]))
        if number == 0:
            classify(comparison, additions, *times)
        else:
            require(comparison["equal"] is True, "Private runtime no longer matches exactly")
        comparisons.append(dict(manifest=manifest_ref, comparison=comparison))
    # The hook replays retained descriptors. Its temporary output is NOT fresh evidence.
    return dict(schema=SCHEMA, status="historical_brackets_and_current_exact_additions_revalidated",
        historical_phases=historical["phases"],
        historical_first_terminal_inventory=historical["fresh_terminal_inventory"],
        historical_first_terminal_rechecked=True, current_inventories=comparisons,
        current_inventory_check_wall_s=time.monotonic() - started,
        current_original_inventory_equality=False, current_original_entries_unchanged=True,
        current_private_inventory_equality=True, current_additions=additions,
        historical_check_time_is_upper_bound=True,
        native_end=times[0].isoformat(), prior_inventory_latest_at=times[1].isoformat(),
        installed_at=times[2].isoformat(), continuous_runtime_integrity_established=False,
        original_ordinary_full_review_success=False, full_review_admitted=False,
        terminal_reviewed=False, next_identity_authorized=False, accuracy_evaluated=False,
        resource_measurements_admitted=False, automatic_retry=False, publication_ready=False,
        limitations=["This component does not replace full resource/environment/output review.",
            "Current equality to the original OS inventory is false, with exactly the bound late additions.",
            "Retained snapshots and package chronology do not establish continuous integrity or installer identity.",
            "Original review failures remain failures; no segmentation-fault cause or fix is established."])


def review(destination, expected_source_sha256):
    destination = Path(destination)
    require(destination.is_absolute() and destination.resolve() == destination
        and destination == DEFAULT_DESTINATION and not destination.exists(),
        "Require one fresh fixed component destination")
    source = record(__file__)
    require(source["sha256"] == expected_source_sha256, "Prospective runtime reviewer changed")
    request_ref, classification_ref = record(REQUEST), record(CLASSIFICATION)
    require(request_ref["sha256"] == REQUEST_SHA
        and classification_ref["sha256"] == CLASSIFICATION_SHA, "Prospective bindings differ")
    evidence = Evidence()
    classification = read(classification_ref)
    require(classification.get("schema") == "native11_postterminal_runtime_addition_classification_v1"
        and classification.get("classification_only") is True
        and classification.get("private_runtime_comparison_equal") is True
        and all(classification.get(k) is False for k in ("current_full_runtime_revalidated",
            "original_ordinary_full_review_success", "full_review_admitted", "next_identity_authorized")),
        "Require the original nonadmitting classification")
    for ref in [classification_ref, *classification["checked_records"],
                classification["prior_inventory_source"]]:
        check(ref)
        evidence.bind(ref["path"])
    failure = read(classification["failed_review_observation"])
    times = chronology(classification, failure)
    request = read(request_ref)
    execution, plan = amendment(request["amendment"])
    validate_request(request, request["amendment"], execution, 23985)
    run = plan["runs"][request["index"]]
    require(run["index"] == 11 and run["cell"] == "p1_c1_r0", "Wrong native identity")
    terminal = verify_terminal(23985)
    fields = scheduler_fields(terminal)
    require(fields.get("JobState", fields.get("State")) == "COMPLETED"
        and fields.get("ExitCode") == "0:0", "Native outcome is not successful")
    if terminal["source"] == "live_controller":
        require(fields.get("Comment") == request_ref["sha256"], "Native comment differs")
    for path, digest in SOURCES.items():
        ref = record(ROOT / path)
        require(ref["sha256"] == digest, "Frozen runtime kernel changed")
        evidence.bind(ref["path"])
    for ref in [source, request_ref, request["plan"], request["amendment"], plan["baseline"],
                *execution["new_sources"], *plan["helper_sources"], *plan["evidence"]]:
        check(ref)
        evidence.bind(ref["path"])
    session = Path(plan["panel_root"]) / "sessions" / "run_11"
    result = evidence.json(session / "result.json")
    verification = evidence.json(Path(run["output_root"]) / "verification.json")
    bind_session(result, verification, request_ref, request, request["amendment"], run)
    prior_ref = classification["original_runtime"]
    prior = read(prior_ref)
    evidence.bind(prior_ref["path"])
    destination.mkdir(parents=True, exist_ok=False)
    save(destination / "scheduler.json", terminal)
    scheduler_ref = record(destination / "scheduler.json")
    try:
        report = replay_runtime(plan, run, session, verification, prior, evidence,
            classification["result"]["additional_entries"], times)
        report.update(source=source, request=request_ref, classification=classification_ref,
            prior_runtime=prior_ref, plan=request["plan"], amendment=request["amendment"],
            scheduler=scheduler_ref, index=11, job_id=23985, cell=run["cell"], observed_unix_ns=time.time_ns())
        report["evidence"] = evidence.finish()
        check(source)
        check(scheduler_ref)
        save(destination / "runtime.json", report)
        return record(destination / "runtime.json")
    except Exception as error:
        save(destination / "failure.json", dict(schema=SCHEMA, status="runtime_component_failed",
            error_type=type(error).__name__, error=str(error), source=source,
            request=request_ref, classification=classification_ref, full_review_admitted=False,
            next_identity_authorized=False, automatic_retry=False, publication_ready=False))
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-sha256", required=True)
    args = parser.parse_args()
    print(json.dumps(review(DEFAULT_DESTINATION, args.source_sha256), sort_keys=True))
