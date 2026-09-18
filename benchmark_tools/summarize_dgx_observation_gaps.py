"""Describe retained observation gaps without reclassifying timing eligibility."""

import argparse
from collections import Counter
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def resource_errors(rows, summary):
    types, events, affected, observations, peak_affected = Counter(), 0, 0, 0, 0
    for row in rows:
        observations += 1
        errors = row["snapshot"]["sampling_errors"]
        affected += bool(errors)
        events += len(errors)
        types.update(e["error_type"] for e in errors)
        if (errors and row["snapshot"]["sampled_sum_process_rss_bytes"]
                == summary["maximum_sampled_sum_rss_bytes"]):
            peak_affected += 1
    if observations != summary["observations"] or affected != summary["samples_with_process_errors"]:
        raise ValueError("Resource error counts differ from reproduced summary")
    return {"observations": observations, "error_samples": affected, "error_events": events,
        "error_types": dict(types), "fraction_of_samples_with_errors": affected / observations if observations else None,
        "error_samples_at_recorded_rss_maximum": peak_affected}


def assemble(results):
    paths = [results / name for name in ("dgx_completed_panel_review_20260918.json", "dgx_resource_replay_20260918.json")]
    sources = [record(p) for p in paths]
    host, resource = [json.loads(p.read_text()) for p in paths]
    if (host["review_failures"] != 0 or resource["failures"] != 0
            or len(host["runs"]) != 27 or len(resource["runs"]) != 27):
        raise ValueError("Require complete successful host and resource replays")
    rows, raw_sources = [], []
    for index, (h, r) in enumerate(zip(host["runs"], resource["runs"])):
        if (h["index"] != index or r["index"] != index
                or h["host_review"]["job_id"] != r["job_id"]
                or r["status"] != "resource_accounting_reproduced_not_timing_admission"):
            raise ValueError("Host/resource run identities differ")
        matches = [item for item in r["evidence"] if Path(item["path"]).name == "samples.jsonl"]
        if len(matches) != 1:
            raise ValueError("Missing or ambiguous resource evidence")
        raw = matches[0]
        check(raw)
        with Path(raw["path"]).open() as stream:
            errors = resource_errors((json.loads(line) for line in stream), r["summary"])
        check(raw)
        raw_sources.append(raw)
        names = Counter()
        for event in h["host_review"]["diagnostics"]["unmatched_events"]:
            if not event["name"].startswith("kworker/"):
                names[event["name"]] += event["count"]
        summary = h["host_review"]["retained_host_summary"]
        rows.append({"index": index, "method": h["method"], "proteomes": h["proteomes"], "repeat": h["repeat"],
            "resource_errors": errors, "host_status": summary["status"],
            "observed_persistent_foreign_core_maximum": summary["maximum_observed_foreign_average_cores"],
            "fixed_monitor_core_threshold": summary["threshold_average_cores"],
            "unmatched_non_kworker_names": dict(names),
            "unmatched_identity_events": h["host_review"]["diagnostics"]["unmatched_identity_events"],
            "kworker_named_events": h["host_review"]["diagnostics"]["kworker_named_identity_events"],
            "host_snapshot_errors": h["host_review"]["diagnostics"]["snapshot_error_types"],
            "host_observation_errors": summary["observation_errors"],
            "command_bracketed": summary["command_bracketed_by_samples"],
            "scientific_timing_admitted": False})
    error_types = Counter()
    for row in rows:
        error_types.update(row["resource_errors"]["error_types"])
    for source in sources:
        check(source)
    return {"status": "observation_gaps_described_no_reclassification", "source": record(__file__),
        "inputs": sources, "raw_sources": raw_sources, "runs": rows,
        "resource_error_event_types": dict(error_types), "scientific_timings_admitted": 0,
        "limitations": ["NoSuchProcess records a process-read failure, not its unobserved RSS or CPU consumption.",
            "No error at the sampled maximum does not rule out missing between-sample peaks.",
            "Process names neither authenticate kernel threads nor establish benign user-space work.",
            "Persistent foreign CPU below threshold does not bound unmatched or entirely unsampled work.",
            "No original status, threshold, value, or timing inclusion rule has been changed."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = assemble(args.results.resolve())
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps(report["resource_error_event_types"], sort_keys=True))
