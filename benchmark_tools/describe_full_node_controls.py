"""Prespecified descriptive control comparisons; no independent-interval inference."""

import argparse
import gzip
import json
from pathlib import Path

from benchmark_tools.describe_dual_native_flags import describe, summarize_rows
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_full_node_controls import ORDER
from benchmark_tools.validate_full_node_control import common_intervals

AUDIT_SHA = "bd8e11ead6ee2e979784c8cd97730d3cf82ee62e4343c9cdf18c0475d2221cd3"


def inventory_check(audit):
    expected = [(i, block, mode) for i, (block, mode) in enumerate(
        (block, mode) for block, modes in enumerate(ORDER) for mode in modes)]
    actual = [(row["index"], row["block"], row["mode"]) for row in audit["trials"]]
    if actual != expected or audit["validated_trials"] != 9 or any(
            row["status"] != "validated" for row in audit["trials"]):
        raise ValueError("Require all nine validated controls in frozen order")


def witness(row):
    files = [item for item in row["replay"]["evidence"] if Path(item["path"]).name == "workload_done.json"]
    if len(files) != 1:
        raise ValueError("Missing or ambiguous workload completion witness")
    check(files[0])
    done = json.loads(Path(files[0]["path"]).read_text())
    workers = done["workers"]
    result = dict(worker_count=len(workers), creations=sum(w["creations"] for w in workers),
        self_cpu_s=sum(w["self_cpu_s"] for w in workers),
        waited_child_cpu_s=sum(w["waited_child_user_s"]+w["waited_child_system_s"] for w in workers))
    common = row["replay"]["trial"]["validation"]
    result["common_work_s"] = (common["common_finished_ns"]-common["common_started_ns"])/1e9
    # Whole-work totals and common-window counters have different boundaries.
    result["worker_creation_rates_per_s"] = [w["creations"] / ((w["finished_ns"]-w["started_ns"])/1e9)
                                            for w in workers]
    check(files[0])
    return result, files[0]


def block_differences(rows):
    blocks = []
    for block in range(3):
        selected = {row["mode"]: row for row in rows if row["block"] == block}
        if set(selected) != {"steady", "churn", "contended"} or sum(row["block"] == block for row in rows) != 3:
            raise ValueError("Incomplete or duplicated block")
        baseline = selected["steady"]["common"]["distributions"]
        contrasts = {}
        for mode in ("churn", "contended"):
            values = selected[mode]["common"]["distributions"]
            if set(values) != set(baseline):
                raise ValueError("Inconsistent block measurements")
            contrasts[mode+"_minus_steady"] = {key: values[key]["median"]-baseline[key]["median"] for key in values}
        blocks.append(dict(block=block, common_interval_median_differences=contrasts))
    return blocks


def report(path):
    source = record(path)
    if source["sha256"] != AUDIT_SHA:
        raise ValueError("Wrong pinned full-node audit")
    audit = json.loads(gzip.decompress(path.read_bytes()))
    inventory_check(audit)
    rows, evidence = [], [source]
    for row in audit["trials"]:
        measurement = row["replay"]["measurement"]
        measured = measurement["measured"]
        adapted = dict(index=row["index"], method=row["mode"], job_id=measured["job_id"],
            started_ns=measured["native"]["started_ns"], screening=measurement["screening"],
            narrow_flagged_intervals=measurement["narrow_flagged_intervals"])
        description = describe(adapted, measured, [])
        common = row["replay"]["trial"]["common_intervals"]
        if common != common_intervals(measured["points"], row["replay"]["trial"]["validation"]):
            raise ValueError("Common-work subset does not reproduce")
        selected = [item for item in description["intervals"] if item["index"] in common]
        if len(selected) != len(common) or not selected:
            raise ValueError("Invalid common-work interval inventory")
        work, item = witness(row)
        evidence.append(item)
        rows.append(dict(index=row["index"], block=row["block"], mode=row["mode"],
            all=description["all"], common=summarize_rows(selected),
            intervals=description["intervals"], common_interval_indices=common, workload=work,
            common_flags=sum(item["flagged"] for item in selected),
            positive_control_detected=row["positive_control_detected"]))
    for item in evidence:
        check(item)
    return dict(status="full_node_controls_described", trials=rows, blocks=block_differences(rows),
        evidence=evidence, source=record(__file__), scientific_timings_admitted=False, publication_ready=False,
        limitations=["Prespecified descriptive comparisons; no interval-level significance or confidence intervals.",
        "Host, frontier, pressure and process witnesses span distinct windows; do not subtract them causally.",
        "Whole-worker creation counts are not counts confined to the common observation subset.",
        "Neither flags nor root/direct counts identify foreign activity; native-tool historical flags remain unresolved.",
        "No threshold changes, selective repeats, overhead correction or scientific timing admission."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = report(args.audit.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
