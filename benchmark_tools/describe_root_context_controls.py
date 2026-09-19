"""Describe fixed-block root-context controls without independent-interval inference."""

import argparse
import gzip
import json
import math
from pathlib import Path
from statistics import median

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.root_context_control_design import ORDER


def distributions(rows):
    if not rows:
        raise ValueError("No observations for distribution")
    keys = set(rows[0])
    if any(set(row) != keys for row in rows):
        raise ValueError("Distribution fields differ")
    result = {}
    for key in sorted(keys):
        values = [row[key] for row in rows]
        if any(type(value) not in (int, float) or not math.isfinite(value) for value in values):
            raise ValueError("Nonfinite or nonnumeric observation")
        result[key] = dict(n=len(values), minimum=min(values), median=median(values), maximum=max(values))
    return result


def subset(intervals):
    return dict(intervals=len(intervals), distributions=distributions([row["values"] for row in intervals]),
        original_flags=sum(row["original_flagged"] for row in intervals),
        narrow_flags=sum(row["narrow_flagged"] for row in intervals),
        root_membership_changes=sum(row["root_membership_changed"] for row in intervals))


def describe(audit):
    expected = [(i, b, m) for i, (b, m) in enumerate((b, m) for b, modes in enumerate(ORDER) for m in modes)]
    actual = [(row["index"], row["retained"]["block"], row["retained"]["mode"]) for row in audit["trials"]]
    if (audit["status"] != "root_context_controls_audited" or actual != expected
            or audit["scientific_timings_admitted"] is not False
            or audit["validated_trials"] != sum(row["status"] == "validated" for row in audit["trials"])):
        raise ValueError("Audit identity or fixed inventory differs")
    rows = []
    for row in audit["trials"]:
        result = dict(index=row["index"], block=row["retained"]["block"], mode=row["retained"]["mode"], status=row["status"])
        if row["status"] != "validated":
            result.update(all=None, common=None, retained=row["retained"])
            rows.append(result)
            continue
        replay = row["replay"]
        native = replay["measurement"]["lineage"]
        points = native["measured"]["points"]
        comparisons = replay["measurement"]["context"]["intervals"]
        if len(comparisons) != len(points)-1:
            raise ValueError("Observation interval inventory differs")
        common = replay["trial"]["common_intervals"]
        if (not common or common != sorted(set(common))
                or any(type(i) is not int or i < 0 or i >= len(comparisons) for i in common)):
            raise ValueError("Invalid common-work subset")
        intervals = []
        for index, context in enumerate(comparisons):
            left, right = points[index]["root_context"], points[index+1]["root_context"]
            values = {"scope_cpu_usec:" + key: value for key, value in context["scope_cpu_usec"].items()}
            values.update({key: context[key] for key in ("root_minus_system_cpu_usec", "root_minus_three_named_children_cpu_usec")})
            values.update({"host_ticks:" + key: value for key, value in context["enclosing_host_category_ticks"].items()})
            values["host_enclosing_window_s"] = (right["host_after"]["finished_ns"] - left["host_before"]["started_ns"])/1e9
            for a, b in zip(left["rows"], right["rows"]):
                if a["scope"] != b["scope"]:
                    raise ValueError("Scope order differs")
                values["scope_enclosing_window_s:" + a["scope"]] = (b["finished_ns"]-a["started_ns"])/1e9
            intervals.append(dict(index=index, values=values, original_flagged=index in native["original_flagged_intervals"],
                narrow_flagged=index in native["narrow_flagged_intervals"],
                root_membership_changed=context["observed_root_membership_changed"], ticks=context["ticks"]))
        result.update(all=subset(intervals), common=subset([intervals[i] for i in common]), intervals=intervals,
            common_interval_indices=common, enclosing_user_cpu_s=replay["trial"]["enclosing_user_cpu_s"],
            positive_control_response=replay["trial"]["positive_control_response"])
        rows.append(result)
    blocks = []
    for block in range(3):
        selected = {row["mode"]: row for row in rows if row["block"] == block}
        contrasts = {}
        for mode in ("idle", "churn", "user-contended"):
            baseline, candidate = selected["steady"]["common"], selected[mode]["common"]
            if baseline is None or candidate is None:
                contrasts[mode + "_minus_steady"] = None
                continue
            a, b = baseline["distributions"], candidate["distributions"]
            if set(a) != set(b):
                raise ValueError("Block distribution fields differ")
            contrasts[mode + "_minus_steady"] = {key: b[key]["median"]-a[key]["median"] for key in a}
        blocks.append(dict(block=block, common_interval_median_differences=contrasts))
    return dict(status="root_context_controls_described", trials=rows, blocks=blocks,
        scientific_timings_admitted=False, publication_ready=False,
        limitations=["Fixed-block descriptive differences, not confidence intervals or independent-interval tests.",
            "Scope and host counters have distinct non-atomic enclosing windows; no causal subtraction or timing correction.",
            "Guest host ticks overlap user/nice ticks and are reported separately, never summed into a total.",
            "Signed residuals are preserved; root membership changes do not identify task CPU or unobserved transients.",
            "Failures and unrun conditions remain missing observations; no selective replacement or native timing admission."])


def report(path, expected_sha):
    source = record(path)
    if source["sha256"] != expected_sha:
        raise ValueError("Wrong pinned root-control audit")
    raw = path.read_bytes()
    audit = json.loads(gzip.decompress(raw) if path.suffix == ".gz" else raw)
    result = describe(audit)
    check(source)
    return dict(result, audit=source, source=record(__file__))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--audit-sha", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = report(args.audit.resolve(), args.audit_sha)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
