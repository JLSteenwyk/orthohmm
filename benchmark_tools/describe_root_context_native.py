"""Describe native context intervals and flag subsets without timing admission."""

import argparse
import gzip
import json
from pathlib import Path

from benchmark_tools.describe_root_context_controls import distributions
from benchmark_tools.measure_native_root_context import evaluate
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.verify_root_context_native_provenance import same

METHODS = ["orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full"]


def summary(rows):
    return dict(intervals=len(rows), distributions=distributions([r["values"] for r in rows]) if rows else None,
        original_flags=sum(r["original_flagged"] for r in rows),
        narrow_flags=sum(r["narrow_flagged"] for r in rows),
        root_membership_changes=sum(r["root_membership_changed"] for r in rows))


def intervals(row, measured):
    points = measured["points"]
    context = evaluate(points, row["job_id"])
    if not same(context, row["root_context"]):
        raise ValueError("Raw context differs from audited context")
    count = len(points)-1
    for key in ("original_flagged_intervals", "narrow_flagged_intervals"):
        values = row[key]
        if values != sorted(set(values)) or any(type(i) is not int or i < 0 or i >= count for i in values):
            raise ValueError("Invalid flag inventory")
    result = []
    for index, comparison in enumerate(context["intervals"]):
        left, right = points[index]["root_context"], points[index+1]["root_context"]
        values = {"scope_cpu_usec:" + k: v for k, v in comparison["scope_cpu_usec"].items()}
        values.update({k: comparison[k] for k in ("root_minus_system_cpu_usec", "root_minus_three_named_children_cpu_usec")})
        values.update({"host_ticks:" + k: v for k, v in comparison["enclosing_host_category_ticks"].items()})
        values["host_enclosing_window_s"] = (right["host_after"]["finished_ns"]-left["host_before"]["started_ns"])/1e9
        for a, b in zip(left["rows"], right["rows"]):
            if a["scope"] != b["scope"]:
                raise ValueError("Scope order differs")
            values["scope_enclosing_window_s:" + a["scope"]] = (b["finished_ns"]-a["started_ns"])/1e9
            values["left_scope_read_s:" + a["scope"]] = (a["finished_ns"]-a["started_ns"])/1e9
            values["right_scope_read_s:" + b["scope"]] = (b["finished_ns"]-b["started_ns"])/1e9
        result.append(dict(index=index, values=values, ticks=comparison["ticks"],
            original_flagged=index in row["original_flagged_intervals"],
            narrow_flagged=index in row["narrow_flagged_intervals"],
            root_membership_changed=comparison["observed_root_membership_changed"]))
    return result


def describe(audit):
    if (audit["status"] != "root_context_native_audited_not_scientific_admission"
            or audit["scientific_timings_admitted"] is not False
            or [(r["index"], r["method"]) for r in audit["runs"]] != list(enumerate(METHODS))):
        raise ValueError("Wrong native audit or task inventory")
    rows = []
    for row in audit["runs"]:
        result = {k: row[k] for k in ("index", "method", "status")}
        if row["status"] not in {"validated", "output_mismatch"}:
            rows.append(dict(result, subsets=None, retained=row))
            continue
        suffix = f"/root_context_native_v1/run_{row['index']:02d}/measurement/lineage_report.json"
        records = [r for r in audit["inventory"] if r["path"].endswith(suffix)]
        if len(records) != 1:
            raise ValueError("Missing or duplicate audited lineage file")
        source = records[0]
        check(source)
        measured = json.loads(Path(source["path"]).read_text())
        values = intervals(row, measured)
        check(source)
        subsets = dict(all=summary(values))
        for kind in ("original", "narrow"):
            for flagged in (True, False):
                name = kind + ("_flagged" if flagged else "_unflagged")
                subsets[name] = summary([r for r in values if r[kind + "_flagged"] is flagged])
        result.update(subsets=subsets, intervals=values, lineage=source,
            native_wall_s=row["native_wall_s"], output_equivalent=row["output_equivalent"])
        rows.append(result)
    return dict(status="native_root_context_described", runs=rows, panel_issues=audit["issues"],
        scientific_timings_admitted=False, publication_ready=False,
        limitations=["One diagnostic per method; interval distributions are not independent replicates or confidence intervals.",
            "Flagged and unflagged subsets retain all original flags; no exclusion or timing correction.",
            "Signed residuals are not causal task attribution; named children are not an exhaustive partition.",
            "Host guest categories overlap user/nice and are reported separately, not summed.",
            "Scope and host enclosing windows differ; empty subsets remain missing, not zero-valued distributions.",
            "Output mismatches and panel issues remain visible; no overhead, isolation or speed ranking is established."])


def report(path, expected_sha):
    source = record(path)
    if source["sha256"] != expected_sha:
        raise ValueError("Wrong pinned native audit")
    raw = path.read_bytes()
    result = describe(json.loads(gzip.decompress(raw) if path.suffix == ".gz" else raw))
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
