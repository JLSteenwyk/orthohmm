"""Independently check initial-edge trace arithmetic against admitted pair rows."""

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path


def record(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest.hexdigest())


def indexed(rows):
    result = {}
    for row in rows:
        key = tuple(row[k] for k in ("refog", "left", "right"))
        if key in result or key[1] >= key[2]:
            raise ValueError("Duplicate or unordered pair")
        result[key] = row
    return result


def audit_rows(original, traced):
    original, traced = indexed(original), indexed(traced)
    if original.keys() != traced.keys():
        raise ValueError("Pair universe differs")
    families, summary = {}, Counter()
    thresholds = {}
    for key, row in traced.items():
        source = original[key]
        scores = []
        for new, old in (("forward_score", "forward_normalized_hit"),
                         ("reverse_score", "reverse_normalized_hit")):
            value = None if row[new] == "NA" else float(row[new])
            expected = None if source[old] == "NA" else float(source[old])
            if value != expected or (value is not None and (not math.isfinite(value) or value <= 0)):
                raise ValueError("Hit score differs or is invalid")
            if value is not None:
                scores.append(value)
        ends = []
        for gene, column in ((key[1], "left_threshold"), (key[2], "right_threshold")):
            value = float(row[column])
            if math.isnan(value) or value <= 0:
                raise ValueError("Invalid threshold")
            if gene in thresholds and thresholds[gene] != value:
                raise ValueError("Inconsistent gene threshold")
            thresholds[gene] = value
            ends.append(value)
        if not scores:
            decision = "no_direct_hit"
        elif any(score >= endpoint for score in scores for endpoint in ends):
            decision = "initial_edge"
        elif all(math.isinf(endpoint) for endpoint in ends):
            decision = "no_finite_endpoint_threshold"
        else:
            decision = "below_endpoint_threshold"
        if row["decision"] != decision:
            raise ValueError("Wrong threshold classification")
        for stage in ("multipass_refined", "root_hogs"):
            if row[stage] not in ("True", "False") or row[stage] != source[stage]:
                raise ValueError("Grouping flag differs or is invalid")
        label = decision + "/root_" + row["root_hogs"]
        summary[label] += 1
        families.setdefault(key[0], Counter())[label] += 1
    return dict(pair_memberships=len(traced), families=families, summary=summary,
                threshold_genes=len(thresholds))


def run(report_path, output):
    if output.exists():
        raise FileExistsError(output)
    report_record = record(report_path)
    report = json.loads(report_path.read_text())
    if report["status"] != "frozen_initial_rbnh_reference_pairs_traced" or report["job_id"] != "21797":
        raise ValueError("Unexpected trace execution")
    records = [report_record, report["source"], report["table"], *report["checked_records"]]
    for item in records:
        if record(item["path"]) != item:
            raise ValueError("Evidence file changed")
    sources = [item for item in report["checked_records"] if Path(item["path"]).name == "reference_pair_trace.tsv"]
    if len(sources) != 1:
        raise ValueError("Missing unique source pair table")
    with Path(sources[0]["path"]).open(newline="") as original, Path(report["table"]["path"]).open(newline="") as traced:
        result = audit_rows(csv.DictReader(original, delimiter="\t"), csv.DictReader(traced, delimiter="\t"))
    for field in ("pair_memberships", "families", "summary"):
        if result[field] != report[field]:
            raise ValueError("Reported counts differ: " + field)
    for item in records:
        if record(item["path"]) != item:
            raise ValueError("Evidence changed during audit")
    result.update(status="initial_edge_pair_arithmetic_verified", checked_records=records,
                  auditor=record(__file__),
                  limitations=["Checks trace arithmetic and source identity, not independent native graph reconstruction.",
                               "Scheduler terminal state is recorded separately in the execution report.",
                               "Raw family pair memberships are not official weighted recall or causal effects."])
    with output.open("x") as stream:
        json.dump(result, stream, sort_keys=True, indent=2)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.report, args.output)
