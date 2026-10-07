"""Independently check accepted-event joins and rational feature summaries."""

import argparse
import csv
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LEDGER_COLUMNS = ["protein_left", "protein_right", "native_gene_left", "native_gene_right", "baseline_state", "candidate_state",
                  "baseline_group_left", "baseline_group_right", "candidate_group", "first_connected_iteration",
                  "first_direct_event", "connection_path"]
METRICS = ["source_size", "target_size", "source_seed_families", "target_seed_families", "support", "margin",
           "species_overlap_fraction", "forward_hits", "reverse_hits", "forward_average", "reverse_average",
           "forward_maximum", "reverse_maximum", "forward_coverage", "reverse_coverage", "forward_normalized_support",
           "reverse_normalized_support"]
LABELS = ("TP-only", "FP-only", "mixed", "no changed scored VGNC pair")


def check(condition, message):
    if not condition:
        raise ValueError(message)


def fingerprint(path):
    path = Path(path).resolve(strict=True)
    check(path.is_file(), "Nonfile evidence")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(1048576)
            if not block:
                break
            digest.update(block)
    return {"path": str(path), "sha256": digest.hexdigest(), "bytes": path.stat().st_size}


def unchanged(reference):
    check(fingerprint(reference["path"]) == reference, "Evidence changed: " + reference["path"])


def json_data(path):
    def object_from_pairs(pairs):
        check(len(pairs) == len({k for k, _ in pairs}), "Duplicate JSON member")
        return dict(pairs)
    with Path(path).open() as handle:
        return json.load(handle, object_pairs_hook=object_from_pairs)


def table_data(path, columns):
    with Path(path).open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        check(next(reader, None) == columns, "Incorrect table columns")
        rows = list(reader)
    check(all(len(row) == len(columns) for row in rows), "Partial table row")
    return rows


def validated_events(trace):
    check(type(trace) is list and len(trace) > 0, "No accepted events")
    result, semantic_keys = {}, set()
    for index, event in enumerate(trace):
        check(type(event) is dict and type(event.get("iteration")) is int and event["iteration"] in {0, 1}, "Bad event iteration")
        ends = []
        for side in ("source", "target"):
            genes = event.get(side + "_genes")
            check(type(genes) is list and len(genes) > 0 and all(type(g) is str and g.strip() == g and g
                  and "\t" not in g and "\n" not in g and "\r" not in g for g in genes)
                  and len(set(genes)) == len(genes), "Bad event membership")
            check(type(event.get(side + "_size")) is int and event[side + "_size"] == len(genes), "Bad declared size")
            check(type(event.get(side + "_cluster")) is int and event[side + "_cluster"] >= 0, "Bad cluster name")
            ends.append(frozenset(genes))
        check(ends[0].isdisjoint(ends[1]), "Endpoints overlap")
        semantic = (event["iteration"], *ends)
        check(semantic not in semantic_keys, "Repeated event identity")
        semantic_keys.add(semantic)
        for metric in METRICS:
            value = event.get(metric)
            check(type(value) in {int, float} and value >= 0 and (math.isfinite(value)
                  or metric == "margin" and value == float("inf")), "Invalid event numeric field: " + metric)
            if metric in {"forward_hits", "reverse_hits", "source_seed_families", "target_seed_families"}:
                check(type(value) is int, "Noninteger event count")
        result[index] = (event, *ends)
    return result


def rational_features(events):
    summary = {}
    for metric in METRICS:
        finite = sorted(Fraction(str(event[metric])) for event in events if math.isfinite(event[metric]))
        n = len(finite)
        median = (finite[(n - 1) // 2] + finite[n // 2]) / 2 if n else None
        summary[metric] = dict(finite=n, missing=0, positive_infinity=len(events) - n,
            minimum=float(finite[0]) if n else None, median=float(median) if n else None,
            mean=float(sum(finite, Fraction(0)) / n) if n else None, maximum=float(finite[-1]) if n else None)
    return summary


def reconstruct(trace, pair_rows):
    events = validated_events(trace)
    check(type(pair_rows) is list and pair_rows, "No changed pairs")
    joined = {index: {"TP": [], "FP": []} for index in events}
    identities, native_identities, path_counts = set(), set(), {}
    for row in pair_rows:
        check(len(row) == len(LEDGER_COLUMNS) and all(type(v) is str for v in row), "Invalid ledger row")
        entry = dict(zip(LEDGER_COLUMNS, row))
        scored = (entry["protein_left"], entry["protein_right"])
        genes = frozenset((entry["native_gene_left"], entry["native_gene_right"]))
        check(scored[0] and scored[0] < scored[1] and scored not in identities and len(genes) == 2
              and "" not in genes and genes not in native_identities, "Duplicate or invalid pair")
        identities.add(scored); native_identities.add(genes)
        state = entry["candidate_state"]
        check((entry["baseline_state"], state) in {("FN", "TP"), ("not_scored", "FP")}
              and entry["baseline_group_left"] and entry["baseline_group_right"] and entry["candidate_group"]
              and entry["baseline_group_left"] != entry["baseline_group_right"], "Invalid scored change")
        check(entry["first_connected_iteration"] in {"0", "1"}, "Invalid connection round")
        round_value, path = int(entry["first_connected_iteration"]), entry["connection_path"]
        check(path in {"transitive_union", "direct_cross_endpoint"}, "Unknown path")
        key = (state, round_value, path)
        path_counts[key] = path_counts.get(key, 0) + 1
        event_name = entry["first_direct_event"]
        if path == "transitive_union":
            check(event_name == "", "Transitive path has direct event")
        else:
            check(event_name.isdecimal() and str(int(event_name)) == event_name and int(event_name) in events,
                  "Bad event reference")
            event, source, target = events[int(event_name)]
            check(event["iteration"] == round_value and len(genes & source) == 1 and len(genes & target) == 1,
                  "Event membership or iteration mismatch")
            joined[int(event_name)][state].append(entry)
    labels = {}
    for index, outcome in joined.items():
        labels[index] = ("mixed" if outcome["TP"] and outcome["FP"] else "TP-only" if outcome["TP"]
                         else "FP-only" if outcome["FP"] else "no changed scored VGNC pair")
    def metric_strings(event):
        return ["positive_infinity" if event[m] == float("inf") else str(event[m]) for m in METRICS]
    annotated = []
    for row in pair_rows:
        if row[11] == "transitive_union":
            annotated.append(row + ["transitive_no_direct_event"] + [""] * len(METRICS))
        else:
            index = int(row[10])
            annotated.append(row + [labels[index]] + metric_strings(events[index][0]))
    event_rows = []
    for index in events:
        event = events[index][0]
        if joined[index]["TP"] or joined[index]["FP"]:
            event_rows.append([str(index), str(event["iteration"]), str(event["source_cluster"]), str(event["target_cluster"]),
                labels[index], str(len(joined[index]["TP"])), str(len(joined[index]["FP"]))] + metric_strings(event))
    summaries = []
    for iteration in (None, 0, 1):
        for label in LABELS:
            selected = [event for index, (event, _, _) in events.items() if labels[index] == label
                        and (iteration is None or event["iteration"] == iteration)]
            summaries.append(dict(iteration=iteration, cohort=label, events=len(selected), features=rational_features(selected)))
    total = dict(accepted_events=len(trace), changed_pairs=len(pair_rows), implicated_events=len(event_rows),
        direct_tp_pairs=sum(len(value["TP"]) for value in joined.values()), direct_fp_pairs=sum(len(value["FP"]) for value in joined.values()),
        transitive_tp_pairs=sum(n for (state, _, path), n in path_counts.items() if state == "TP" and path == "transitive_union"),
        transitive_fp_pairs=sum(n for (state, _, path), n in path_counts.items() if state == "FP" and path == "transitive_union"))
    localized = [dict(candidate_state=state, iteration=iteration, connection_path=path, pairs=n)
                 for (state, iteration, path), n in sorted(path_counts.items())]
    return annotated, event_rows, dict(totals=total, cohort_summaries=summaries, localized_summary=localized)


def equivalent(actual, expected):
    if type(expected) is dict:
        check(type(actual) is dict and actual.keys() == expected.keys(), "Summary dictionary differs")
        for key in expected:
            equivalent(actual[key], expected[key])
    elif type(expected) is list:
        check(type(actual) is list and len(actual) == len(expected), "Summary list differs")
        for left, right in zip(actual, expected):
            equivalent(left, right)
    elif type(expected) is float:
        check(type(actual) in {int, float} and math.isfinite(actual)
              and math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-12), "Rational numerical summary differs")
    else:
        check(type(actual) is type(expected) and actual == expected, "Exact count or identity differs")


def review(report_path, report_sha):
    report_ref = fingerprint(report_path)
    check(report_ref["sha256"] == report_sha, "Exported report anchor differs")
    report = json_data(report_path)
    check(report.get("schema") == "native_qfo_candidate_accepted_support_v1"
          and report.get("status") == "descriptive_accepted_event_support_exported", "Wrong report schema")
    for name in ("grouping_replayed", "aliases_rederived", "accuracy_rescored", "native_inference_reexecuted",
                 "new_scoring_or_admission", "uncertainty_admitted", "defaults_changed", "calibrated_confidence",
                 "causal_mechanism_established", "publication_ready"):
        check(report.get(name) is False, "Wrong scope: " + name)
    check(report.get("failed_r1_timing_remains_ineligible") is True, "Failed timing scope differs")
    refs = [report[k] for k in ("alias_report", "alias_readback", "native_trace", "ledger", "source", "reader_source", "protocol")]
    check(refs == report["checked_records"] and len({r["path"] for r in refs}) == 7, "Checked evidence inventory differs")
    anchors = ["3b3ab8f1132e70362e1450ac990caf7a9dad536e528bc54d6af39db16183cd53",
               "4afc4aa6853a69027b9dc710602309f44923315b1e68dbae4cd1e61866f4b92d",
               "ec867165230af4b80e8df7d0227ccd5719910a41dcc83211022a7c265240dc48",
               "3f0598f7e00cfc0942146604ab1863ad9f03f2d5017efd5244c7e76a29d43642"]
    check([r["sha256"] for r in refs[:4]] == anchors, "Original selected evidence differs")
    check(report["reader_source"] == fingerprint(__file__)
          and report["source"] == fingerprint(Path(__file__).with_name("export_native_qfo_candidate_support.py"))
          and report["protocol"] == fingerprint(ROOT / "benchmark_tools/results/NATIVE_QFO_CANDIDATE_SUPPORT_PROTOCOL_20261006.md"),
          "Source or protocol binding differs")
    for ref in refs + [report_ref, report["annotated_pairs"], report["implicated_events"]]:
        unchanged(ref)
    prior, old_reader = json_data(refs[0]["path"]), json_data(refs[1]["path"])
    check(prior.get("schema") == "native_qfo_candidate_alias_group_trace_v1"
          and prior.get("status") == "original_protein_identity_join_and_all_changed_pairs_localized"
          and old_reader.get("schema") == "native_qfo_candidate_alias_group_readback_v1"
          and old_reader.get("status") == "original_protein_identity_join_and_all_paths_independently_verified"
          and old_reader["report"] == refs[0] and prior["native_trace"] == refs[2] and prior["ledger"] == refs[3]
          and refs[2] in prior["checked_records"], "Original report/readback mismatch")
    for field in ("genes", "baseline_groups", "candidate_groups", "accepted_merges", "changed_pairs", "bridges", "localized_summary"):
        check(prior[field] == old_reader[field], "Original independent evidence differs")
    for previous in (prior, old_reader):
        check(all(previous.get(k) is True for k in ("whole_candidate_partition_reconstructed", "complete_pair_localization",
              "original_scored_identifiers_preserved", "failed_r1_timing_remains_ineligible"))
              and all(previous.get(k) is False for k in ("accuracy_rescored", "native_inference_reexecuted",
              "new_scoring_or_admission", "uncertainty_admitted", "publication_ready")), "Original scope mismatch")
    trace = json_data(refs[2]["path"])
    ledger = table_data(refs[3]["path"], LEDGER_COLUMNS)
    pair_rows, event_rows, analysis = reconstruct(trace, ledger)
    check(len(trace) == prior["accepted_merges"] and len(ledger) == prior["changed_pairs"]
          and analysis["localized_summary"] == prior["localized_summary"], "Original complete counts differ")
    check(table_data(report["annotated_pairs"]["path"], LEDGER_COLUMNS + ["event_cohort"] + METRICS) == pair_rows,
          "Complete annotated pair table differs")
    check(table_data(report["implicated_events"]["path"], ["event_index", "iteration", "source_cluster", "target_cluster",
          "cohort", "recovered_tp_pairs", "added_fp_pairs"] + METRICS) == event_rows, "Unique implicated event table differs")
    for key, value in analysis.items():
        equivalent(report[key], value)
    for ref in refs + [report_ref, report["annotated_pairs"], report["implicated_events"]]:
        unchanged(ref)
    return dict(schema="native_qfo_candidate_accepted_support_readback_v1", status="accepted_event_features_independently_verified",
        report=report_ref, source=fingerprint(__file__), totals=analysis["totals"], localized_summary=analysis["localized_summary"],
        cohort_counts=[{k: row[k] for k in ("iteration", "cohort", "events")} for row in analysis["cohort_summaries"]],
        exporter_or_union_kernel_imported=False, rational_summaries_verified=True, numeric_absolute_tolerance=1e-12,
        numeric_relative_tolerance=1e-12, original_scored_identifiers_preserved=True, uncertainty_admitted=False,
        defaults_changed=False, calibrated_confidence=False, causal_mechanism_established=False, publication_ready=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    check(not args.output.exists(), "Existing readback output must not be overwritten")
    result = review(args.report, args.report_sha256)
    with args.output.open("x") as handle:
        json.dump(result, handle, sort_keys=True, indent=2, allow_nan=False); handle.write("\n")
    print(json.dumps(dict(status=result["status"], **result["totals"]), sort_keys=True))
