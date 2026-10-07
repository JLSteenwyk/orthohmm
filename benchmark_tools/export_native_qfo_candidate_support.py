"""Describe serialized accepted-event features without rescoring or tuning."""

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import statistics
import sys

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "benchmark_tools/results"
REPORT_SHA = "3b3ab8f1132e70362e1450ac990caf7a9dad536e528bc54d6af39db16183cd53"
READBACK_SHA = "4afc4aa6853a69027b9dc710602309f44923315b1e68dbae4cd1e61866f4b92d"
TRACE_SHA = "ec867165230af4b80e8df7d0227ccd5719910a41dcc83211022a7c265240dc48"
LEDGER_SHA = "3f0598f7e00cfc0942146604ab1863ad9f03f2d5017efd5244c7e76a29d43642"
HEADER = ["protein_left", "protein_right", "native_gene_left", "native_gene_right",
          "baseline_state", "candidate_state", "baseline_group_left", "baseline_group_right",
          "candidate_group", "first_connected_iteration", "first_direct_event", "connection_path"]
FEATURES = ["source_size", "target_size", "source_seed_families", "target_seed_families",
            "support", "margin", "species_overlap_fraction", "forward_hits", "reverse_hits",
            "forward_average", "reverse_average", "forward_maximum", "reverse_maximum",
            "forward_coverage", "reverse_coverage", "forward_normalized_support", "reverse_normalized_support"]
COHORTS = ["TP-only", "FP-only", "mixed", "no changed scored VGNC pair"]
PAIR_HEADER = HEADER + ["event_cohort"] + FEATURES
EVENT_HEADER = ["event_index", "iteration", "source_cluster", "target_cluster", "cohort",
                "recovered_tp_pairs", "added_fp_pairs"] + FEATURES
PROTOCOL = RESULTS / "NATIVE_QFO_CANDIDATE_SUPPORT_PROTOCOL_20261006.md"


def need(condition, message):
    if not condition:
        raise ValueError(message)


def identity(path):
    path = Path(path).resolve(strict=True)
    need(path.is_file(), "Not a regular file")
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest.hexdigest())


def verify(ref):
    need(identity(ref["path"]) == ref, "Changed current bytes: " + ref["path"])


def read_json(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            need(key not in result, "Duplicate JSON key")
            result[key] = value
        return result
    with Path(path).open() as stream:
        return json.load(stream, object_pairs_hook=unique)


def scalar(value):
    return "positive_infinity" if value == math.inf else str(value)


def validate_trace(trace):
    need(isinstance(trace, list) and trace, "Empty or malformed accepted trace")
    seen = set()
    for event in trace:
        need(isinstance(event, dict) and type(event.get("iteration")) is int
             and event["iteration"] in (0, 1), "Invalid event round")
        memberships = []
        for role in ("source", "target"):
            members = event.get(role + "_genes")
            need(isinstance(members, list) and members and all(type(g) is str and g
                 and g.strip() == g and not any(c in g for c in "\t\n\r") for g in members)
                 and len(set(members)) == len(members), "Invalid event members")
            need(type(event.get(role + "_size")) is int and event[role + "_size"] == len(members),
                 "Declared event size differs")
            need(type(event.get(role + "_cluster")) is int and event[role + "_cluster"] >= 0,
                 "Invalid event cluster")
            memberships.append(tuple(sorted(members)))
        need(not set(memberships[0]).intersection(memberships[1]), "Overlapping event endpoints")
        key = (event["iteration"], *memberships)
        need(key not in seen, "Duplicate semantic accepted event")
        seen.add(key)
        for feature in FEATURES:
            value = event.get(feature)
            need(type(value) in (int, float) and value >= 0
                 and (math.isfinite(value) or feature == "margin" and value == math.inf),
                 "Invalid event feature: " + feature)
            if feature.endswith("_hits") or feature.endswith("_seed_families"):
                need(type(value) is int, "Invalid event count: " + feature)


def summaries(events):
    result = {}
    for feature in FEATURES:
        values = [event[feature] for event in events if math.isfinite(event[feature])]
        result[feature] = dict(finite=len(values), missing=0,
            positive_infinity=sum(event[feature] == math.inf for event in events),
            minimum=min(values) if values else None, median=statistics.median(values) if values else None,
            mean=math.fsum(values) / len(values) if values else None, maximum=max(values) if values else None)
    return result


def analyze(trace, rows):
    validate_trace(trace)
    need(isinstance(rows, list) and rows, "Empty changed-pair ledger")
    counts = [[0, 0] for _ in trace]
    pairs, native_pairs, localized = set(), set(), Counter()
    for row in rows:
        need(len(row) == len(HEADER) and all(type(s) is str for s in row), "Malformed pair row")
        a, b, ga, gb, old, new, gl, gr, group, round_text, event_text, path = row
        need(a and a < b and (a, b) not in pairs and ga and gb and ga != gb
             and tuple(sorted((ga, gb))) not in native_pairs, "Duplicate or invalid pair identity")
        pairs.add((a, b)); native_pairs.add(tuple(sorted((ga, gb))))
        need((old, new) in (("FN", "TP"), ("not_scored", "FP")) and gl and gr and gl != gr and group,
             "Wrong pair transition or group identity")
        need(round_text in ("0", "1"), "Wrong pair connection round")
        need(path in ("direct_cross_endpoint", "transitive_union"), "Wrong connection path")
        localized[int(round_text), new, path] += 1
        if path == "transitive_union":
            need(event_text == "", "Transitive pair invents a direct event")
            continue
        need(event_text.isdecimal() and str(int(event_text)) == event_text
             and int(event_text) < len(trace), "Invalid direct event reference")
        event = trace[int(event_text)]
        need(event["iteration"] == int(round_text), "Direct event round differs")
        left, right = set(event["source_genes"]), set(event["target_genes"])
        need(ga in left and gb in right or gb in left and ga in right, "Direct event endpoints differ")
        counts[int(event_text)][0 if new == "TP" else 1] += 1
    cohorts = [("mixed" if tp and fp else "TP-only" if tp else "FP-only" if fp
                else "no changed scored VGNC pair") for tp, fp in counts]
    annotated = []
    for row in rows:
        if row[-1] == "transitive_union":
            annotated.append(row + ["transitive_no_direct_event"] + [""] * len(FEATURES))
        else:
            index = int(row[-2])
            annotated.append(row + [cohorts[index]] + [scalar(trace[index][f]) for f in FEATURES])
    implicated = [[str(index), str(event["iteration"]), str(event["source_cluster"]), str(event["target_cluster"]),
                   cohorts[index], str(counts[index][0]), str(counts[index][1])] + [scalar(event[f]) for f in FEATURES]
                  for index, event in enumerate(trace) if sum(counts[index])]
    cohort_summary = []
    for round_value in (None, 0, 1):
        for cohort in COHORTS:
            selected = [event for index, event in enumerate(trace) if cohorts[index] == cohort
                        and (round_value is None or event["iteration"] == round_value)]
            cohort_summary.append(dict(iteration=round_value, cohort=cohort, events=len(selected), features=summaries(selected)))
    totals = dict(accepted_events=len(trace), changed_pairs=len(rows), implicated_events=len(implicated),
                  direct_tp_pairs=sum(tp for tp, _ in counts), direct_fp_pairs=sum(fp for _, fp in counts),
                  transitive_tp_pairs=sum(r[-1] == "transitive_union" and r[5] == "TP" for r in rows),
                  transitive_fp_pairs=sum(r[-1] == "transitive_union" and r[5] == "FP" for r in rows))
    localized_summary = [dict(iteration=r, candidate_state=s, connection_path=p, pairs=n)
                         for (r, s, p), n in sorted(localized.items(), key=lambda item: (item[0][1], item[0][0], item[0][2]))]
    return annotated, implicated, dict(totals=totals, cohort_summaries=cohort_summary, localized_summary=localized_summary)


def inputs(report_path, readback_path, report_sha=REPORT_SHA, readback_sha=READBACK_SHA,
           trace_sha=TRACE_SHA, ledger_sha=LEDGER_SHA):
    report_ref, readback_ref = identity(report_path), identity(readback_path)
    need(report_ref["sha256"] == report_sha and readback_ref["sha256"] == readback_sha, "Prior evidence anchor differs")
    prior, readback = read_json(report_path), read_json(readback_path)
    need(prior.get("schema") == "native_qfo_candidate_alias_group_trace_v1"
         and prior.get("status") == "original_protein_identity_join_and_all_changed_pairs_localized"
         and readback.get("schema") == "native_qfo_candidate_alias_group_readback_v1"
         and readback.get("status") == "original_protein_identity_join_and_all_paths_independently_verified"
         and readback.get("report") == report_ref, "Prior evidence scope differs")
    for key in ("accepted_merges", "changed_pairs", "genes", "baseline_groups", "candidate_groups", "bridges", "localized_summary"):
        need(prior[key] == readback[key], "Prior independent agreement differs: " + key)
    for evidence in (prior, readback):
        need(all(evidence.get(k) is True for k in ("whole_candidate_partition_reconstructed", "complete_pair_localization",
             "original_scored_identifiers_preserved", "failed_r1_timing_remains_ineligible"))
             and all(evidence.get(k) is False for k in ("accuracy_rescored", "native_inference_reexecuted",
             "new_scoring_or_admission", "uncertainty_admitted", "publication_ready")), "Prior evidence flags differ")
    refs = [report_ref, readback_ref, prior["native_trace"], prior["ledger"], identity(__file__),
            identity(Path(__file__).with_name("readback_native_qfo_candidate_support.py")), identity(PROTOCOL)]
    need(prior["native_trace"]["sha256"] == trace_sha and prior["ledger"]["sha256"] == ledger_sha
         and prior["native_trace"] in prior["checked_records"],
         "Trace/ledger anchor differs")
    for ref in refs:
        verify(ref)
    trace = read_json(prior["native_trace"]["path"])
    with Path(prior["ledger"]["path"]).open(newline="") as stream:
        table = csv.reader(stream, delimiter="\t")
        need(next(table, None) == HEADER, "Wrong complete ledger header")
        rows = list(table)
    return prior, refs, trace, rows


def write_table(path, header, rows):
    with path.open("x", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(header); writer.writerows(rows)


def export(report_path, readback_path, output):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    try:
        prior, refs, trace, rows = inputs(report_path, readback_path)
        annotated, implicated, analysis = analyze(trace, rows)
        need(len(trace) == prior["accepted_merges"] and len(rows) == prior["changed_pairs"]
             and analysis["localized_summary"] == prior["localized_summary"], "Complete ledger counts differ")
        write_table(output / "annotated_pairs.tsv", PAIR_HEADER, annotated)
        write_table(output / "implicated_events.tsv", EVENT_HEADER, implicated)
        for ref in refs:
            verify(ref)
        report = dict(schema="native_qfo_candidate_accepted_support_v1", status="descriptive_accepted_event_support_exported",
            source=identity(__file__), reader_source=refs[5], protocol=refs[6], alias_report=refs[0], alias_readback=refs[1],
            native_trace=refs[2], ledger=refs[3], checked_records=refs,
            annotated_pairs=identity(output / "annotated_pairs.tsv"), implicated_events=identity(output / "implicated_events.tsv"),
            runtime=dict(python=sys.version, platform=platform.platform()), **analysis,
            unit="each accepted event once; dependent pair counts only descriptive",
            missing_feature_policy="all required features present; missing count zero",
            infinite_margin_policy="positive_infinity kept separate from finite summaries",
            grouping_replayed=False, aliases_rederived=False, accuracy_rescored=False, native_inference_reexecuted=False,
            new_scoring_or_admission=False, uncertainty_admitted=False, defaults_changed=False,
            calibrated_confidence=False, causal_mechanism_established=False, publication_ready=False,
            failed_r1_timing_remains_ineligible=True,
            limitations=["Development-exposed association conditional on accepted events, not rejected alternatives.",
                         "No changed VGNC pair is unlabeled, not a true negative or correct event.",
                         "A group-level direct event is not evidence of a direct pairwise HMM hit.",
                         "Membership, sizes and reference scope confound descriptive differences.",
                         "Shared-host postprocessing is not native inference timing or isolated-speed evidence."])
        with (output / "report.json").open("x") as stream:
            json.dump(report, stream, sort_keys=True, indent=2, allow_nan=False); stream.write("\n")
        return report
    except Exception as exc:
        with (output / "failure.json").open("x") as stream:
            json.dump(dict(schema="native_qfo_candidate_accepted_support_failure_v1", error=str(exc),
                           automatic_retry=False, publication_ready=False), stream, indent=2); stream.write("\n")
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alias-report", type=Path, default=RESULTS / "native_qfo_candidate_alias_group_20261006_v1/report.json")
    parser.add_argument("--alias-readback", type=Path, default=RESULTS / "native_qfo_candidate_alias_group_readback_20261006_v1.json")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = export(args.alias_report, args.alias_readback, args.output)
    print(json.dumps(dict(status=result["status"], **result["totals"]), sort_keys=True))
