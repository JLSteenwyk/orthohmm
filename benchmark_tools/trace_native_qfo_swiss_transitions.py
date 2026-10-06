"""Trace actual SwissTrees pair decisions in the admitted native P0/C0 cells."""

import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import audit_native_qfo_swiss_counts as ordinary
from benchmark_tools import audit_recovered_native_qfo_swiss_counts as recovered
from benchmark_tools import export_native_qfo_scientific_scores as reporter
from benchmark_tools.audit_qfo_swiss_counts import HEADER, LABELS, statistics
from benchmark_tools.bootstrap_qfo_swiss_stages import METRICS, aggregate
from benchmark_tools.export_native_factorial_progress import load, require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

import numpy as np

CELLS = ("p0_c0_r0", "p0_c0_r1")
AUDITS = (
    ("native_qfo_swiss_family_count_audit_v1", "supplied_native_swiss_family_counts_verified", ordinary),
    ("recovered_native_qfo_swiss_family_count_audit_v1",
     "supplied_recovered_native_swiss_family_counts_verified", recovered),
)
TRANSITIONS = tuple(a + "->" + b for a in LABELS for b in LABELS)


def read_labels(path, families):
    require(len(families) == 18 and len(set(families)) == 18, "Require18 distinct reference families")
    labels = {}
    with gzip.open(path, "rt", newline="") as handle:
        require(handle.readline().rstrip("\r\n") == HEADER, "Changed native raw header")
        for row in csv.reader(handle, delimiter="\t"):
            require(len(row) == 4, "Invalid native raw width")
            family, a, b, label = row
            require(family in families and a and b and a != b and label in LABELS,
                    "Invalid family/pair/label")
            key = (family, *sorted((a, b)))
            require(key not in labels, "Duplicate canonical family pair")
            labels[key] = label
    require({key[0] for key in labels} == set(families), "Missing reference family")
    return labels


def verified_counts(labels, families, audited):
    require([row["family"] for row in audited["families"]] == families, "Changed audited family order")
    counts = {f: Counter({label: 0 for label in LABELS}) for f in families}
    members = {f: set() for f in families}
    for (family, a, b), label in labels.items():
        counts[family][label] += 1
        members[family].update((a, b))
    values, seen = [], set()
    for row in audited["families"]:
        family = row["family"]
        require(dict(counts[family]) == row["counts_without_prior"]
                and sorted(members[family]) == row["represented_genes"]
                and len(members[family]) > 5 and not seen.intersection(members[family]),
                "Raw counts/members differ from disjoint audited universe")
        seen.update(members[family])
        scores = statistics(counts[family])
        require(scores == row["statistics_with_prior"], "Changed native family arithmetic")
        values.append([scores["PPV"], scores["TPR"]])
    point = dict(zip(METRICS, aggregate(np.asarray(values).mean(axis=0)).tolist()))
    require(all(abs(point[m] - audited["aggregate"][m]) <= 1e-12 for m in METRICS),
            "Raw counts differ from audited macro statistic")
    return counts


def compare_labels(before, after, families, before_audit, after_audit):
    require(before.keys() == after.keys(), "Changed reference relation universe")
    left = verified_counts(before, families, before_audit)
    right = verified_counts(after, families, after_audit)
    tables = {f: Counter({t: 0 for t in TRANSITIONS}) for f in families}
    changes = []
    for key in sorted(before):
        a, b = before[key], after[key]
        require((a in ("TP", "FN")) == (b in ("TP", "FN")), "Changed reference truth label")
        tables[key[0]][a + "->" + b] += 1
        if a != b:
            changes.append((*key, a, b))
    total = Counter({t: 0 for t in TRANSITIONS})
    rows = []
    for f in families:
        table = tables[f]
        require(all(sum(table[a + "->" + b] for b in LABELS) == left[f][a] for a in LABELS)
                and all(sum(table[a + "->" + b] for a in LABELS) == right[f][b] for b in LABELS),
                "Transition marginals differ from raw counts")
        total.update(table)
        rows.append(dict(family=f, relations=sum(table.values()), transitions=dict(table),
            before_counts=dict(left[f]), after_counts=dict(right[f]),
            removed_true_positives=table["TP->FN"], removed_false_positives=table["FP->TN"],
            added_true_positives=table["FN->TP"], added_false_positives=table["TN->FP"]))
    return dict(families=rows, reference_relations=len(before), transitions=dict(total),
        changed_relations=len(changes), removed_true_positives=total["TP->FN"],
        removed_false_positives=total["FP->TN"], added_true_positives=total["FN->TP"],
        added_false_positives=total["TN->FP"],
        after_predictions_subset_on_reference=not (total["FN->TP"] or total["TN->FP"])), changes


def run(snapshot_path, snapshot_sha, audits, output, changes_path):
    require(output != changes_path and all(not p.exists() and not p.is_symlink()
            for p in (output, changes_path)), "Require two distinct fresh output paths")
    evidence = []
    snapshot, snapshot_ref = load(snapshot_path, snapshot_sha, evidence)
    require(snapshot["schema"] == "native_qfo_scientific_reporting_snapshot_v1"
            and snapshot["source"] == record(reporter.__file__) and snapshot["publication_ready"] is False,
            "Wrong scientific snapshot/source")
    normal, repairs = [], []
    for row in snapshot["rows"]:
        destination = (normal if row["status"] == "supplied_native_admission" else repairs
                       if row["status"] == "supplied_recovered_scientific_admission" else None)
        if destination is not None:
            destination.append((row["admission"]["path"], row["admission"]["sha256"]))
    replay = reporter.collect(snapshot["plan"]["path"], snapshot["plan"]["sha256"], normal, repairs)
    require(all(snapshot[k] == value for k, value in replay.items()), "Scientific snapshot replay differs")
    rows = {row["cell"]: row for row in snapshot["rows"]}
    require(len(rows) == len(snapshot["rows"]) and all(rows[c]["accuracy_admitted"] is True for c in CELLS),
            "Require two distinct admitted native cells")
    require(rows[CELLS[1]]["resources"] is None and rows[CELLS[1]]["timing_admitted"] is False
            and rows[CELLS[1]]["timing_eligible"] is False, "Recovered science relabels failed timing")
    require(len(audits) == 2, "Require ordinary and recovered count audits")
    labels, count_rows, bindings, families, reference = [], [], [], None, None
    for cell, (path, digest), (schema, status, source) in zip(CELLS, audits, AUDITS):
        audit, ref = load(path, digest, evidence)
        require(audit["schema"] == schema and audit["status"] == status
                and audit["source"] == record(source.__file__) and audit["publication_ready"] is False
                and audit["new_accuracy_or_resource_admission"] is False
                and audit["historical_intervals_attached"] is False
                and audit["independent_confirmation"] is False
                and type(audit["new_bootstrap_draws"]) is int and audit["new_bootstrap_draws"] == 0,
                "Changed count audit source/scope")
        candidates = [r for r in audit["cells"] if r["cell"] == cell]
        require(len(candidates) == 1, "Missing or duplicate audited cell")
        counts, native = candidates[0], rows[cell]
        require(all(counts[k] == native[k] for k in ("index", "native_job_id", "admission"))
                and counts["native_endpoint_f1"] == native["scores"]["SwissTrees"],
                "Count audit/native admission mismatch")
        if families is None:
            families, reference = audit["families"], audit["reference"]
        require(audit["families"] == families and audit["reference"] == reference,
                "Changed reference families or source")
        admission, _ = load(native["admission"]["path"], native["admission"]["sha256"], evidence)
        execution_ref = admission["execution_report"]
        require(execution_ref in admission["checked_records"], "Unadmitted execution evidence")
        execution, observed = load(execution_ref["path"], execution_ref["sha256"], evidence)
        require(observed == execution_ref, "Execution metadata differs")
        raw = counts["raw_file"]
        require(raw in admission["checked_records"] and raw in execution["outputs"]
                and raw in audit["checked_inputs"] and Path(raw["path"]).parent.name == "SwissTrees"
                and raw["path"].endswith("raw.txt.gz"), "Unadmitted SwissTrees raw evidence")
        check(raw)
        evidence.extend([raw, reference, audit["source"], *audit["helpers"]])
        values = read_labels(Path(raw["path"]), families)
        require(len(values) == audit["reference_relation_count"], "Incomplete audited reference universe")
        labels.append(values)
        count_rows.append(counts)
        bindings.append(dict(cell=cell, admission=native["admission"], count_audit=ref, raw=raw,
            native_endpoint_f1=native["scores"]["SwissTrees"], macro_statistics=counts["aggregate"],
            **{k: native[k] for k in ("timing_eligible", "timing_admitted", "resources") if k in native}))
    comparison, changes = compare_labels(*labels, families, *count_rows)
    evidence.extend([snapshot["source"], record(__file__), record(ordinary.__file__),
        record(recovered.__file__), record(statistics.__code__.co_filename), record(aggregate.__code__.co_filename)])
    for ref in evidence:
        check(ref)
    with changes_path.open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("family", "protein_a", "protein_b", "before", "after"))
        writer.writerows(changes)
    result = dict(schema="native_qfo_swiss_pair_transitions_v1", snapshot=snapshot_ref,
        source=record(__file__), cells=bindings, reference=reference, comparison=comparison,
        changed_relations_ledger=record(changes_path), checked_inputs=evidence,
        new_scoring_or_admission=False, new_bootstrap_draws=0, uncertainty_admitted=False,
        scientific_timings_admitted=False, independent_confirmation=False, publication_ready=False,
        limitations=[
            "Only the18 SwissTrees families and native P0/C0 cells; not all predicted relations or selected defaults.",
            "Transitions describe actual assessed pairs, not causal proof of reconciliation events or tree correctness.",
            "Initial HMM search remains on; R also changes clique to resolved-pair prediction semantics.",
            "Pooled unsmoothed pair counts are descriptive, not the native prior-adjusted macro F1 endpoint.",
            "All changed relations are retained in a deterministic ledger, not selected examples.",
            "Family dependence and development exposure persist; no uncertainty or generalization admission.",
            "Recovered inference timing remains failed/null/ineligible; contention effects unknown/tool-dependent.",
            "Direct admitted evidence is checked; no expensive transitive admission/resource replay or new job."])
    for ref in evidence:
        check(ref)
    with output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("snapshot", "output", "changes"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--snapshot-sha256", required=True)
    parser.add_argument("--counts-audit", nargs=2, action="append", required=True)
    args = parser.parse_args()
    result = run(args.snapshot, args.snapshot_sha256, args.counts_audit, args.output, args.changes)
    print(json.dumps({k: v for k, v in result["comparison"].items() if k != "families"}, sort_keys=True))
