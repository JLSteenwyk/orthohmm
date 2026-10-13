"""Independently read back fresh search arrays and compare diagnostic links."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.run_search_recovery import link_changes, resolve_settings
from benchmark_tools.search_decision_trace import classify_search_result
from benchmark_tools.trace_ob_search_decisions import load_watched


def read_raw(path):
    with np.load(path, allow_pickle=False) as raw:
        if set(raw.files) != {"query_ids", "target_ids", "query_indices", "target_indices",
                              "scores", "evalues", "candidate_count"}:
            raise ValueError("Unknown search snapshot fields")
        values = {name: raw[name].copy() for name in raw.files}
    for name in ("query_ids", "target_ids", "query_indices", "target_indices", "scores", "evalues"):
        if values[name].ndim != 1:
            raise ValueError("Require one-dimensional search arrays")
    count = values["candidate_count"]
    if count.shape != () or count.dtype.kind not in "iu" or count.item() < 0:
        raise ValueError("Invalid native candidate count")
    result = SimpleNamespace(**{name: values[name] for name in
                                ("query_indices", "target_indices", "scores", "evalues")},
                             candidate_count=int(count.item()))
    return result, values["query_ids"].tolist(), values["target_ids"].tolist()


def audit_direction(item, query_ids, target_ids, pairs, watched, thresholds):
    check(item["raw"])
    result, queries, targets = read_raw(item["raw"]["path"])
    if queries != query_ids or targets != target_ids:
        raise ValueError("Snapshot sequence universe/order differs from FASTA")
    if (item["candidate_count"] != result.candidate_count or item["query_count"] != len(queries)
            or item["target_count"] != len(targets) or item["watched_pairs"] != len(pairs)):
        raise ValueError("Reported search dimensions differ from arrays")
    if [row["threshold"] for row in item["thresholds"]] != thresholds:
        raise ValueError("Wrong significance thresholds")
    output = {}
    for entry in item["thresholds"]:
        threshold = entry["threshold"]
        rows = classify_search_result(result, queries, targets, pairs, threshold)
        counts, recovered = Counter(), 0
        expected_table = []
        for row in rows:
            historical = watched[row["query"], row["target"]]["historical_present"]
            counts[row["decision"]] += 1
            recovered += row["decision"] == "accepted" and not historical
            expected_table.append({key: "" if value is None else str(value)
                                   for key, value in {**row, "historical_present": historical}.items()})
        check(entry["table"])
        with Path(entry["table"]["path"]).open(newline="") as stream:
            observed_table = list(csv.DictReader(stream, delimiter="\t"))
        if observed_table != expected_table:
            raise ValueError("Watched table differs from independent raw-array classification")
        significant = int(np.count_nonzero(result.evalues < threshold))
        if (dict(counts) != entry["decisions"] or recovered != entry["historically_absent_accepted"]
                or significant != entry["all_significant_hits"]):
            raise ValueError("Reported search counts differ from raw arrays")
        output[threshold] = {"rows": rows, "decisions": dict(counts),
                            "historically_absent_accepted": recovered,
                            "all_significant_hits": significant}
    return output


def compare_historical_control(historical, current):
    if historical is None:
        return {"status": "unavailable", "scope": "Historical comparison is supporting evidence, not a prerequisite."}
    check(historical)
    previous = json.loads(Path(historical["path"]).read_text())
    if previous.get("forced_watched_candidates"):
        raise ValueError("Cannot compare ordinary search to forced candidates")
    old = {(item["query_species"], item["target_species"]): item for item in previous["directions"]}
    new = {(item["query_species"], item["target_species"]): item for item in current["directions"]}
    if len(old) != len(previous["directions"]) or set(old) != set(new):
        raise ValueError("Historical/current direction universes differ")
    comparisons = []
    for key in sorted(new):
        check(old[key]["raw"])
        check(new[key]["raw"])
        a, qa, ta = read_raw(old[key]["raw"]["path"])
        b, qb, tb = read_raw(new[key]["raw"]["path"])
        comparisons.append({"query_species": key[0], "target_species": key[1],
            "identities_equal": qa == qb and ta == tb,
            **{field + "_equal": bool(np.array_equal(getattr(a, field), getattr(b, field)))
               for field in ("query_indices", "target_indices", "scores", "evalues")}})
    exact = all(all(value for key, value in row.items() if key.endswith("_equal"))
                for row in comparisons)
    return {"status": "exact_persisted_cpu_array_agreement" if exact else "historical_cpu_array_differences",
            "directions": comparisons,
            "scope": "Contemporary CPU observation compared to persisted audited CPU arrays, not authentication of the original publication execution."}


def summarize(execution_path, output):
    execution = json.loads(execution_path.read_text())
    for item in [*execution["inputs"], *execution["sources"], *execution["runtime"]["files"]]:
        check(item)
    watched = load_watched(execution["pair_trace"]["path"])
    species, owners = {}, {}
    for item in execution["fastas"]:
        name = Path(item["path"]).name
        species[name] = [sequence.id for sequence in SeqIO.parse(item["path"], "fasta")]
        for gene in species[name]:
            if gene in owners:
                raise ValueError("Duplicate FASTA identity")
            owners[gene] = name
    directions = {}
    for pair in watched:
        directions.setdefault((owners[pair[0]], owners[pair[1]]), []).append(pair)
    expected = sorted(directions)
    report = {"status": "partial_independent_readback", "execution_plan": record(execution_path),
              "audit_source": record(__file__), "arms": [], "accuracy_evaluated": False,
              "defaults_promoted": False,
              "scope": "Reference-conditioned search links, not full-pipeline accuracy or comparative efficiency."}
    audited = {}
    for arm in execution["plan"]["arms"]:
        label = arm["label"]
        path = Path(execution["output"]) / label / "report.json"
        outcome = {"label": label}
        if not path.exists():
            outcome["status"] = "report_missing_execution_state_unknown"
            report["arms"].append(outcome)
            continue
        raw_report = json.loads(path.read_text())
        outcome["retained_report"] = record(path)
        if raw_report["status"] == "failed":
            outcome.update(status="failed", error=raw_report.get("error"))
            report["arms"].append(outcome)
            continue
        if (raw_report["execution_plan"] != report["execution_plan"]
                or raw_report["source_revision"] != execution["source_revision"]
                or raw_report["label"] != label
                or raw_report["settings"] != resolve_settings(execution["plan"], label)
                or raw_report["status"] != "fresh_search_completed_pending_independent_audit"
                or raw_report["explicit_control_fixture_exact"] is not True):
            raise ValueError("Search identity/status differs from frozen execution")
        if [(row["query_species"], row["target_species"]) for row in raw_report["directions"]] != expected:
            raise ValueError("Missing, duplicated, or reordered search directions")
        observed = {}
        for key, item in zip(expected, raw_report["directions"]):
            pairs = sorted(directions[key])
            wanted = {q for q, _ in pairs}
            queries = [gene for gene in species[key[0]] if gene in wanted]
            observed[key] = audit_direction(item, queries, species[key[1]], pairs, watched,
                                           execution["plan"]["significance_thresholds"])
        audited[label] = observed
        totals = []
        for threshold in execution["plan"]["significance_thresholds"]:
            counts = Counter()
            for entries in observed.values():
                counts.update(entries[threshold]["decisions"])
            totals.append({"threshold": threshold, "decisions": dict(counts),
                "historically_absent_accepted": sum(entries[threshold]["historically_absent_accepted"]
                                                     for entries in observed.values()),
                "all_significant_hits": sum(entries[threshold]["all_significant_hits"]
                                            for entries in observed.values())})
        outcome.update(status="independently_verified", thresholds=totals,
            candidates=sum(item["candidate_count"] for item in raw_report["directions"]),
            resources=raw_report["resources"])
        if label == "control":
            outcome["historical_cpu_comparison"] = compare_historical_control(
                execution.get("historical_control"), raw_report)
        report["arms"].append(outcome)
    if "control" in audited:
        for outcome in report["arms"]:
            if outcome["status"] != "independently_verified":
                continue
            label = outcome["label"]
            for total in outcome["thresholds"]:
                threshold = total["threshold"]
                changes, baseline_changes = Counter(), Counter()
                for key in expected:
                    current = audited[label][key][threshold]["rows"]
                    changes.update(link_changes(audited["control"][key][threshold]["rows"], current))
                    baseline_changes.update(link_changes(audited["control"][key][0.0001]["rows"], current))
                total["versus_control_same_threshold"] = dict(changes)
                total["versus_control_evalue_1e4"] = dict(baseline_changes)
    if len(audited) == len(execution["plan"]["arms"]):
        report["status"] = "all_search_arms_independently_verified"
    write_json(output, report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execution-plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    summarize(args.execution_plan.resolve(), args.output.resolve())
