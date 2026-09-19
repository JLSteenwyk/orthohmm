"""Join admitted search decisions to retained final-group pair observations."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.trace_ob_initial_edges import TRACE_SHA

AUDIT_SHA = "8081c8ff005c911b0122ba3a797aaee51c1ac65a1c3f17f551512d3cb47ce730"
DECISIONS = {"accepted", "not_selected_by_prefilter", "scored_not_significant"}


def category(left, right):
    if left not in DECISIONS or right not in DECISIONS:
        raise ValueError("Unknown search decision")
    if "accepted" in (left, right):
        return "at_least_one_accepted_direction"
    if left == right:
        return "both_prefilter_excluded" if left == "not_selected_by_prefilter" else "both_scored_not_significant"
    return "one_prefilter_excluded_one_scored_not_significant"


def join(rows, decisions):
    counts, families, joined, seen = Counter(), {}, [], set()
    for row in rows:
        identity = (row["refog"], row["left"], row["right"])
        if identity in seen or row["left"] == row["right"]:
            raise ValueError("Duplicate family pair or self pair")
        seen.add(identity)
        forward = decisions[row["left"], row["right"]]
        reverse = decisions[row["right"], row["left"]]
        if ((forward == "accepted") != (row["forward_normalized_hit"] != "NA")
                or (reverse == "accepted") != (row["reverse_normalized_hit"] != "NA")):
            raise ValueError("Historical hit presence disagrees with observed diagnostic")
        if row["root_hogs"] not in ("True", "False") or row["same_species"] not in ("True", "False"):
            raise ValueError("Invalid grouping/species flag")
        label = category(forward, reverse)
        state = "grouped" if row["root_hogs"] == "True" else "separated"
        species = "same_species" if row["same_species"] == "True" else "cross_species"
        key = ":".join((label, state, species))
        counts[key] += 1
        families.setdefault(row["refog"], Counter())[key] += 1
        joined.append(dict(refog=row["refog"], left=row["left"], right=row["right"],
                           forward=forward, reverse=reverse, search_category=label,
                           final_grouping=state, species_relation=species))
    return joined, dict(counts), {f: dict(c) for f, c in sorted(families.items())}


def run(root, output):
    if output.exists():
        raise FileExistsError(output)
    audit_path = root / "benchmark_tools/results/ob_search_decisions_audit_20260918.json"
    audit = read_frozen(audit_path, AUDIT_SHA)
    check(audit["source_report"])
    report = json.loads(Path(audit["source_report"]["path"]).read_text())
    trace_path = root / "benchmark_tools/results/ob_family_trace_verified_20260916.json"
    trace = read_frozen(trace_path, TRACE_SHA)
    checked = [record(audit_path), audit["source_report"], record(trace_path), trace["pair_trace"]]
    decisions = {}
    for direction in report["directions"]:
        item = direction["table"]
        if item not in audit["checked_records"]:
            raise ValueError("Unadmitted direction table")
        check(item)
        checked.append(item)
        with Path(item["path"]).open(newline="") as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                pair = row["query"], row["target"]
                if pair in decisions:
                    raise ValueError("Duplicate directed observation")
                decisions[pair] = row["decision"]
    for item in checked:
        check(item)
    if len(decisions) != audit["directed_pairs"]:
        raise ValueError("Incomplete directed inventory")
    with Path(trace["pair_trace"]["path"]).open(newline="") as stream:
        rows, counts, families = join(csv.DictReader(stream, delimiter="\t"), decisions)
    if len(rows) != 40733 or len(families) != 70 or sum(counts.values()) != len(rows):
        raise ValueError("Wrong family-pair inventory")
    output.mkdir(parents=True)
    table = output / "pairs.tsv"
    with table.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for item in checked:
        check(item)
    summary = dict(status="admitted_search_grouping_observations_joined", source=record(__file__),
                   checked_records=checked, pair_memberships=len(rows), families=families,
                   counts=counts, table=record(table), accuracy_evaluated=False,
                   publication_ready=False, limitations=[
                       "Descriptive reference pair memberships, not official weighted recall or independent units.",
                       "Observed search presence matches historical trace; numerical/runtime equivalence not established.",
                       "Co-membership can arise through indirect paths and later stages; no causal F1 attribution.",
                       "Includes within-species and low-certainty pairs; no counterfactual rescoring or tuning."])
    with (output / "summary.json").open("x") as stream:
        json.dump(summary, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.root.resolve(), args.output.resolve())
