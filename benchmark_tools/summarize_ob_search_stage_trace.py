"""Joint direct-hit and grouping outcomes over all admitted OrthoBench trace pairs."""

import argparse
from collections import Counter
import csv
import itertools
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen

TRACE_SHA = "bda00fb593b357bc8f07e43544feae598150ae891ed82c988fb339a92349de0c"
STAGES = ("multipass", "multipass_refined", "strict_profiles", "strict_profiles_refined", "candidates", "root_hogs")


def boolean(value):
    if value not in ("True", "False"):
        raise ValueError("Invalid trace boolean")
    return value == "True"


def hit(value):
    if value == "NA":
        return 0
    score = float(value)
    if not math.isfinite(score) or score <= 0:
        raise ValueError("Invalid trace hit score")
    return 1


def summarize(rows, families, transitions):
    seen = {name: set() for name in families}
    genes = {name: set().union(*(set(group["family_genes"]) for group in family["stages"]["multipass"]["groups"]))
             for name, family in families.items()}
    counts = {name: Counter() for name in families}
    for row in rows:
        name = row["refog"]
        if name not in families:
            raise ValueError("Unexpected reference family")
        left, right = row["left"], row["right"]
        if left >= right or not {left, right} <= genes[name] or (left, right) in seen[name]:
            raise ValueError("Invalid or duplicate reference pair")
        seen[name].add((left, right))
        direction = hit(row["forward_normalized_hit"]) + hit(row["reverse_normalized_hit"])
        scope = "within_species" if boolean(row["same_species"]) else "cross_species"
        flags = {stage: boolean(row[stage]) for stage in STAGES}
        for category in ("all", scope):
            counts[name][f"{category}/hits{direction}/pairs"] += 1
            for stage in STAGES:
                outcome = "together" if flags[stage] else "separated"
                counts[name][f"{category}/hits{direction}/{stage}/{outcome}"] += 1
            for before, after, _ in transitions:
                outcome = ("retained" if flags[before] and flags[after] else "lost" if flags[before]
                           else "gained" if flags[after] else "absent_both")
                counts[name][f"{category}/hits{direction}/{before}_to_{after}/{outcome}"] += 1
    for name, family in families.items():
        if seen[name] != set(itertools.combinations(sorted(genes[name]), 2)):
            raise ValueError("Incomplete reference pair universe")
        observed = counts[name]
        if sum(observed[f"all/hits{d}/pairs"] for d in range(3)) != family["possible_pairs"]:
            raise ValueError("Family pair total differs")
        if (observed["all/hits2/pairs"] != family["search"]["both_direction_pairs"]
                or sum(observed[f"all/hits{d}/pairs"] for d in (1, 2)) != family["search"]["either_direction_pairs"]):
            raise ValueError("Search marginals differ")
        for stage in STAGES:
            if sum(observed[f"all/hits{d}/{stage}/together"] for d in range(3)) != family["stages"][stage]["within_family_pairs"]:
                raise ValueError("Stage marginals differ")
    aggregate = Counter()
    for values in counts.values():
        aggregate.update(values)
    return {"families": {name: dict(sorted(values.items())) for name, values in counts.items()},
            "aggregate_membership_counts": dict(sorted(aggregate.items()))}


def run(source, output):
    if output.exists():
        raise FileExistsError(output)
    report = read_frozen(source, TRACE_SHA)
    if report["status"] != "native_pair_trace_verified_branch_labels_corrected":
        raise ValueError("Trace not admitted")
    records = [record(source), report["pair_trace"]]
    for item in records:
        check(item)
    with Path(report["pair_trace"]["path"]).open(newline="") as stream:
        result = summarize(csv.DictReader(stream, delimiter="\t"), report["families"], report["transition_specification"])
    for item in records:
        check(item)
    result.update(status="admitted_ob_search_stage_joint_counts", source=record(__file__), checked_records=records,
                  reference_inventory=report["reference_inventory"], transitions=report["transition_specification"],
                  limitations=["Development-exposed descriptive counts, not an independent accuracy test.",
                      "Includes low-certainty and within-species pairs; not official weighted benchmark recall.",
                      "Counts sum family memberships; overlapping reference labels are not deduplicated.",
                      "No direct hit does not identify a prefilter versus scoring rejection or exclude an indirect graph path.",
                      "A retained direct hit need not be an accepted RBNH edge; grouping outcomes do not identify causal edges.",
                      "Profile branch comparison is not a sequential processing step; no new parameter tuning or uncertainty claim."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.source.resolve(), args.output.resolve())
