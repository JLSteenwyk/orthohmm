"""Validate SwissTrees family sufficient statistics for all retained comparators."""

import argparse
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_swiss_counts import read_raw, statistics
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

COMPARISON_SHA = "094f842ad8211450519d4238c0b3a5020a7969049d4b7306f5465d5acbf914a3"
COUNTS_SHA = "546bb5bd6957c8ea990324b79fa31f22b0ed721bc7d6b94b609ab15258f97183"
METHODS = ("orthohmm_high_sensitivity", "orthohmm_phylogeny_satellite_v2",
           "orthofinder_3_1_5_full", "orthofinder_3_1_5_sequence_only",
           "sonicparanoid_2_0_9", "proteinortho_6_3_6", "fastoma_0_3_5", "orthomcl_1_4")


def aggregate_verified(counts, participant, reported_f1):
    values = {name: statistics(row) for name, row in counts.items()}
    means = {metric: sum(v[metric] for v in values.values()) / len(values) for metric in ("PPV", "TPR")}
    means["F1"] = 2 * means["PPV"] * means["TPR"] / (means["PPV"] + means["TPR"])
    expected = {"PPV": participant["metric_y"], "TPR": participant["metric_x"], "F1": reported_f1}
    if any(not math.isclose(means[m], expected[m], rel_tol=0, abs_tol=5e-8) for m in means):
        raise ValueError("Raw family counts do not reproduce comparator endpoint")
    return values, means


def audit(repo):
    results = repo / "benchmark_tools/results"
    comparison_path = results / "publication_comparison_orthomcl_complete_20260916.json"
    count_path = results / "qfo_swiss_counts_20260917.json"
    comparison_identity, count_identity = record(comparison_path), record(count_path)
    if comparison_identity["sha256"] != COMPARISON_SHA or count_identity["sha256"] != COUNTS_SHA:
        raise ValueError("Changed frozen comparison or reference audit")
    comparison, baseline = json.loads(comparison_path.read_text()), json.loads(count_path.read_text())
    if tuple(r["key"] for r in comparison["methods"]) != METHODS:
        raise ValueError("Unexpected comparator inventory")
    families = baseline["families"]
    if len(families) != 18 or len(set(families)) != 18 or baseline["shared_represented_genes"]:
        raise ValueError("Require 18 disjoint represented families")
    checked = [comparison_identity, count_identity, baseline["reference"], baseline["native_scorer"],
               baseline["stages"][0]["raw_file"]]
    _, expected_truth, expected_members = read_raw(Path(checked[-1]["path"]), families)
    methods = []
    for row in comparison["methods"]:
        native = row["qfo"]["metric_details"]["SwissTrees F"]
        if native["axes"]["x_axis"] != "TPR" or native["axes"]["y_axis"] != "PPV":
            raise ValueError("Changed native axis semantics")
        checked.append(native["source"])
        paths = list(Path(native["source"]["path"]).parent.glob("*raw.txt.gz"))
        if len(paths) != 1:
            raise ValueError("Missing or ambiguous raw evidence")
        raw_identity = record(paths[0])
        checked.append(raw_identity)
        counts, truth, members = read_raw(paths[0], families)
        if truth != expected_truth or members != expected_members:
            raise ValueError("Comparator reference truth or family members differ")
        values, aggregate = aggregate_verified(counts, native["participant"], row["qfo"]["scores"]["SwissTrees F"])
        methods.append({"method": row["key"], "display_name": row["method"],
                        "output_semantics": row["qfo"]["output_semantics"],
                        "native_aggregate": native["source"], "raw_file": raw_identity,
                        "aggregate": aggregate,
                        "families": [{"family": f, "counts_without_prior": dict(counts[f]),
                                      "statistics_with_prior": values[f],
                                      "represented_genes": sorted(members[f])} for f in families]})
    for item in checked:
        check(item)
    return {"status": "eight_comparator_swiss_family_counts_verified", "source": record(__file__),
            "reader": record(Path(__file__).with_name("audit_qfo_swiss_counts.py")),
            "checked_inputs": checked, "families": families, "methods": methods,
            "reference_relation_count": len(expected_truth), "shared_represented_genes": {},
            "count_conversion": "native confusion count = raw one-direction relation count / 2 + 1",
            "limitations": ["Retained historical comparison outputs, not scores transferred from recovered ablations.",
                            "No paired intervals computed by this audit.",
                            "Shared history and merged predictions can correlate disjoint families.",
                            "Aggregate and reference truth verified; original execution provenance not reconstructed.",
                            "OrthoFinder sequence-only is an MCL-checkpoint diagnostic; FastOMA uses a supplied tree."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.repo)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
