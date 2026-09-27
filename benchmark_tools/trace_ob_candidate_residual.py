"""Localize residual historical/replay candidate differences without inference."""

import argparse
import json
import math
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.compare_installed_ob_search import partition
from benchmark_tools.audit_installed_orthobench import compare_partitions
from benchmark_tools.probe_installed_ob_clustering import write_json

SOURCES = {
    "historical": ("ob_native_replay_verification_20260916.json", "b5dcac4012a84659d6c52f179671d3f4fc10c9e04cc1283cf48635ce7fb09b79"),
    "prepared": ("orthobench_factorial_prepared_20260916.json", "5c325f4d77865e0c7571fe4bb4df0be0959977a49e1d22169f192518700f9382"),
    "replay": ("ob_dependency_replay_readback_22320.json", "973d1bf10bb0e4e0028e0762cbee36def2c57e4897aa1308c25e231b87f55df4"),
}


def event_key(event):
    sides = []
    for name in ("source_genes", "target_genes"):
        values = event[name]
        if not isinstance(values, list) or not values or any(not isinstance(g, str) or not g for g in values):
            raise ValueError("Invalid merge gene list")
        if len(set(values)) != len(values):
            raise ValueError("Duplicate merge gene")
        sides.append(tuple(sorted(values)))
    if set(sides[0]) & set(sides[1]) or type(event["iteration"]) is not int or event["iteration"] < 0:
        raise ValueError("Invalid directed merge")
    return (event["iteration"], *sides)


def compare_events(left, right):
    maps = []
    for events in (left, right):
        keyed = {event_key(e): e for e in events}
        if len(keyed) != len(events):
            raise ValueError("Duplicate semantic merge event")
        maps.append(keyed)
    a, b = maps
    common = a.keys() & b.keys()
    numerical = {}
    for field in ("support", "margin", "forward_average", "reverse_average",
                  "forward_normalized_support", "reverse_normalized_support"):
        changed, maximum = 0, 0.0
        for key in common:
            x, y = a[key][field], b[key][field]
            if math.isnan(x) or math.isnan(y):
                raise ValueError("NaN trace value")
            if x != y:
                changed += 1
                maximum = max(maximum, abs(x - y))
        numerical[field] = dict(changed_events=changed, max_absolute_difference=maximum)
    first = next((i for i, (x, y) in enumerate(zip(left, right)) if event_key(x) != event_key(y)), None)
    if first is None and len(left) != len(right):
        first = min(len(left), len(right))
    return dict(left_events=len(left), right_events=len(right), common_events=len(common),
                left_only=[a[k] for k in sorted(a.keys() - b.keys())],
                right_only=[b[k] for k in sorted(b.keys() - a.keys())],
                first_semantic_order_difference_zero_based=first,
                common_event_numerical_differences=numerical)


def trace(repo):
    records, reports = [], {}
    for key, (name, sha) in SOURCES.items():
        path = repo / "benchmark_tools/results" / name
        item = record(path)
        if item["sha256"] != sha:
            raise ValueError("Changed retained report")
        records.append(item)
        reports[key] = json.loads(path.read_text())
    pinned = {r["path"]: r for r in reports["replay"]["checked_records"]}
    def read(path):
        item = pinned[str(path)]
        check(item)
        records.append(item)
        return json.loads(path.read_text())
    root = repo / "benchmarks/work/ob_dependency_replay_v2_20260926/leiden011"
    replay, stage = read(root / "replay.json"), read(root / "stage_report.json")
    names_path = repo / "benchmarks/work/publication_installed_orthobench_20260926/inference/orthohmm_working_res/high_sensitivity_checkpoint/gene_names.txt"
    records.append(pinned[str(names_path)])
    check(records[-1])
    universe = set(names_path.read_text().splitlines())
    stages = []
    for entry in replay["stages"]:
        old = reports["historical"]["stages"][entry["label"]]["output"]
        new = entry["output"]
        if pinned[new["path"]] != new:
            raise ValueError("Unbound replay stage")
        for item in (old, new):
            check(item)
            records.append(item)
        comparison = compare_partitions(partition(Path(old["path"]), universe), partition(Path(new["path"]), universe))
        stages.append(dict(stage=entry["label"], byte_equal=old["sha256"] == new["sha256"], **comparison))
    old_arm = reports["prepared"]["candidate_arms"]["p1_c1"]
    if old_arm["seed_partition"] != reports["historical"]["stages"]["strict_profiles_refined"]["output"]:
        raise ValueError("Historical candidate seed identity differs")
    if old_arm["expansion"]["parameters"] != stage["candidates"]["parameters"]:
        raise ValueError("Candidate parameters differ")
    old_candidate, new_candidate = old_arm["candidate_partition"], stage["candidate_partition"]
    old_trace = old_arm["membership_constraints"]
    new_trace = pinned[stage["candidates"]["merge_trace_sidecar"]]
    for item in (old_candidate, new_candidate, old_trace, new_trace):
        check(item)
        records.append(item)
    a, b = (partition(Path(item["path"]), universe) for item in (old_candidate, new_candidate))
    differences = dict(historical_only_groups=[sorted(g) for g in sorted(set(a)-set(b), key=lambda g: tuple(sorted(g)))],
                       replay_only_groups=[sorted(g) for g in sorted(set(b)-set(a), key=lambda g: tuple(sorted(g)))])
    changed = sorted(set().union(*(set(a)-set(b))))
    event_comparison = compare_events(json.loads(Path(old_trace["path"]).read_text()),
                                      json.loads(Path(new_trace["path"]).read_text()))
    for item in records:
        check(item)
    return dict(status="candidate_residual_localized", source=record(__file__), checked_records=records,
                stages=stages, candidates=compare_partitions(a, b), changed_genes=changed,
                changed_groups=differences, merge_events=event_comparison,
                first_observed_partition_difference="candidate_expansion" if all(s["byte_equal"] for s in stages) and set(a)!=set(b) else "unresolved",
                accuracy_evaluated=False, native_inference_run=False,
                limitations=["Post-hoc stage localization, not a causal score/order intervention.",
                    "Semantic merge identity includes iteration and directed source/target gene sets, not numeric cluster labels.",
                    "First different event order is not necessarily the first causally consequential event.",
                    "No final phylogeny/F1 attribution or scientific default change."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    write_json(args.output, trace(args.repo.resolve()))
