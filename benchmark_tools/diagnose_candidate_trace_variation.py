"""Read retained satellite traces and test cap/tie sensitivity in the frozen engine."""

import argparse
from collections import Counter, defaultdict
import importlib.util
import json
import math
from pathlib import Path
import platform
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.collect_factorial_resources import FIXED_INPUTS, require
from benchmark_tools.prepare_ob_candidate_neighborhood import record


LINKAGE = ("benchmark_tools/results/factorial_native_resource_linkage_20261004/linkage.json",
           "896063b253c06f8691068231d241f8cc6eb2c7f588baee4efa1a6e72a42d4cac")
TRACE_SHA = {8: "203629b9e360d6c7b49d6af2a1255d2c0b74320aae6694e0c4fde462a6006393",
             16: "335084b5f7f0352e46ba29c05e61457ab752b80e03e0ed1d56dd62ab41ac9135",
             24: "c9eb544b595f8426b747fa699cfc848e7c02e6dd11518dbc22957dcf19b47c5c"}
FEATURES = ("support", "margin", "species_overlap_fraction", "forward_hits", "reverse_hits",
            "forward_average", "reverse_average", "forward_maximum", "reverse_maximum",
            "forward_coverage", "reverse_coverage", "forward_normalized_support", "reverse_normalized_support")


def key(row):
    return row["iteration"], tuple(sorted(row["source_genes"])), tuple(sorted(row["target_genes"]))


def validate_trace(trace):
    require(isinstance(trace, list) and bool(trace), "Empty or malformed trace")
    observed = {}
    for row in trace:
        require(type(row["iteration"]) is int and row["iteration"] in (0, 1), "Wrong trace iteration")
        for role in ("source", "target"):
            genes = row[role + "_genes"]
            require(isinstance(genes, list) and bool(genes)
                    and all(isinstance(g, str) and g and g.strip() == g for g in genes)
                    and len(set(genes)) == len(genes)
                    and type(row[role + "_size"]) is int and row[role + "_size"] == len(genes),
                    "Invalid trace membership or size")
            require(type(row[role + "_cluster"]) is int and row[role + "_cluster"] >= 0,
                    "Invalid cluster identifier")
        require(not set(row["source_genes"]) & set(row["target_genes"]), "Overlapping trace groups")
        for name in FEATURES:
            value = row[name]
            require(not isinstance(value, bool) and isinstance(value, (float, int)) and value >= 0
                    and (math.isfinite(value) or name == "margin" and value == math.inf),
                    "Invalid trace feature: " + name)
        identity = key(row)
        require(identity not in observed, "Duplicate semantic merge identity")
        observed[identity] = row
    return observed


def scalar(value):
    return "positive_infinity" if value == math.inf else value


def selection(row):
    return {k: scalar(row[k]) for k in ("source_genes", "source_cluster", "target_cluster", "support", "margin")}


def trace_comparison(original, current, cap):
    old, new = validate_trace(original), validate_trace(current)
    common = sorted(old.keys() & new.keys())
    require(bool(common), "No comparable accepted merges")
    features = {}
    for name in FEATURES:
        delta = [0.0 if old[k][name] == new[k][name] else abs(old[k][name] - new[k][name]) for k in common]
        features[name] = {"changed_values": sum(v != 0 for v in delta),
                          "maximum_absolute_delta": scalar(max(delta))}
    def anchors(trace):
        groups = defaultdict(list)
        for row in trace:
            groups[(row["iteration"], tuple(sorted(row["target_genes"])))].append(row)
        return groups
    old_anchors, new_anchors = anchors(original), anchors(current)
    changed = []
    for identity in sorted(old_anchors.keys() & new_anchors.keys()):
        a, b = old_anchors[identity], new_anchors[identity]
        sa = {tuple(sorted(row["source_genes"])) for row in a}
        sb = {tuple(sorted(row["source_genes"])) for row in b}
        if sa == sb:
            continue
        only_a = [row for row in a if tuple(sorted(row["source_genes"])) not in sb]
        only_b = [row for row in b if tuple(sorted(row["source_genes"])) not in sa]
        supports = [row["support"] for row in only_a + only_b]
        margins = [row["margin"] for row in only_a + only_b]
        changed.append({"iteration": identity[0], "target_genes": list(identity[1]),
                        "original_attachments": len(a), "native_attachments": len(b),
                        "both_at_attachment_cap": len(a) == len(b) == cap,
                        "original_only_sources": sorted(sa - sb), "native_only_sources": sorted(sb - sa),
                        "changed_selection_support_spread": max(supports) - min(supports),
                        "changed_selection_margin_spread": scalar(0.0 if max(margins) == min(margins)
                                                                    else max(margins) - min(margins)),
                        "original_selected": [selection(row) for row in a],
                        "native_selected": [selection(row) for row in b]})
    def differences(keys):
        return {str(k): v for k, v in sorted(Counter(identity[0] for identity in keys).items())}
    return {"original_merges": len(old), "native_merges": len(new), "common_semantic_merges": len(common),
            "original_only_merges_by_iteration": differences(old.keys() - new.keys()),
            "native_only_merges_by_iteration": differences(new.keys() - old.keys()),
            "common_source_cluster_ids_changed": sum(old[k]["source_cluster"] != new[k]["source_cluster"] for k in common),
            "common_target_cluster_ids_changed": sum(old[k]["target_cluster"] != new[k]["target_cluster"] for k in common),
            "common_feature_differences": features, "common_anchors_with_changed_sources": changed,
            "original_only_anchor_identities": len(old_anchors.keys() - new_anchors.keys()),
            "native_only_anchor_identities": len(new_anchors.keys() - old_anchors.keys())}


def load_engine(path):
    spec = importlib.util.spec_from_file_location("candidate_diagnostic_frozen_refinement", path)
    module = importlib.util.module_from_spec(spec)
    previous = sys.modules.get(spec.name)
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        if previous is None:
            del sys.modules[spec.name]
        else:
            sys.modules[spec.name] = previous
    return module


def controlled_fixture(engine, parameters):
    # Nine equally supported satellites exceed the four-per-round, two-round cap.
    clusters = [[gene] for gene in range(9)] + [list(range(9, 19))]
    queries, targets, scores = [], [], []
    for satellite in range(9):
        for anchor in range(9, 19):
            queries.extend((satellite, anchor))
            targets.extend((anchor, satellite))
            scores.extend((1.0, 1.0))
    perturbed = [math.nextafter(score, math.inf) if q == 8 or t == 8 else score
                 for q, t, score in zip(queries, targets, scores)]
    cases = (("baseline", clusters, scores), ("cluster_order_only", list(reversed(clusters)), scores),
             ("one_ulp_scores_only", clusters, perturbed))
    observations = []
    for name, order, values in cases:
        trace = []
        groups, merges, relations, iterations = engine.merge_supported_satellite_candidate_clusters(
            order, queries, targets, values, [i % 3 for i in range(19)], merge_trace=trace, **parameters)
        canonical = sorted([sorted(group) for group in groups])
        observations.append({"case": name, "partition": canonical, "merges": merges,
                             "directed_relations": relations, "iterations": iterations,
                             "unattached_satellites": [group[0] for group in canonical if len(group) == 1],
                             "selections": [{"iteration": row["iteration"], "source_genes": list(row["source_genes"]),
                                             "support": row["support"], "margin": scalar(row["margin"])} for row in trace],
                             "input_score_maximum_change": max(abs(a - b) for a, b in zip(scores, values)),
                             "partition_equal_to_baseline": None if not observations else canonical == observations[0]["partition"]})
    observations[0]["partition_equal_to_baseline"] = True
    return {"genes": 19, "initial_groups": 10, "directed_hits": len(scores),
            "parameters": parameters, "cases": observations,
            "scope": "Synthetic sufficient mechanism using frozen production function; not historical-input causal replay or accuracy evidence"}


def collect(root):
    root = Path(root).resolve()
    inputs = {}
    def check(ref):
        require(record(ref["path"]) == ref, "Input checksum mismatch: " + ref["path"])
        inputs[ref["path"]] = ref
        return Path(ref["path"])
    def fixed(spec):
        ref = record(root / spec[0])
        require(ref["sha256"] == spec[1], "Wrong fixed diagnostic input")
        return json.loads(check(ref).read_text())
    linkage = fixed(LINKAGE)
    prep = fixed(FIXED_INPUTS["ob_preparation"])
    arm = prep["candidate_arms"]["p1_c1"]
    original = json.loads(check(arm["membership_constraints"]).read_text())
    engine_ref = next(ref for ref in prep["core_sources"] if ref["path"].endswith("/orthohmm/refinement.py"))
    engine = load_engine(check(engine_ref))
    points = []
    for index, digest in TRACE_SHA.items():
        path = root / ("benchmarks/results/threadripper_scaling_v1/run_%02d/orthohmm_satellite_v2/"
                       "orthohmm_working_res/phylogeny_candidate_merges.json" % index)
        ref = record(path)
        require(ref["sha256"] == digest, "Native trace identity changed")
        trace = json.loads(check(ref).read_text())
        point = next(row for row in linkage["points"] if row["index"] == index)
        compared = trace_comparison(original, trace, arm["expansion"]["parameters"]["max_satellites_per_anchor"])
        require(compared["original_merges"] == compared["native_merges"] == arm["expansion"]["merges"],
                "Trace count differs from recorded expansion")
        points.append({"index": index, "trace": ref, "comparison": compared,
                       "candidate_genes_in_changed_groups": point["candidate_partition"]["genes_in_changed_groups"],
                       "final_genes_in_changed_groups": point["partition"]["genes_in_changed_groups"]})
    fixture = controlled_fixture(engine, arm["expansion"]["parameters"])
    for ref in list(inputs.values()):
        check(ref)
    return {"schema": "candidate_trace_variation_v1", "status": "retained_traces_and_frozen_fixture_checked",
            "inputs": sorted(inputs.values(), key=lambda ref: ref["path"]), "points": points,
            "fixture": fixture, "diagnostic_runtime": {"python": platform.python_version(), "numpy": np.__version__},
            "source": record(Path(__file__)), "native_inference_or_scoring_repeated": False,
            "frozen_method_modified": False, "publication_ready": False,
            "limitations": [
                "All accepted merge records are compared by iteration and named memberships, not only group counts or numeric labels. Rejected candidates and complete hit arrays are not replayed.",
                "Common accepted cluster IDs remain unchanged; the cluster-order fixture is a sensitivity result, not proof of historical relabeling.",
                "Capped selection with tied/near-tied evidence is consistent with the retained differences. The fixture proves this is a sufficient mechanism, not a causal explanation of every historical difference or of score-bit provenance.",
                "The 19-gene fixture uses the frozen source and original satellite parameters under the recorded reporting runtime, not full native benchmark execution or the original native deployment.",
                "Structural reference exclusion/equality is inherited from the pinned linkage report, not new benchmark scoring. Native pair files, trees and alignments remain unexamined here.",
                "Do not change the frozen scientific method, round scores opportunistically, choose matching repeats or claim universal invariance. A future determinism repair needs prospective numerical/ordering tests and new independent validation if assignments change.",
                "Existing timing observations retain unknown potentially method-dependent shared-host contention; this diagnostic neither corrects timing nor attributes variation to contention.",
            ]}


def render(report):
    lines = ["# Candidate Trace Variation", "", "Retained accepted merges; no new benchmark inference or accuracy scoring.", "",
             "| Run | Common merges / 8440 | Original-only merges, round 0 / 1 | Changed common anchors | All changed common anchors at cap | Max common support delta |",
             "| ---: | ---: | --- | ---: | --- | ---: |"]
    for point in report["points"]:
        c = point["comparison"]
        counts = c["original_only_merges_by_iteration"]
        lines.append("| %d | %d | %d / %d | %d | %s | %.3g |" % (
            point["index"], c["common_semantic_merges"], counts.get("0", 0), counts.get("1", 0),
            len(c["common_anchors_with_changed_sources"]),
            all(a["both_at_attachment_cap"] for a in c["common_anchors_with_changed_sources"]),
            c["common_feature_differences"]["support"]["maximum_absolute_delta"]))
    lines += ["", "## Frozen-Engine Fixture", "",
              "Nineteen genes; nine equally supported singleton satellites and a ten-gene anchor; original satellite parameters.",
              "| Case | Unattached satellite | Merges | Rounds | Largest score perturbation |",
              "| --- | ---: | ---: | ---: | ---: |"]
    for row in report["fixture"]["cases"]:
        lines.append("| %s | %s | %d | %d | %.3g |" % (row["case"], row["unattached_satellites"],
                     row["merges"], row["iterations"], row["input_score_maximum_change"]))
    lines += ["", "## Limits", "", *["- " + text for text in report["limitations"]]]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output-directory", required=True, type=Path)
    args = parser.parse_args()
    require(not args.output_directory.exists(), "Refusing existing diagnostic output")
    report = collect(args.root)
    args.output_directory.mkdir(parents=True)
    (args.output_directory / "diagnostic.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    (args.output_directory / "diagnostic.md").write_text(render(report))
    print(json.dumps({"status": report["status"], "traces": 4, "fixture_cases": 3,
                      "native_inference_repeated": False}))


if __name__ == "__main__":
    main()
