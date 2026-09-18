"""Quantify reference overlap and cross-family scored dependencies, without CIs."""

import argparse
from collections import Counter, defaultdict
import gzip
import json
from pathlib import Path
import sys

import igraph as ig

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_vgnc_mapping import reference_data
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen

INVENTORY_SHA = "436daaa78988eb70f1ad6aa99fc71ba87155bd360d66c689f583ce8903be2472"
REFERENCE_SHA = "3f9898487b1b3dab866eaf93255bb25ab7b7699f32753828ba7d2851eafce5dc"


def blocks(truth):
    labels = sorted(set(truth.values()))
    indices = {label: i for i, label in enumerate(labels)}
    owners = defaultdict(set)
    for pair, family in truth.items():
        for protein in pair:
            owners[protein].add(family)
    edges = set()
    for families in owners.values():
        members = sorted(indices[f] for f in families)
        edges.update((members[0], other) for other in members[1:])
    components = ig.Graph(n=len(labels), edges=sorted(edges), directed=False).connected_components()
    grouped = [sorted(labels[i] for i in component) for component in components]
    mapping = {family: group[0] for group in grouped for family in group}
    return mapping, dict(family_labels=len(labels), reference_proteins=len(owners),
        shared_proteins=sum(len(v) > 1 for v in owners.values()), reference_blocks=len(grouped),
        merged_label_groups=sorted(group for group in grouped if len(group) > 1))


def summarize_raw(path, mapping):
    block_names = sorted(set(mapping.values()))
    indices = {label: i for i, label in enumerate(block_names)}
    counts, within, between = Counter(), Counter(), Counter()
    links, seen, annotations = set(), set(), {}
    with gzip.open(path, "rt") as stream:
        for line in stream:
            a, b, category, fa, fb, sa, sb = line.rstrip("\n").split("\t")
            if category not in ("TP", "FP", "FN") or a == b or fa not in mapping or fb not in mapping:
                raise ValueError("Invalid scored row")
            key = (category, *sorted((a, b)))
            if key in seen:
                raise ValueError("Duplicate scored row")
            seen.add(key)
            for protein, annotation in ((a, (fa, sa)), (b, (fb, sb))):
                if protein in annotations and annotations[protein] != annotation:
                    raise ValueError("Conflicting protein annotation")
                annotations[protein] = annotation
            left, right = mapping[fa], mapping[fb]
            counts[category] += 1
            (within if left == right else between)[category] += 1
            if left != right:
                links.add(tuple(sorted((indices[left], indices[right]))))
                if category != "FP":
                    raise ValueError("Asserted truth crosses reference-defined overlap blocks")
    components = ig.Graph(n=len(block_names), edges=sorted(links), directed=False).connected_components()
    sizes = sorted(components.sizes(), reverse=True)
    touched = {i for edge in links for i in edge}
    return dict(counts={c: counts[c] for c in ("TP", "FP", "FN")},
        within_reference_block={c: within[c] for c in ("TP", "FP", "FN")},
        between_reference_blocks={c: between[c] for c in ("TP", "FP", "FN")},
        distinct_cross_block_links=len(links), blocks_incident_to_cross_links=len(touched),
        prediction_link_components=len(sizes), largest_prediction_link_component_blocks=sizes[0] if sizes else 0,
        ten_largest_prediction_link_components=sizes[:10])


def audit(root):
    inventory_path = root / "benchmark_tools/results/vgnc_raw_label_inventory_20260917.json"
    inventory = read_frozen(inventory_path, INVENTORY_SHA)
    reference_inventory = inventory["reference_inventory"]
    refdata = read_frozen(Path(reference_inventory["path"]), reference_inventory["sha256"])
    reference = refdata["reference"]
    if reference["sha256"] != REFERENCE_SHA:
        raise ValueError("Changed frozen reference")
    checked = [record(inventory_path), reference_inventory, reference, inventory["native_scorer"],
               *[r["raw"] for r in inventory["records"]]]
    for item in checked:
        check(item)
    if [r["stage_index"] for r in inventory["records"]] != list(range(4)):
        raise ValueError("Changed stage inventory")
    truth, _ = reference_data(Path(reference["path"]))
    mapping, summary = blocks(truth)
    stages = []
    for row in inventory["records"]:
        observed = summarize_raw(Path(row["raw"]["path"]), mapping)
        if observed["counts"] != row["counts"]:
            raise ValueError("Counts differ from audited native inventory")
        stages.append(dict(stage_index=row["stage_index"], raw=row["raw"], **observed))
    for item in checked:
        check(item)
    return dict(status="vgnc_dependency_structure_described", source=record(__file__), checked_inputs=checked,
        reference=summary, reference_pairs=len(truth), stages=stages, igraph_version=ig.__version__,
        publication_ready=False, uncertainty_admitted=False,
        limitations=["Historical development-exposed stages; no corrected-release or competitor substitution.",
            "Reference overlap blocks use labels/proteins only; prediction-link components depend on each method.",
            "Neither kind of component is established as an independent biological sampling unit.",
            "No bootstrap, confidence interval or independence claim is produced.",
            "Uses previously audited raw categories; does not independently requery prediction databases."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = audit(args.root.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
