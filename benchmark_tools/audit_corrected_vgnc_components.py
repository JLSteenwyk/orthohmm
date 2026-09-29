"""Describe corrected VGNC cross-block dependency graphs without inferring CIs."""

import argparse
from collections import Counter
import csv
from itertools import combinations
import json
from pathlib import Path

import igraph as ig

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

SOURCE_SHA = "01a51a62a536bcf759f2a0ef83c56196a92022b0bb8290448eacbfefa47d6005"
CATEGORIES = ("TP", "FP", "FN")


def read_reference(path):
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["block", "proteins", "asserted_pairs", "species"]:
            raise ValueError("Unexpected reference columns")
        result = {}
        for row in reader:
            if not row["block"] or row["block"] in result:
                raise ValueError("Empty or duplicate reference block")
            values = {k: int(row[k]) for k in ("proteins", "asserted_pairs")}
            if min(values.values()) <= 0:
                raise ValueError("Empty reference block counts")
            result[row["block"]] = values
    if not result:
        raise ValueError("Empty reference inventory")
    return result


def read_cells(path, reference):
    totals, diagonal = Counter(), Counter()
    links, seen = set(), set()
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["block_left", "block_right", *CATEGORIES]:
            raise ValueError("Unexpected count columns")
        for row in reader:
            left, right = row["block_left"], row["block_right"]
            pair = left, right
            if left not in reference or right not in reference or left > right or pair in seen:
                raise ValueError("Invalid or repeated block cell")
            seen.add(pair)
            counts = {k: int(row[k]) for k in CATEGORIES}
            if min(counts.values()) < 0 or sum(counts.values()) <= 0:
                raise ValueError("Invalid sparse count cell")
            totals.update(counts)
            if left == right:
                diagonal[left] += counts["TP"] + counts["FN"]
            else:
                if counts["TP"] or counts["FN"] or not counts["FP"]:
                    raise ValueError("Cross-block truth is not allowed")
                links.add(pair)
    if dict(diagonal) != {k: v["asserted_pairs"] for k, v in reference.items()}:
        raise ValueError("Reference truth counts changed")
    return links, {k: totals[k] for k in CATEGORIES}


def components(reference, links):
    labels = sorted(reference)
    indices = {name: i for i, name in enumerate(labels)}
    graph = ig.Graph(n=len(labels), edges=[(indices[a], indices[b]) for a, b in sorted(links)], directed=False)
    groups = [sorted(labels[i] for i in group) for group in graph.connected_components()]
    groups.sort(key=lambda g: (-len(g), g))
    counts = Counter(map(len, groups))
    touched = sum(len(g) for g in groups if len(g) > 1)
    total_truth = sum(v["asserted_pairs"] for v in reference.values())
    return dict(reference_blocks=len(labels), distinct_cross_block_links=len(links),
        connected_components=len(groups), isolated_blocks=counts[1], linked_blocks=touched,
        component_size_histogram=dict(sorted(counts.items())),
        largest_component_blocks=len(groups[0]),
        largest_component_reference_pairs=sum(reference[k]["asserted_pairs"] for k in groups[0]),
        linked_reference_pairs=sum(reference[k]["asserted_pairs"] for g in groups if len(g) > 1 for k in g),
        total_reference_pairs=total_truth,
        largest_components=[dict(blocks=len(g), reference_pairs=sum(reference[k]["asserted_pairs"] for k in g),
            proteins=sum(reference[k]["proteins"] for k in g), first_label=g[0]) for g in groups[:10]])


def audit(source):
    source_ref = record(source)
    if source_ref["sha256"] != SOURCE_SHA:
        raise ValueError("Changed corrected VGNC mapping")
    data = json.loads(source.read_text())
    refs = [source_ref, data["reference_table"], *data["checked_records"], *[m["table"] for m in data["methods"]]]
    for ref in refs:
        check(ref)
    reference = read_reference(Path(data["reference_table"]["path"]))
    if len(reference) != data["reference"]["reference_blocks"] or len(data["methods"]) != 8:
        raise ValueError("Unexpected corrected panel")
    links_by_method, methods = {}, []
    for method in data["methods"]:
        links, counts = read_cells(Path(method["table"]["path"]), reference)
        if counts != method["counts"] or len(links) != method["nonzero_cross_block_cells"]:
            raise ValueError("Counts or links differ from admitted mapping")
        if method["key"] in links_by_method:
            raise ValueError("Repeated method key")
        links_by_method[method["key"]] = links
        methods.append(dict(key=method["key"], label=method["label"], counts=counts,
                            dependency_graph=components(reference, links)))
    comparisons = [dict(methods=[a, b], dependency_graph=components(reference, links_by_method[a] | links_by_method[b]))
                   for a, b in combinations(links_by_method, 2)]
    union = components(reference, set.union(*links_by_method.values()))
    for ref in refs:
        check(ref)
    return dict(status="corrected_vgnc_component_structure_described", source=record(Path(__file__)),
        checked_records=refs, igraph_version=ig.__version__, methods=methods, pairwise_unions=comparisons,
        all_method_union=union, uncertainty_admitted=False, publication_ready=False,
        limitations=["Components use observed method-specific false-positive links, not prespecified independent sampling units.",
                     "Shared evolutionary history and dependence without scored links are not modeled.",
                     "Unioning methods preserves shared scored links but does not establish a joint sampling law.",
                     "No bootstrap, confidence interval, new scoring or changed native category is produced."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.source)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(result["all_method_union"], indent=2))
