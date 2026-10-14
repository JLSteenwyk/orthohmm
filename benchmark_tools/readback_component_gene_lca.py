"""Independently rebuild LCA buckets from native node ancestry and raw split support."""

import argparse
from collections import defaultdict
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.compare_qfo_scored_pairs import read_scores
from benchmark_tools.diagnose_hillclimb_pair_confidence import read_confidence
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json
from orthohmm.phylogeny import parse_gene_tree
from qfo_benchmark.og_to_pairwise import _strip_to_uniprot


def split_key(genes, universe):
    side = frozenset(genes)
    if not side or not side < universe:
        raise ValueError("Require a nonempty proper branch split")
    return tuple(sorted((tuple(sorted(side)), tuple(sorted(universe - side)))))


def raw_support(raw, genes):
    tree = parse_gene_tree(raw)
    names = {f"g{index:08d}": gene for index, gene in enumerate(sorted(genes))}
    for leaf in tree.leaf_node_iter():
        leaf.taxon.label = names[str(leaf.taxon.label)]
    observed = [str(leaf.taxon.label) for leaf in tree.leaf_node_iter()]
    universe = frozenset(genes)
    if set(observed) != universe or len(observed) != len(universe):
        raise ValueError("Raw tree membership differs from native family")
    values = {}
    for node in tree.preorder_node_iter():
        if node.parent_node is None or node.is_leaf():
            continue
        try:
            support = float(node.label)
        except (TypeError, ValueError):
            continue
        if 1 < support <= 100:
            support /= 100
        if not math.isfinite(support) or not 0 <= support <= 1:
            continue
        key = split_key((str(leaf.taxon.label) for leaf in node.leaf_iter()), universe)
        if key in values and values[key] != support:
            raise ValueError("Conflicting measured raw split support")
        values[key] = support
    return universe, values


def lca_from_parent_table(pair, nodes):
    if pair[0] == pair[1] or any(gene not in nodes or nodes[gene]["event"] != "leaf" for gene in pair):
        raise ValueError("Require distinct native leaf genes")
    def chain(gene):
        seen = set()
        while gene:
            if gene not in nodes or gene in seen:
                raise ValueError("Invalid native ancestry")
            seen.add(gene)
            yield gene
            gene = nodes[gene]["parent_node_id"]
    left = set(chain(pair[0]))
    for node in chain(pair[1]):
        if node in left:
            return nodes[node]
    raise ValueError("Disconnected native gene ancestry")


def reconstructed_feature(pair, nodes, universe, values):
    node = lca_from_parent_table(pair, nodes)
    is_root = not node["parent_node_id"]
    support = None if is_root else values.get(split_key(node["genes"].split(","), universe))
    overlap = int(node["species_overlap_count"])
    if overlap < 0 or node["mapping_conflict"] not in {"true", "false"}:
        raise ValueError("Invalid native event evidence")
    conflict = node["mapping_conflict"] == "true"
    event = "duplication" if overlap or conflict else "speciation"
    pair_event = "duplication" if overlap else "uncertain" if conflict else "speciation"
    if (node["event"], node["pair_event"]) != (event, pair_event):
        raise ValueError("Native node table differs from fixed reconciliation rules")
    confidence = "high" if pair_event == "speciation" or (pair_event == "duplication" and (overlap >= 2 or (support or 0.) >= .9)) else "medium"
    if support is None:
        label = "unavailable"
    else:
        labels = ("below_0_5", "0_5_to_0_7", "0_7_to_0_9", "at_least_0_9")
        label = labels[sum(support >= cut for cut in (.5, .7, .9))]
    return ("reconciled", event, pair_event, conflict, "at_least_2" if overlap >= 2 else str(overlap), confidence, is_root, label)


def run(report_path, output):
    if output.exists() or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Require new CPU-only independent LCA readback")
    report = json.loads(report_path.read_text())
    if report["status"] != "admitted_component_functional_changes_localized_by_corrected_gene_lca":
        raise ValueError("LCA diagnostic is not complete")
    scheduler = subprocess.check_output(["sacct", "-j", report["job_id"], "-X", "--format=JobID,State,ExitCode", "--parsable2"], text=True)
    if scheduler.splitlines() != ["JobID|State|ExitCode", report["job_id"] + "|COMPLETED|0:0"]:
        raise ValueError("Diagnostic producer is not successfully terminal")
    evidence = [record(report_path), *report["checked"], *report["sources"]]
    for pin in evidence:
        check(pin)
    pins = {Path(pin["path"]): pin for pin in report["checked"]}
    functional_path = next(path for path in pins if path.name == "large_copy_components_functional_diagnostic_20261013.json")
    functional = json.loads(functional_path.read_text())
    values = {}
    for metric in ("GO", "EC"):
        values[metric] = []
        for cell in ("control", "components500"):
            item = next(row for row in functional["raw_summaries"] if row["metric"] == metric and row["cell"] == cell)
            if Path(item["raw"]["path"]) not in pins:
                raise ValueError("Unpinned functional score source")
            scores = read_scores(Path(item["raw"]["path"]), metric)
            if {"pairs": len(scores), "sum_millionths": sum(scores.values())} != report["raw_counts"][metric + ":" + cell]:
                raise ValueError("Functional score total differs from diagnostic")
            values[metric].append(scores)
    gene_path = next(path for path in pins if path.name == "gene_names.txt")
    names = gene_path.read_text().splitlines()
    normalized = {_strip_to_uniprot(name): name for name in names}
    if len(normalized) != len(names) or len(set(names)) != len(names):
        raise ValueError("Non-injective gene/accession universe")
    watched = set().union(*(set(scores) for pair in values.values() for scores in pair))
    full_pairs = {pair: tuple(sorted(normalized[accession] for accession in pair)) for pair in watched}
    needed = set().union(*(set(pair) for pair in full_pairs.values()))
    candidate_path = next(path for path in pins if path.name == "orthohmm_edges_clustered.txt")
    membership, selected_genes, visited = {}, {}, set()
    with candidate_path.open() as stream:
        for index, line in enumerate(stream):
            genes = line.split()
            if not genes or len(genes) != len(set(genes)) or visited.intersection(genes):
                raise ValueError("Invalid frozen candidate partition")
            visited.update(genes)
            if needed.intersection(genes):
                family = f"Family{index:07d}"
                selected_genes[family] = tuple(sorted(genes))
                membership.update((gene, family) for gene in genes if gene in needed)
    if visited != set(names) or set(membership) != needed:
        raise ValueError("Candidate partition lost required genes")
    selected_pairs, features = defaultdict(list), {}
    for pair, (a,b) in full_pairs.items():
        if membership[a] == membership[b]:
            selected_pairs[membership[a]].append(pair)
        else:
            features[pair] = ("cross_candidate", "unavailable", "unavailable", False, "unavailable", "unavailable", False, "unavailable")
    node_path = next(path for path in pins if path.name == "orthohmm_reconciliation_nodes.tsv")
    phylo = node_path.parent
    nodes = defaultdict(dict)
    with node_path.open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            family = row["source_family"]
            if family in selected_pairs:
                if row["node_id"] in nodes[family]:
                    raise ValueError("Duplicate native node identifier")
                nodes[family][row["node_id"]] = row
    family_rows = []
    for family, pairs in sorted(selected_pairs.items()):
        genes = selected_genes[family]
        if family not in nodes:
            for pair in pairs:
                features[pair] = ("native_bypass", "unavailable", "unavailable", False, "unavailable", "unavailable", False, "unavailable")
            status = "native_bypass"
        else:
            raw_path = phylo / "gene_trees" / (family + ".raw.nwk")
            checkpoint_path = phylo / "checkpoints" / (family + ".json")
            if raw_path not in pins or checkpoint_path not in pins:
                raise ValueError("Relevant native family has no pinned checkpoint")
            checkpoint = json.loads(checkpoint_path.read_text())
            if tuple(sorted(checkpoint["genes"])) != genes or checkpoint["raw_tree_sha256"] != pins[raw_path]["sha256"]:
                raise ValueError("Relevant raw tree differs from active family checkpoint")
            universe, supports = raw_support(raw_path.read_text(), genes)
            if sum(not node["parent_node_id"] for node in nodes[family].values()) != 1:
                raise ValueError("Native family has multiple roots")
            for pair in pairs:
                features[pair] = reconstructed_feature(full_pairs[pair], nodes[family], universe, supports)
            status = "support_and_rooting_verified"
        family_rows.append({"family": family, "genes": len(genes), "watched_pairs": len(pairs), "status": status})
    if family_rows != report["family_rows"] or len(watched) != report["watched_pairs"]:
        raise ValueError("Independent family/route inventory differs")
    native_path = next(path for path in pins if path.name == "orthohmm_pairwise_orthologs_confidence.tsv")
    right_pairs = set().union(*(set(pair[1]) for pair in values.values()))
    expected_count = sum(json.loads(functional_path.read_text())["native_confidence_counts"][1].values())
    native_confidence, _ = read_confidence(native_path, right_pairs, expected_count)
    rebuilt = []
    for metric, (left, right) in values.items():
        buckets = defaultdict(lambda: [0,0])
        for pair in left.keys() | right.keys():
            if pair in left and pair in right:
                if left[pair] != right[pair]:
                    raise ValueError("Shared functional score changed")
                category, value = "shared", right[pair]
            elif pair in right:
                category, value = "added", right[pair]
            else:
                category, value = "removed", left[pair]
            key = (category, *features[pair], native_confidence.get(pair, "not_assessed_in_component_scored_sets"))
            buckets[key][0] += 1
            buckets[key][1] += value
        for key, (count, total) in sorted(buckets.items()):
            category, route, event, pair_event, conflict, overlap, confidence, root, support, native = key
            rebuilt.append({"metric": metric, "category": category, "family_route": route, "lca_event": event,
                "lca_pair_event": pair_event, "lca_mapping_conflict": conflict, "lca_species_overlap": overlap,
                "lca_event_confidence": confidence, "lca_is_root": root, "corrected_lca_support_bin": support,
                "final_native_scored_pair_confidence": native, "scored_relations": count,
                "score_sum_millionths": total, "serialized_mean": total / count / 1e6})
    if rebuilt != report["rows"]:
        raise ValueError("Independent native-ancestry/raw-support buckets differ from diagnostic")
    for pin in evidence:
        check(pin)
    write_json(output, {"status": "all_gene_lca_buckets_independently_reconstructed_from_native_ancestry",
        "producer": record(report_path), "reader": record(__file__), "scheduler": scheduler,
        "families": len(family_rows), "watched_pairs": len(watched), "rows": rebuilt,
        "scope": "Independent LCA reconstruction from retained native parent-node records and original raw split supports; no rooting or reconciliation calls. Same exposed data, not independent biological validation or new accuracy scores.",
        "defaults_promoted": False, "goal_completed": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(args.report.resolve(), args.output.resolve())
