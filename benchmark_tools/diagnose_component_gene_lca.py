"""Describe component GO/EC changes by corrected native gene-LCA evidence."""

import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.compare_qfo_scored_pairs import read_scores
from benchmark_tools.diagnose_component_pair_structure import partition_membership
from benchmark_tools.diagnose_hillclimb_pair_confidence import read_confidence
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.validate_gene_support_transport import independent_tree_inventory, matching_lengths
from orthohmm.phylogeny import parse_gene_tree, reconcile_gene_tree
from orthohmm.phylogeny_pipeline import _gene_tokens, _root_external_gene_tree
from qfo_benchmark.og_to_pairwise import _strip_to_uniprot


def support_bin(value):
    if value is None:
        return "unavailable"
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError("Invalid normalized LCA support")
    return "below_0_5" if value < .5 else "0_5_to_0_7" if value < .7 else "0_7_to_0_9" if value < .9 else "at_least_0_9"


def node_index(nodes):
    indexed = {node.node_id: node for node in nodes}
    if len(indexed) != len(nodes) or sum(node.parent_node_id is None for node in nodes) != 1:
        raise ValueError("Invalid reconciled node identities/root")
    ancestors = {}
    for node in nodes:
        chain, current = [], node.node_id
        while current is not None:
            if current not in indexed or current in chain:
                raise ValueError("Invalid reconciled ancestry")
            chain.append(current)
            current = indexed[current].parent_node_id
        ancestors[node.node_id] = tuple(chain)
    return indexed, ancestors


def pair_features(pair, indexed, ancestors):
    a, b = pair
    if a == b or a not in indexed or b not in indexed or indexed[a].event != "leaf" or indexed[b].event != "leaf":
        raise ValueError("Require two distinct reconciled leaf genes")
    left = set(ancestors[a])
    common = next(node for node in ancestors[b] if node in left)
    node = indexed[common]
    root = node.parent_node_id is None
    return (node.event, node.pair_event, node.mapping_conflict,
            "at_least_2" if node.species_overlap_count >= 2 else str(node.species_overlap_count),
            node.event_confidence, root, support_bin(None if root else node.branch_support))


def category_and_score(pair, left, right):
    if pair in left and pair in right:
        if left[pair] != right[pair]:
            raise ValueError("Shared serialized functional score changed")
        return "shared", right[pair]
    return ("added", right[pair]) if pair in right else ("removed", left[pair])


def verify_totals(rows, raw_counts):
    for metric in ("GO", "EC"):
        for cell, categories in (("control", {"shared", "removed"}), ("components500", {"shared", "added"})):
            selected = [row for row in rows if row["metric"] == metric and row["category"] in categories]
            observed = {"pairs": sum(row["scored_relations"] for row in selected),
                        "sum_millionths": sum(row["score_sum_millionths"] for row in selected)}
            if observed != raw_counts[metric + ":" + cell]:
                raise ValueError("LCA buckets changed original functional counts/sums")


def run(repo, output, source_revision):
    if output.exists() or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Require new CPU-only LCA diagnostic output")
    started = time.monotonic()
    executor = Path(__file__).resolve().parent.parent
    revision = subprocess.check_output(["git", "rev-parse", "--verify", source_revision + "^{commit}"], cwd=repo, text=True).strip()
    source_names = ("orthohmm/phylogeny.py", "orthohmm/phylogeny_pipeline.py",
        "benchmark_tools/diagnose_component_gene_lca.py", "benchmark_tools/validate_gene_support_transport.py",
        "benchmark_tools/hillclimb_component_gene_lca_diagnostic_20261014.json")
    sources = [record(executor / name) for name in source_names]
    for name, pin in zip(source_names, sources):
        if hashlib.sha256(subprocess.check_output(["git", "show", revision + ":" + name], cwd=repo)).hexdigest() != pin["sha256"]:
            raise ValueError("Wrong frozen LCA executor source")
    base = repo / "benchmark_tools/results"
    repair_path = base / "gene_support_transport_repair_admission_20261014.json"
    repair = json.loads(repair_path.read_text())
    if repair["status"] != "opt_in_support_repair_verified_on_prespecified_native_sample":
        raise ValueError("Support transport has not been admitted")
    check(repair["primary"])
    check(repair["independent_repeat"])
    admission_path = repo / "benchmarks/results/component_recovery_admission_r1_v1_20261013/report.json"
    admission = json.loads(admission_path.read_text())
    if admission["status"] != "terminal_component_candidates_and_raw_scores_independently_verified":
        raise ValueError("Component native result not admitted")
    checked = [record(repair_path), repair["primary"], repair["independent_repeat"], record(admission_path),
               admission["execution_plan"], admission["original_failed_producer"]]
    for pin in checked:
        check(pin)
    repeated_repair = json.loads(Path(repair["independent_repeat"]["path"]).read_text())
    plan = json.loads(Path(admission["execution_plan"]["path"]).read_text())
    producer = json.loads(Path(admission["original_failed_producer"]["path"]).read_text())
    candidate_pin = producer["structural_controls"]["components500"]["candidate"]["candidate_partition"]
    gene_path = Path(plan["data"]["checkpoint"]) / "gene_names.txt"
    gene_pin = next(pin for pin in plan["inputs"] if Path(pin["path"]) == gene_path)
    checked.extend((candidate_pin, gene_pin))
    check(candidate_pin)
    check(gene_pin)
    names = gene_path.read_text().splitlines()
    if len(names) != 984137 or names != sorted(set(names)):
        raise ValueError("Changed original gene universe")
    accessions = {_strip_to_uniprot(gene): index for index, gene in enumerate(names)}
    if len(accessions) != len(names):
        raise ValueError("Non-injective gene accession mapping")
    membership, _ = partition_membership(Path(candidate_pin["path"]), names)
    functional_receipt_path = base / "large_copy_components_functional_diagnostic_receipt_20261013.json"
    receipt = json.loads(functional_receipt_path.read_text())
    check(receipt["diagnostic"])
    functional = json.loads(Path(receipt["diagnostic"]["path"]).read_text())
    checked.extend((record(functional_receipt_path), receipt["diagnostic"]))
    values, raw_counts = {}, {}
    for metric in ("GO", "EC"):
        values[metric] = []
        for cell in ("control", "components500"):
            row = next(row for row in functional["raw_summaries"] if row["metric"] == metric and row["cell"] == cell)
            check(row["raw"])
            checked.append(row["raw"])
            scores = read_scores(Path(row["raw"]["path"]), metric)
            if len(scores) != row["pairs"] or abs(sum(scores.values()) / len(scores) / 1e6 - row["admitted_mean"]) > 5.051e-7:
                raise ValueError("Raw functional endpoint changed")
            values[metric].append(scores)
            raw_counts[metric + ":" + cell] = {"pairs": len(scores), "sum_millionths": sum(scores.values())}
    watched = set().union(*(set(scores) for pair in values.values() for scores in pair))
    pair_genes = {pair: tuple(sorted((names[accessions[pair[0]]], names[accessions[pair[1]]]))) for pair in watched}
    families = defaultdict(list)
    cross_family = set()
    for pair in sorted(watched):
        a, b = (accessions[gene] for gene in pair)
        if membership[a] == membership[b]:
            families[f"Family{int(membership[a]):07d}"].append(pair)
        else:
            cross_family.add(pair)
    selected_genes = {}
    with Path(candidate_pin["path"]).open() as stream:
        for index, line in enumerate(stream):
            if not line.strip():
                raise ValueError("Unexpected blank native candidate row")
            family = f"Family{index:07d}"
            if family in families:
                selected_genes[family] = tuple(sorted(line.split()))
    if selected_genes.keys() != families.keys():
        raise ValueError("Watched candidate family missing")
    phylo = Path(admission["products"]["orthohmm_root_hogs.tsv"]["path"]).parent
    node_path = phylo / "orthohmm_reconciliation_nodes.tsv"
    support_audit_path = base / "gene_support_transport_audit_20261014.json"
    support_audit_pin = next(pin for pin in repeated_repair["checked"] if Path(pin["path"]) == support_audit_path)
    check(support_audit_pin)
    checked.append(support_audit_pin)
    support_audit = json.loads(support_audit_path.read_text())
    node_pin = next(pin for pin in support_audit["checked"] if Path(pin["path"]) == node_path)
    check(node_pin)
    checked.append(node_pin)
    owners, active = defaultdict(dict), set()
    with node_path.open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            family = row["source_family"]
            active.add(family)
            if family in families and row["event"] == "leaf":
                gene = row["genes"]
                if gene in owners[family] or gene != row["node_id"]:
                    raise ValueError("Invalid native leaf ownership")
                owners[family][gene] = row["species"]
    manifest_path = phylo / "provenance_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if manifest["root_duplication_rule"] != "species_overlap" or manifest["pair_orthology_rule"] != "positive_paralogy":
        raise ValueError("Native rules changed")
    species_pin = record(phylo / "species_tree.rooted.nwk")
    if species_pin["sha256"] != manifest["species_tree_sha256"]:
        raise ValueError("Species tree differs from native provenance")
    checked.extend((record(manifest_path), species_pin))
    species_tree = parse_gene_tree(Path(species_pin["path"]).read_text())
    native_pin = admission["products"]["orthohmm_pairwise_orthologs_confidence.tsv"]
    check(native_pin)
    checked.append(native_pin)
    right_watched = set().union(*(set(pair[1]) for pair in values.values()))
    native_confidence, _ = read_confidence(Path(native_pin["path"]), right_watched, admission["coverage"]["pairs"])
    features = {pair: ("cross_candidate", "unavailable", "unavailable", False, "unavailable", "unavailable", False, "unavailable") for pair in cross_family}
    family_rows = []
    for family, pairs in sorted(families.items()):
        genes = selected_genes[family]
        if family not in active:
            for pair in pairs:
                features[pair] = ("native_bypass", "unavailable", "unavailable", False, "unavailable", "unavailable", False, "unavailable")
            family_rows.append({"family": family, "genes": len(genes), "watched_pairs": len(pairs), "status": "native_bypass"})
            continue
        checkpoint_path = phylo / "checkpoints" / (family + ".json")
        raw_path = phylo / "gene_trees" / (family + ".raw.nwk")
        root_path = phylo / "gene_trees" / (family + ".rooted.nwk")
        checkpoint = json.loads(checkpoint_path.read_text())
        pins = [record(path) for path in (checkpoint_path, raw_path, root_path)]
        if (checkpoint["status"] != "complete" or checkpoint["family_id"] != family
                or tuple(sorted(checkpoint["genes"])) != genes or set(owners[family]) != set(genes)
                or pins[1]["sha256"] != checkpoint["raw_tree_sha256"]
                or pins[2]["sha256"] != checkpoint["rooted_tree_sha256"]):
            raise ValueError("Relevant active native checkpoint mismatch: " + family)
        checked.extend(pins)
        tokens = _gene_tokens(genes)
        raw = parse_gene_tree(raw_path.read_text())
        for leaf in raw.leaf_node_iter():
            leaf.taxon.label = tokens[str(leaf.taxon.label)]
        fixed = _root_external_gene_tree(raw_path.read_text(), tokens, species_tree, owners[family], gene_tree_support_mode="bipartition")
        original, repaired, prior = (independent_tree_inventory(tree) for tree in (raw, fixed, parse_gene_tree(root_path.read_text())))
        if (original[0] != set(genes) or repaired[0] != original[0] or repaired[1] != original[1]
                or repaired[2] != original[2] or not matching_lengths(original[3], repaired[3])
                or not matching_lengths(prior[4], repaired[4])):
            raise ValueError("Relevant native rooting/support invariant failed: " + family)
        result = reconcile_gene_tree(fixed, species_tree, owners[family], root_duplication_rule="species_overlap", pair_orthology_rule="positive_paralogy")
        indexed, ancestry = node_index(result.nodes)
        confidence = {(a,b): label for a,b,label in result.ortholog_pair_confidence}
        for pair in pairs:
            full_pair = pair_genes[pair]
            evidence = pair_features(full_pair, indexed, ancestry)
            if pair in right_watched and confidence.get(full_pair) != native_confidence[pair]:
                raise ValueError("Current watched native pair/confidence changed")
            features[pair] = ("reconciled", *evidence)
        family_rows.append({"family": family, "genes": len(genes), "watched_pairs": len(pairs), "status": "support_and_rooting_verified"})
    if features.keys() != watched:
        raise ValueError("Diagnostic lost watched functional pairs")
    rows = []
    for metric, (left, right) in values.items():
        buckets = defaultdict(lambda: [0,0])
        for pair in left.keys() | right.keys():
            category, score = category_and_score(pair, left, right)
            if category != "removed" and pair in cross_family:
                raise ValueError("Predicted functional pair crosses native candidate families")
            key = (category, *features[pair], native_confidence.get(pair, "not_assessed_in_component_scored_sets"))
            buckets[key][0] += 1
            buckets[key][1] += score
        for key, (count, total) in sorted(buckets.items()):
            category, route, event, pair_event, conflict, overlap, confidence, root, support, native = key
            rows.append({"metric": metric, "category": category, "family_route": route, "lca_event": event,
                "lca_pair_event": pair_event, "lca_mapping_conflict": conflict, "lca_species_overlap": overlap,
                "lca_event_confidence": confidence, "lca_is_root": root, "corrected_lca_support_bin": support,
                "final_native_scored_pair_confidence": native, "scored_relations": count,
                "score_sum_millionths": total, "serialized_mean": total / count / 1e6})
    verify_totals(rows, raw_counts)
    for pin in (*checked, *sources):
        check(pin)
    usage = resource.getrusage(resource.RUSAGE_SELF)
    write_json(output, {"status": "admitted_component_functional_changes_localized_by_corrected_gene_lca",
        "source_revision": revision, "job_id": os.environ.get("SLURM_JOB_ID"), "checked": checked, "sources": sources,
        "raw_counts": raw_counts, "watched_pairs": len(watched), "family_rows": family_rows, "rows": rows,
        "wall_s": time.monotonic() - started, "cpu_s": usage.ru_utime + usage.ru_stime,
        "peak_rss_kib": usage.ru_maxrss,
        "scope": "Prespecified descriptive GO/EC eligible scored-pair analysis using saved tree reroot/reconciliation. No filtering, reference changes, new benchmark scoring, orthology posterior calibration or biological validation. Shared-host diagnostic resource use is not full inference timing.",
        "native_inference_rerun": False, "new_accuracy_scores": False, "defaults_promoted": False, "goal_completed": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--source-revision", required=True)
    args = parser.parse_args()
    run(args.repo.resolve(), args.output.resolve(), args.source_revision)
