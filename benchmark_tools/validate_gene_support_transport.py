"""Replay the prespecified native tree sample without alignment or tree inference."""

import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json
from orthohmm.phylogeny import parse_gene_tree, reconcile_gene_tree
from orthohmm.phylogeny_pipeline import _gene_tokens, _root_external_gene_tree


def independent_tree_inventory(tree):
    """Read split/support/length identities independently of the repair helpers."""
    leaves = [str(n.taxon.label) for n in tree.leaf_node_iter()]
    universe = set(leaves)
    if len(leaves) != len(universe) or not all(leaves):
        raise ValueError("Invalid native tree leaf universe")
    splits, labels, lengths, rooted = set(), defaultdict(set), defaultdict(float), {}
    for node in tree.preorder_node_iter():
        side = tuple(sorted(str(n.taxon.label) for n in node.leaf_iter()))
        rooted[side] = node.edge_length
        if node.parent_node is None:
            continue
        complement = tuple(sorted(universe - set(side)))
        split = tuple(sorted((side, complement)))
        splits.add(split)
        if node.edge_length is None or not math.isfinite(node.edge_length):
            raise ValueError("Unavailable/nonfinite native tree branch length")
        lengths[split] += node.edge_length
        if node.is_leaf():
            continue
        try:
            value = float(node.label)
        except (TypeError, ValueError):
            continue
        if 1 < value <= 100:
            value /= 100
        if math.isfinite(value) and 0 <= value <= 1:
            labels[split].add(value)
    return universe, splits, dict(labels), dict(lengths), rooted


def matching_lengths(left, right):
    if left.keys() != right.keys():
        return False
    return all(a == b if a is None or b is None else math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-12)
               for key, a in left.items() for b in (right[key],))


def validate_family(raw_text, saved_rooted_text, genes, species_tree, owners):
    tokens = _gene_tokens(tuple(genes))
    raw = parse_gene_tree(raw_text)
    for leaf in raw.leaf_node_iter():
        leaf.taxon.label = tokens[str(leaf.taxon.label)]
    saved = parse_gene_tree(saved_rooted_text)
    legacy = _root_external_gene_tree(raw_text, tokens, species_tree, owners)
    fixed = _root_external_gene_tree(raw_text, tokens, species_tree, owners,
                                     gene_tree_support_mode="bipartition")
    original = independent_tree_inventory(raw)
    prior = independent_tree_inventory(saved)
    repeated = independent_tree_inventory(legacy)
    repaired = independent_tree_inventory(fixed)
    if any(item[0] != set(genes) or item[1] != original[1] for item in (prior, repeated, repaired)):
        raise ValueError("Native support replay changed leaf universe or unrooted topology")
    if (prior[2] != repeated[2] or not matching_lengths(prior[4], repeated[4])
            or not matching_lengths(repeated[4], repaired[4])):
        raise ValueError("Legacy replay differs from saved rooting or repair changed rooted topology/lengths")
    if repaired[2] != original[2]:
        raise ValueError("Repair did not preserve original measured branch support")
    if not matching_lengths(original[3], repaired[3]):
        raise ValueError("Repair changed original unrooted branch lengths")
    results = [reconcile_gene_tree(tree, species_tree, owners,
                root_duplication_rule="species_overlap", pair_orthology_rule="positive_paralogy")
               for tree in (saved, legacy, fixed)]
    for result in results[1:]:
        for attribute in ("ortholog_pairs", "ortholog_pair_confidence", "root_groups"):
            if getattr(result, attribute) != getattr(results[0], attribute):
                raise ValueError("Support replay changed current family " + attribute)
    return {"genes": len(genes), "supported_splits": len(original[2]),
            "legacy_support_map_changed": prior[2] != original[2],
            "repaired_support_map_preserved": True,
            "rooted_topology_and_lengths_unchanged": True,
            "ortholog_pairs": len(results[2].ortholog_pairs),
            "root_groups": len(results[2].root_groups),
            "current_pairs_confidence_and_groups_unchanged": True}


def run(repo, output, source_revision):
    if output.exists() or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Require new CPU-only support validation output")
    revision = subprocess.check_output(
        ["git", "rev-parse", "--verify", source_revision + "^{commit}"], cwd=repo, text=True).strip()
    executor = Path(__file__).resolve().parent.parent
    source_names = (
        "orthohmm/phylogeny.py", "orthohmm/phylogeny_pipeline.py",
        "benchmark_tools/validate_gene_support_transport.py",
        "benchmark_tools/hillclimb_gene_support_transport_repair_20261014.json")
    sources = [record(executor / name) for name in source_names]
    for name, pin in zip(source_names, sources):
        blob = subprocess.check_output(["git", "show", revision + ":" + name], cwd=repo)
        if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
            raise ValueError("Executor differs from frozen source revision: " + name)
    started = time.monotonic()
    audit_path = repo / "benchmark_tools/results/gene_support_transport_audit_20261014.json"
    audit = json.loads(audit_path.read_text())
    if audit["status"] != "active_native_raw_rooted_support_positions_audited" or len(audit["rows"]) != 32:
        raise ValueError("Unexpected prespecified support sample")
    checked = [record(audit_path), *audit["checked"]]
    for pin in checked:
        check(pin)
    pins = {Path(pin["path"]): pin for pin in audit["checked"]}
    node_path = next(path for path in pins if path.name == "orthohmm_reconciliation_nodes.tsv")
    phylo = node_path.parent
    manifest_path = phylo / "provenance_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if manifest["root_duplication_rule"] != "species_overlap" or manifest["pair_orthology_rule"] != "positive_paralogy":
        raise ValueError("Native reconciliation rules changed")
    species_path = phylo / "species_tree.rooted.nwk"
    species_pin = record(species_path)
    if species_pin["sha256"] != manifest["species_tree_sha256"]:
        raise ValueError("Species tree differs from native provenance")
    checked.extend((record(manifest_path), species_pin))
    species_tree = parse_gene_tree(species_path.read_text())
    if {str(n.taxon.label) for n in species_tree.leaf_node_iter()} != set(manifest["species_tree_taxa"]):
        raise ValueError("Native species-tree universe differs from provenance")
    families = {row["family"] for row in audit["rows"]}
    if len(families) != 32:
        raise ValueError("Duplicate prespecified sample family")
    owners = {family: {} for family in families}
    with node_path.open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["source_family"] in families and row["event"] == "leaf":
                family, gene = row["source_family"], row["genes"]
                if gene in owners[family] or row["node_id"] != gene or "," in row["species"]:
                    raise ValueError("Invalid sampled native leaf ownership")
                owners[family][gene] = row["species"]
    rows = []
    for previous in audit["rows"]:
        family = previous["family"]
        checkpoint_path = phylo / "checkpoints" / (family + ".json")
        checkpoint = json.loads(checkpoint_path.read_text())
        raw_path = phylo / "gene_trees" / (family + ".raw.nwk")
        rooted_path = phylo / "gene_trees" / (family + ".rooted.nwk")
        if (checkpoint["status"] != "complete" or checkpoint["family_id"] != family
                or pins[raw_path]["sha256"] != checkpoint["raw_tree_sha256"]
                or pins[rooted_path]["sha256"] != checkpoint["rooted_tree_sha256"]
                or len(checkpoint["genes"]) != previous["genes"]
                or set(checkpoint["genes"]) != set(owners[family])):
            raise ValueError("Sample checkpoint no longer matches active native artifacts")
        row = validate_family(raw_path.read_text(), rooted_path.read_text(), checkpoint["genes"],
                              species_tree, owners[family])
        if row["legacy_support_map_changed"] != (not previous["support_map_preserved"]):
            raise ValueError("Independent legacy support readback differs from original audit")
        rows.append({"family": family, **row})
    if sum(row["legacy_support_map_changed"] for row in rows) != audit["families_with_changed_support_positions"]:
        raise ValueError("Original support audit total changed")
    for pin in (*checked, *sources):
        check(pin)
    write_json(output, {
        "status": "prespecified_native_support_sample_repaired_with_current_family_outputs_unchanged",
        "source_revision": revision,
        "job_id": os.environ.get("SLURM_JOB_ID"), "checked": checked, "sources": sources, "rows": rows,
        "families": len(rows), "families_with_legacy_support_changes": sum(row["legacy_support_map_changed"] for row in rows),
        "all_original_support_maps_restored": True, "all_current_family_pairs_confidence_and_groups_unchanged": True,
        "wall_s": time.monotonic() - started,
        "scope": "The same 32 checkpoint-bound families: raw-tree reroot/reconciliation only, before final satellite membership filtering. Not fresh alignment/tree/HMM inference, full-dataset output equivalence, independent biological validation or new benchmark scores.",
        "native_inference_rerun": False, "new_accuracy_scores": False,
        "defaults_promoted": False, "goal_completed": False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-revision", required=True)
    args = parser.parse_args()
    run(args.repo.resolve(), args.output.resolve(), args.source_revision)
