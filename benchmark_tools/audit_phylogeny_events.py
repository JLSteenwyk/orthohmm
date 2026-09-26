"""Recompute frozen phylogeny events, root groups and complete native pair sets."""

import argparse
from collections import Counter, defaultdict
from itertools import combinations, groupby
import json
from pathlib import Path

import dendropy

from benchmark_tools.audit_phylogeny_structure import rows
from benchmark_tools.derive_phylogeny_events import NODE_COLUMNS, derive, constrain
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_installed_orthobench import fasta_ids


def read_tree(path):
    return dendropy.Tree.get(path=str(path), schema="newick", preserve_underscores=True,
                            rooting="force-rooted")


def audit(directory, structural_report, constraints_path=None):
    evidence = json.loads(structural_report.read_text())
    if evidence["status"] != "phylogeny_structure_verified":
        raise ValueError("Require structural readback")
    prior = {item["path"]: item for item in evidence["checked_records"]}
    watched = [record(structural_report), *prior.values()]
    for item in watched:
        check(item)
    manifest_path = directory / "provenance_manifest.json"
    if prior.get(str(manifest_path)) != record(manifest_path):
        raise ValueError("Wrong structural evidence directory")
    manifest = json.loads(manifest_path.read_text())
    if (manifest["root_duplication_rule"] != "species_overlap"
            or manifest["pair_orthology_rule"] != "positive_paralogy"):
        raise ValueError("Oracle implements only the frozen publication rules")
    membership = manifest["membership_reconciliation"]
    if bool(membership) != bool(constraints_path):
        raise ValueError("Require the recorded membership constraints exactly when applied")
    constraints = []
    if constraints_path:
        watched.append(record(constraints_path))
        constraints = json.loads(constraints_path.read_text())
        if len(constraints) != membership["constraints"]:
            raise ValueError("Constraint count mismatch")
    gene_species = {}
    for item in manifest["input_proteomes"]:
        paths = [Path(p) for p, rec in prior.items()
                 if Path(p).name == item["filename"] and rec["sha256"] == item["sha256"]]
        if len(paths) != 1:
            raise ValueError("Ambiguous input identity")
        gene_species.update(dict.fromkeys(fasta_ids(paths), item["taxon"]))
    groups, gene_family = defaultdict(list), {}
    for row in rows(directory / "orthohmm_root_hogs.tsv", ["root_hog", "source_family", "genes"]):
        members = frozenset(row["genes"].split(","))
        groups[row["source_family"]].append(members)
        gene_family.update(dict.fromkeys(members, row["source_family"]))
    by_family = defaultdict(list)
    for constraint in constraints:
        genes = set(constraint["source_genes"]) | set(constraint["target_genes"])
        if not genes or not genes <= gene_family.keys():
            raise ValueError("Unknown constraint gene")
        families = {gene_family[gene] for gene in genes}
        if len(families) != 1:
            raise ValueError("Cross-family constraint")
        by_family[next(iter(families))].append(constraint)
    actual_pairs = defaultdict(dict)
    pair_columns = ["gene_a", "species_a", "gene_b", "species_b", "confidence"]
    for row in rows(directory / "orthohmm_pairwise_orthologs_confidence.tsv", pair_columns):
        a, b = row["gene_a"], row["gene_b"]
        actual_pairs[gene_family[a]][a, b] = row["confidence"]
    nodes_path = directory / "orthohmm_reconciliation_nodes.tsv"
    watched.append(record(nodes_path))
    node_groups = iter(groupby(rows(nodes_path, NODE_COLUMNS), lambda row: row["source_family"]))
    species_tree = read_tree(directory / "species_tree.rooted.nwk")
    totals, membership_totals = Counter(), Counter()
    for family, actual_groups in sorted(groups.items()):
        genes = set().union(*actual_groups)
        path = directory / "gene_trees" / f"{family}.rooted.nwk"
        if str(path) in prior:
            expected_rows, original_groups, pairs = derive(read_tree(path), species_tree, gene_species, family)
            observed_family, observed_rows = next(node_groups, (None, []))
            if observed_family != family or list(observed_rows) != expected_rows:
                raise ValueError(f"Reconciliation node semantics differ: {family}")
            annotated = read_tree(directory / "gene_trees" / f"{family}.reconciled.nwk")
            observed = {}
            for node in annotated.postorder_node_iter():
                clade = frozenset(leaf.taxon.label for leaf in node.leaf_iter())
                if clade in observed:
                    raise ValueError("Repeated annotated-tree clade")
                observed[clade] = node.taxon.label if node.is_leaf() else node.label
            expected = {}
            for row in expected_rows:
                label = row["node_id"]
                if row["event"] != "leaf":
                    code = {"speciation": "S", "duplication": "D", "uncertain": "U"}[row["pair_event"]]
                    label += f"|{code}@{row['species_tree_node']}"
                expected[frozenset(row["genes"].split(","))] = label
            if observed != expected:
                raise ValueError(f"Annotated tree topology/events differ: {family}")
            totals["reconciled_families"] += 1
            totals["nodes"] += len(expected_rows)
            totals["duplications"] += sum(row["event"] == "duplication" for row in expected_rows)
            totals["speciations"] += sum(row["event"] == "speciation" for row in expected_rows)
            totals["uncertain_events"] += sum(row["pair_event"] == "uncertain" for row in expected_rows)
        else:
            original_groups = [frozenset(genes)]
            pairs = {tuple(sorted((a, b))): "high" for a, b in combinations(genes, 2)
                     if gene_species[a] != gene_species[b]}
            totals["bypassed_families"] += 1
        refined, expected_pairs, counts = constrain(original_groups, pairs, by_family[family], bool(membership))
        if set(refined) != set(actual_groups):
            raise ValueError(f"Final root-group semantics differ: {family}")
        if expected_pairs != actual_pairs.pop(family, {}):
            raise ValueError(f"Native pair completeness/confidence differs: {family}")
        membership_totals.update(counts)
        totals["ortholog_pairs"] += len(expected_pairs)
        totals["root_hogs"] += len(refined)
    if next(node_groups, None) is not None or actual_pairs:
        raise ValueError("Unexpected trailing node or pair families")
    summary = manifest["results"]
    for key in ("reconciled_families", "bypassed_families", "duplications", "speciations",
                "uncertain_events", "ortholog_pairs", "root_hogs"):
        if totals[key] != summary[key]:
            raise ValueError(f"Summary semantics mismatch: {key}")
    if membership and {**membership_totals, "policy": "high_confidence_pair"} != membership:
        raise ValueError("Satellite-membership audit differs")
    for item in watched:
        check(item)
    return dict(status="frozen_phylogeny_event_pair_semantics_verified", totals=dict(totals),
        membership_totals=dict(membership_totals), checked_records=watched,
        sources=[record(__file__), record(Path(__file__).with_name("derive_phylogeny_events.py"))],
        dendropy_version=dendropy.__version__, scientific_scores_admitted=False,
        limitations=["Independent rule implementation, but the same DendroPy parser family is used.",
                     "Conditional on saved rooted trees and supplied constraint sidecar; no biological truth claim.",
                     "Does not prove optimal rooting, HMM search correctness or hierarchy-table correctness.",
                     "Constraint historical provenance must be pinned separately by caller."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--structure", type=Path, required=True)
    parser.add_argument("--constraints", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = audit(args.directory.resolve(), args.structure.resolve(), args.constraints)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
