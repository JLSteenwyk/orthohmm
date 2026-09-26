"""Validate reconciliation selection and the pre-membership hierarchy table."""

import argparse
from collections import Counter, defaultdict
from itertools import groupby, zip_longest
import json
from pathlib import Path

from benchmark_tools.audit_phylogeny_structure import rows
from benchmark_tools.derive_phylogeny_events import NODE_COLUMNS
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_installed_orthobench import fasta_ids


HIERARCHY_COLUMNS = ["hog_id", "parent_hog_id", "species_tree_node", "source_family", "event", "genes"]


def needs_tree(genes, gene_species):
    counts = Counter(gene_species[gene] for gene in genes)
    return len(genes) >= 3 and len(counts) >= 2 and any(count > 1 for count in counts.values())


def expected_hierarchy(families, gene_species, node_rows, totals):
    batches = iter(groupby(node_rows, lambda row: row["source_family"]))
    for family, genes in sorted(families.items()):
        if not needs_tree(genes, gene_species):
            totals["bypassed_families"] += 1
            totals["bypass_rows"] += 1
            yield dict(zip(HIERARCHY_COLUMNS,
                           [family + ".root", "", "", family, "unambiguous", ",".join(sorted(genes))]))
            continue
        name, records = next(batches, (None, []))
        if name != family:
            raise ValueError(f"Wrong reconciliation selection/order: {family}")
        totals["reconciled_families"] += 1
        for row in records:
            if row["event"] == "leaf":
                continue
            totals["internal_node_rows"] += 1
            yield dict(zip(HIERARCHY_COLUMNS, [family + "." + row["node_id"],
                family + "." + row["parent_node_id"] if row["parent_node_id"] else "",
                row["species_tree_node"], family, row["event"], row["genes"]]))
    if next(batches, None) is not None:
        raise ValueError("Unexpected reconciled families after selection")


def audit(directory, event_report):
    evidence = json.loads(event_report.read_text())
    if evidence["status"] != "frozen_phylogeny_event_pair_semantics_verified":
        raise ValueError("Require successful event/pair semantic readback")
    prior = {item["path"]: item for item in evidence["checked_records"]}
    watched = [record(event_report)]

    def checked(path):
        item = record(path)
        if prior.get(item["path"]) != item:
            raise ValueError("Event evidence does not cover current file")
        watched.append(item)

    manifest_path = directory / "provenance_manifest.json"
    checked(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    gene_species = {}
    for item in manifest["input_proteomes"]:
        matches = [Path(path) for path, rec in prior.items()
                   if Path(path).name == item["filename"] and rec["sha256"] == item["sha256"]]
        if len(matches) != 1:
            raise ValueError("Ambiguous input identity")
        checked(matches[0])
        genes = fasta_ids(matches)
        if genes & gene_species.keys():
            raise ValueError("Duplicate gene identity")
        gene_species.update(dict.fromkeys(genes, item["taxon"]))
    roots = directory / "orthohmm_root_hogs.tsv"
    checked(roots)
    families = defaultdict(set)
    for row in rows(roots, ["root_hog", "source_family", "genes"]):
        families[row["source_family"]].update(row["genes"].split(","))
    node_path = directory / "orthohmm_reconciliation_nodes.tsv"
    checked(node_path)
    hierarchy = directory / "orthohmm_hierarchical_orthogroups.tsv"
    watched.append(record(hierarchy))
    totals = Counter()
    expected = expected_hierarchy(families, gene_species, rows(node_path, NODE_COLUMNS), totals)
    for index, (left, right) in enumerate(zip_longest(expected, rows(hierarchy, HIERARCHY_COLUMNS)), 1):
        if left != right:
            raise ValueError(f"Hierarchy differs from validated pre-membership nodes at row {index}")
        totals["hierarchical_groups"] += 1
    for key in ("reconciled_families", "bypassed_families"):
        if totals[key] != evidence["totals"][key] or totals[key] != manifest["results"][key]:
            raise ValueError(f"Selection totals differ: {key}")
    for item in watched:
        check(item)
    return dict(status="phylogeny_selection_and_hierarchy_verified", totals=dict(totals),
        checked_records=watched, reader=record(__file__), scientific_scores_admitted=False,
        hierarchy_semantics="Node hierarchy before satellite-membership filtering, not final root-HOG partition.",
        limitations=["Conditional on prior event/pair readback and frozen family membership.",
                     "Does not establish correctness of initial candidate selection or biological truth.",
                     "Hierarchy is not substituted for the benchmarked root groups or native ortholog pairs."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--events", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = audit(args.directory.resolve(), args.events.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
