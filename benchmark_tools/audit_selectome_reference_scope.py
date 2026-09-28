"""Bound taxonomic coverage without inferring the missing TreeFam mapping."""

import argparse
from bisect import bisect_right
from collections import Counter
import gzip
import json
from pathlib import Path

from Bio import Phylo

from benchmark_tools.audit_qfo_swiss_counts import read_raw
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def summarize(truth, owners, included):
    proteins = {p for _, a, b in truth for p in (a, b)}
    if set(owners) != proteins or any(not s for s in owners.values()):
        raise ValueError("Incomplete or extra protein species assignments")
    pairs = Counter()
    for (_, a, b), ortholog in truth.items():
        category = "both_inside" if owners[a] in included and owners[b] in included else "at_least_one_outside"
        pairs[category] += 1
        pairs[category + ("_ortholog" if ortholog else "_paralog")] += 1
    return dict(incident_proteins=len(proteins), proteins_by_species=dict(Counter(owners.values())),
        incident_proteins_outside=sum(owners[p] not in included for p in proteins), relations=dict(pairs))


def run(repo):
    results = repo / "benchmark_tools/results"
    audit_path = results / "qfo_treefam_counts_20260917.json"
    earlier = json.loads(audit_path.read_text())
    selected = results / "selectome_tf7ab_inventory_20260928.json"
    selected_ref = record(selected)
    inventory = json.loads(selected.read_text())
    if inventory["subtree_taxa"] != {"Euteleostomi": 9850}:
        raise ValueError("Unexpected recovered tree scope")
    base = repo / "qfo_benchmark/benchmark-webservice/reference_data/2020"
    mapping_ref, lineage_ref = record(base / "mapping.json.gz"), record(base / "lineage_tree.phyloxml")
    raw = earlier["stages"][0]["raw"]
    refs = [record(audit_path), selected_ref, mapping_ref, lineage_ref, raw]
    for ref in [*earlier["checked_inputs"], *inventory["inputs"], *refs]:
        check(ref)
    _, truth, members = read_raw(Path(raw["path"]), ["TreeFamA"])
    if len(truth) != earlier["reference_relations"] or len(members["TreeFamA"]) != earlier["relation_incident_members"]:
        raise ValueError("Retained reference counts changed")
    lineage = Phylo.read(lineage_ref["path"], "phyloxml")
    clades = [c for c in lineage.find_clades() if c.taxonomies and c.taxonomy.scientific_name == "Euteleostomi"]
    if len(clades) != 1:
        raise ValueError("Ambiguous lineage clade")
    included = {c.taxonomy.code for c in clades[0].get_terminals()}
    names = {c.taxonomy.code: c.taxonomy.scientific_name for c in lineage.get_terminals() if c.taxonomies}
    with gzip.open(mapping_ref["path"], "rt") as handle:
        mapping = json.load(handle)
    offsets, species = mapping["Goff"], mapping["species"]
    if len(offsets) != len(species) + 1 or offsets[0] != 0 or any(a >= b for a, b in zip(offsets, offsets[1:])):
        raise ValueError("Invalid native species offsets")
    owners, entries = {}, set()
    for protein in members["TreeFamA"]:
        entry = mapping["mapping"][protein]
        if type(entry) is not int or not 1 <= entry <= offsets[-1] or entry in entries:
            raise ValueError("Missing, invalid or aliased reference entry")
        entries.add(entry)
        owners[protein] = species[bisect_right(offsets, entry - 1) - 1]
        if owners[protein] not in names:
            raise ValueError("Reference species missing from retained lineage")
    summary = summarize(truth, owners, included)
    for ref in refs:
        check(ref)
    return dict(status="reference_taxonomic_scope_checked", source=record(__file__), inputs=refs,
        included_species=sorted(included), lineage_title=lineage.name,
        reference_species_names={s: names[s] for s in summary["proteins_by_species"]}, **summary,
        reference_members_without_relations_not_analyzed=earlier["reference_members_without_relations"],
        original_mapping_recovered=False, benchmark_admitted=False,
        limitations=["Uses retained QfO mapping only for existing reference IDs, not Selectome gene mapping.",
            "Outside-clade count is a necessary coverage exclusion, not an accuracy effect or reconstruction.",
            "Inside-clade membership does not establish gene, family or event coverage.",
            "The retained lineage title names QfO 2018; used as taxonomy metadata with observed reference-species coverage checked."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    value = run(args.repo.resolve())
    with args.output.open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
