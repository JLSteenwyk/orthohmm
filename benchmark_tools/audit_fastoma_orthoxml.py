"""Validate FastOMA OrthoXML identifiers, root membership and native pair scope."""

import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.fastoma_to_pairwise import input_owners, iter_pairs
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.inventory_swiss_sequences import PREPARED_SHA

NS = "{http://orthoXML.org/2011/}"
PARENTS = {"orthoXML": {None}, "species": {"orthoXML"}, "database": {"species"},
           "genes": {"database"}, "gene": {"genes"}, "taxonomy": {"orthoXML"},
           "taxon": {"taxonomy", "taxon"}, "scores": {"orthoXML"}, "scoreDef": {"scores"},
           "groups": {"orthoXML"}, "orthologGroup": {"groups", "orthologGroup", "paralogGroup"},
           "paralogGroup": {"orthologGroup", "paralogGroup"}, "geneRef": {"orthologGroup", "paralogGroup"},
           "score": {"orthologGroup", "paralogGroup"}, "property": {"orthologGroup", "paralogGroup"}}


def read_xml(path, owners):
    stack, declarations, proteins, members = [], {}, set(), {}
    species, taxa, group_ids, roots = {}, {}, set(), set()
    counts = Counter()
    expected_species = set(owners.values())
    current_species = current_root = None
    for event, elem in ET.iterparse(path, events=("start", "end")):
        if not elem.tag.startswith(NS):
            raise ValueError("Unexpected OrthoXML namespace")
        tag = elem.tag[len(NS):]
        if event == "end":
            stack.pop()
            if tag == "species":
                current_species = None
            if tag == "orthologGroup" and stack[-1].tag == NS + "groups":
                current_root = None
            if stack:
                stack[-1].remove(elem)
            elem.clear()
            continue
        parent = stack[-1].tag[len(NS):] if stack else None
        if tag not in PARENTS or parent not in PARENTS[tag]:
            raise ValueError("Unexpected OrthoXML element placement: " + tag)
        counts[tag] += 1
        if tag == "orthoXML":
            if parent is not None or elem.get("origin") != "FastOMA 0.3.5" or elem.get("version") != "0.5":
                raise ValueError("Unexpected native FastOMA OrthoXML identity")
        elif tag == "species":
            current_species = elem.get("name")
            if parent != "orthoXML" or current_species not in expected_species or current_species in species:
                raise ValueError("Unexpected or duplicate species declaration")
            species[current_species] = elem.get("taxonId")
        elif tag == "gene":
            gene, protein = elem.get("id"), elem.get("protId")
            if (parent != "genes" or current_species is None or not gene or gene in declarations
                    or protein in proteins or protein not in owners or owners[protein] != current_species):
                raise ValueError("Invalid or ambiguous OrthoXML gene declaration")
            declarations[gene] = protein
            proteins.add(protein)
        elif tag == "taxon":
            ident, name = elem.get("id"), elem.get("name")
            if parent not in ("taxonomy", "taxon") or not ident or ident in taxa or not name:
                raise ValueError("Invalid taxonomy node")
            taxa[ident] = name
        elif tag == "orthologGroup":
            ident = elem.get("id")
            if (parent not in ("groups", "orthologGroup", "paralogGroup") or not ident or ident in group_ids
                    or elem.get("taxonId") not in taxa):
                raise ValueError("Invalid ortholog group identity or taxon")
            group_ids.add(ident)
            if parent == "groups":
                # Native collect_subhogs uses the prefix before '_' in RootHOGs.tsv.
                current_root = ident.split("_")[0]
                if current_root in roots:
                    raise ValueError("Ambiguous root-HOG identifier normalization")
                roots.add(current_root)
        elif tag == "paralogGroup":
            if parent not in ("orthologGroup", "paralogGroup") or current_root is None:
                raise ValueError("Paralog group outside a root HOG")
        elif tag == "geneRef":
            protein = declarations.get(elem.get("id"))
            if (parent not in ("orthologGroup", "paralogGroup") or current_root is None
                    or protein is None or protein in members):
                raise ValueError("Unknown, repeated or misplaced gene reference")
            members[protein] = current_root
        elif tag == "score":
            if not math.isfinite(float(elem.get("value", "nan"))):
                raise ValueError("Nonfinite group score")
        stack.append(elem)
    if (counts["orthoXML"] != 1 or counts["groups"] != 1 or counts["taxonomy"] != 1
            or set(species) != expected_species or not members or not roots):
        raise ValueError("Incomplete OrthoXML species/group scope")
    if any(taxa.get(taxon) != name for name, taxon in species.items()):
        raise ValueError("Species taxonomy identity mismatch")
    covered_roots = set(members.values())
    if covered_roots != roots:
        raise ValueError("Empty root HOG")
    return {"input_proteins": len(owners), "species": len(species), "declared_proteins": len(proteins),
            "input_not_declared": len(owners) - len(proteins), "hog_referenced_proteins": len(members),
            "declared_not_hog_referenced": len(proteins) - len(members), "root_hogs": len(roots),
            "element_counts": dict(counts)}, members


def check_root_table(path, members):
    seen = set()
    with path.open() as stream:
        rows = csv.reader(stream, delimiter="\t")
        if next(rows, None) != ["RootHOG", "Protein", "OMAmerRootHOG"]:
            raise ValueError("Unexpected RootHOG table header")
        for row in rows:
            if len(row) != 3 or not all(row):
                raise ValueError("Malformed RootHOG row")
            root, protein, _ = row
            if protein in seen or members.get(protein) != root:
                raise ValueError("RootHOG table differs from XML membership")
            seen.add(protein)
    if seen != set(members):
        raise ValueError("Incomplete RootHOG table")
    return len(seen)


def check_pair_scope(path, owners, members):
    count, endpoints = 0, set()
    for a, b in iter_pairs(path, owners):
        if a not in members or b not in members or members[a] != members[b]:
            raise ValueError("Native pair outside its XML root-HOG scope")
        count += 1
        endpoints.update((a, b))
    if not count:
        raise ValueError("Empty native pair file")
    return {"native_pair_rows": count, "pair_endpoint_proteins": len(endpoints),
            "hog_members_without_pairs": len(members) - len(endpoints)}


def audit_original(root, destination):
    if destination.exists():
        raise FileExistsError(destination)
    prepared_path = root / "benchmark_tools/results/qfo_factorial_prepared_20260917.json"
    prepared = read_frozen(prepared_path, PREPARED_SHA)
    directory = root / "qfo_benchmark/results/fastoma/output"
    xml, table, pairs = [directory / name for name in ("FastOMA_HOGs.orthoxml", "RootHOGs.tsv", "orthologs.tsv.gz")]
    checked = [record(prepared_path), *prepared["input_fastas"], record(xml), record(table), record(pairs),
               record(__file__), record(Path(__file__).with_name("fastoma_to_pairwise.py"))]
    for item in checked:
        check(item)
    owners = input_owners([Path(r["path"]) for r in prepared["input_fastas"]])
    content, members = read_xml(xml, owners)
    content["root_table_members"] = check_root_table(table, members)
    content.update(check_pair_scope(pairs, owners, members))
    for item in checked:
        check(item)
    report = {"status": "historical_fastoma_xml_membership_and_pair_scope_verified",
              "checked_records": checked, "content": content, "accuracy_evaluated": False,
              "corrected_release_result": False, "publication_ready": False,
              "limitations": ["Structural/identifier consistency, not biological orthology correctness or complete workflow admission.",
                  "Proteins omitted from declarations or HOGs are reported, not silently treated as assigned.",
                  "Root-HOG co-membership is necessary but not sufficient for an orthology relation.",
                  "No pair uniqueness assertion here; distinct-pair audit is separate.",
                  "Reference-relative coverage and causes of filtering require separate analysis."]}
    with destination.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit_original(args.root.resolve(), args.output.resolve())
