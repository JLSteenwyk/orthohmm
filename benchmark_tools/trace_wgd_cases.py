"""Trace all prespecified WGD cases through retained reconciliation artifacts."""

import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.assemble_wgd_application import ADMISSIONS, artifact, unchanged
from benchmark_tools.audit_wgd_results import recalculate
from benchmark_tools.reconstruct_reconciliation_trace import reconstruct_nodes, reconstruct_bypass, apply_logged_constraints
from benchmark_tools.run_wgd_application import pinned
from benchmark_tools.snapshot_orthohmm_input_order import record
from benchmark_tools.validate_scaling_outputs import input_universe


def partition_index(groups):
    index = {}
    for name, genes in groups.items():
        if not genes or len(genes) != len(set(genes)):
            raise ValueError("Empty or duplicated group members")
        for gene in genes:
            if gene in index:
                raise ValueError("Overlapping groups")
            index[gene] = name
    return index


def loss_stage(gene, anchors, candidates, preconstraint, final):
    if any(final[gene] == final[a] for a in anchors):
        return "retained_in_anchor_group"
    if not any(candidates[gene] == candidates[a] for a in anchors):
        return "outside_anchor_candidates"
    if not any(preconstraint[gene] == preconstraint[a] for a in anchors):
        return "root_lineage_split"
    return "satellite_constraint_split"


def run(repo):
    results = repo / "benchmark_tools/results"
    admission_name, admission_sha = ADMISSIONS[1]
    admission_path = results / admission_name
    admission = pinned({"path": str(admission_path), "sha256": admission_sha})
    audit_path = results / "biological_wgd_rescore_audit_20260917.json"
    audit = pinned({"path": str(audit_path), "sha256": "1b8159c30f8ef0073223900164bdc847d089bcf578dd649dacf652336701175e"})
    report = pinned(audit["report"])
    prepared = pinned(report["input_manifest"])
    reference = pinned(prepared["reference"])
    owners, _ = input_universe({**prepared, "proteomes": 4})
    checked = []
    def check(path, expected=None):
        item = record(path)
        if expected is not None and item["sha256"] != expected:
            raise ValueError("Artifact hash mismatch: " + str(path))
        checked.append(item)
        return Path(path)
    admitted = {r["path"]: r for r in admission["verified_artifacts"] + admission["native"]["checked_files"]}
    def admitted_file(basename):
        found = [r for r in admitted.values() if Path(r["path"]).name == basename]
        if len(found) != 1:
            raise ValueError("Missing or ambiguous admitted file: " + basename)
        unchanged(found[0])
        return check(found[0]["path"], found[0]["sha256"])
    for name, digest in (("phylogeny.py", "216d97608e6dede8960f3b06fe7036bd23da8fae5677df88583ba9369ec3c1bf"),
                         ("phylogeny_pipeline.py", "44e00316546b5df78354badd2b1a6bb595b685e98b86f686f1a32543d7f15e4f")):
        check(admitted_file(name), digest)
    manifest = json.loads(admitted_file("provenance_manifest.json").read_text())
    if (manifest["root_duplication_rule"] != "species_overlap" or manifest["pair_orthology_rule"] != "positive_paralogy"
            or manifest["membership_reconciliation"]["policy"] != "high_confidence_pair"):
        raise ValueError("Unsupported reconciliation rules")
    candidate_path = admitted_file("phylogeny_candidate_superfamilies.txt")
    check(candidate_path, manifest["input_cluster_sha256"])
    candidates = {f"Family{i:07d}": line.split() for i, line in enumerate(candidate_path.read_text().splitlines())}
    candidate_index = partition_index(candidates)
    final_path = admitted_file("orthohmm_root_hogs.tsv")
    directory = final_path.parent
    check(directory / "species_tree.rooted.nwk", manifest["species_tree_sha256"])
    final, by_family = {}, defaultdict(list)
    with final_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["root_hog", "source_family", "genes"]:
            raise ValueError("Unexpected root HOG schema")
        for row in reader:
            if row["root_hog"] in final:
                raise ValueError("Repeated HOG identity")
            genes = row["genes"].split(",")
            final[row["root_hog"]] = genes
            by_family[row["source_family"]].append(genes)
    final_index = partition_index(final)
    if set(candidate_index) != set(owners) or set(final_index) != set(owners):
        raise ValueError("Incomplete native universe")
    if set(by_family) != set(candidates) or any(
            set(candidates[f]) != set().union(*map(set, by_family[f])) for f in candidates):
        raise ValueError("Final source families differ from candidates")
    keys = {tuple(e["orf_pair"]) for e in prepared["prespecified_examples"]}
    cohort = [r for r in prepared["cohort_pairs"] if tuple(r["orf_pair"]) in keys]
    selected = {candidate_index[g] for r in cohort for g in r["orf_pair"] + r.get("available_pillar_members", [])}
    nodes = defaultdict(list)
    with check(directory / "orthohmm_reconciliation_nodes.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["source_family"] in selected:
                nodes[row["source_family"]].append(row)
    events = defaultdict(list)
    for index, event in enumerate(json.loads(admitted_file("phylogeny_candidate_merges.json").read_text())):
        families = {candidate_index[g] for g in event["source_genes"] + event["target_genes"]}
        if len(families) != 1:
            raise ValueError("Constraint crosses candidate boundary")
        family = families.pop()
        if family in selected:
            events[family].append((index, event))
    preconstraint, families = {}, {}
    for family in sorted(selected):
        genes = set(candidates[family])
        if nodes[family]:
            checkpoint = json.loads(check(directory / "checkpoints" / (family + ".json")).read_text())
            if (checkpoint["status"] != "complete" or checkpoint["family_id"] != family or set(checkpoint["genes"]) != genes
                    or checkpoint["species_tree_sha256"] != manifest["species_tree_sha256"]):
                raise ValueError("Inconsistent family checkpoint")
            for suffix, key in (("raw", "raw_tree_sha256"), ("rooted", "rooted_tree_sha256"), ("reconciled", "annotated_tree_sha256")):
                check(directory / "gene_trees" / f"{family}.{suffix}.nwk", checkpoint[key])
            reconstructed = reconstruct_nodes(nodes[family], genes, owners)
        else:
            reconstructed = reconstruct_bypass(genes, owners)
        groups, constraints = apply_logged_constraints(reconstructed, events[family], genes)
        if {frozenset(g) for g in groups} != {frozenset(g) for g in by_family[family]}:
            raise ValueError("Reconstructed final membership differs")
        for i, group in enumerate(reconstructed["root_groups"]):
            preconstraint[f"{family}:pre:{i}"] = sorted(group)
        families[family] = {"candidate_members": sorted(genes), "reconciled": bool(nodes[family]),
                            "root_duplication_nodes": reconstructed["root_duplication_nodes"],
                            "propagated_split_nodes": reconstructed["propagated_split_nodes"],
                            "preconstraint_groups": [sorted(g) for g in reconstructed["root_groups"]],
                            "final_groups": [sorted(g) for g in groups], "constraints": constraints,
                            "nodes": nodes[family]}
    pre_index = partition_index(preconstraint)
    stages = {name: {tuple(r["orf_pair"]): r for r in recalculate(cohort, groups, reference, owners)}
              for name, groups in (("candidates", candidates), ("preconstraint", preconstraint), ("final", final))}
    cases = []
    for example in prepared["prespecified_examples"]:
        key = tuple(example["orf_pair"])
        row = stages["final"][key]
        expected = next(r for r in report["methods"]["orthohmm_satellite_v2"]["rows"] if tuple(r["orf_pair"]) == key)
        if row != expected:
            raise ValueError("Case outcome differs from audited application report")
        homologs = [g for g in row.get("available_pillar_members", []) if owners[g] != "Scerevisiae"] if row["reference_eligible"] else []
        cases.append({**example, "stages": {s: rows[key] for s, rows in stages.items()},
                      "homolog_destinations": [{"gene": g, "candidate": candidate_index[g], "final_group": final_index[g],
                                                "stage": loss_stage(g, key, candidate_index, pre_index, final_index)} for g in homologs]})
    for item in checked:
        unchanged(item)
    return {"status": "prespecified_case_native_membership_reconstructed", "application_audit": record(audit_path),
            "admission": record(admission_path), "sources": [record(__file__), record(Path(__file__).with_name("reconstruct_reconciliation_trace.py")),
                                                           record(Path(__file__).with_name("audit_wgd_results.py"))],
            "checked_artifacts": checked, "cases": cases, "families": families,
            "limitations": ["Retained node calls are reconstructed, not independently inferred evolutionary history.",
                            "Node-table and family checkpoint hashes are recorded retrospectively; tree hashes are cross-checked against checkpoints.",
                            "Only families incident to the six prespecified cases are reconstructed.",
                            "Stage categories follow execution order; root-lineage splits are not also counted as constraint splits.",
                            "This trace does not identify prefilter/search causes or validate the biological correctness of tree rooting.",
                            "No cross-species ancestral-copy labels are supplied by the reference pillars."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = run(args.repo.resolve())
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
