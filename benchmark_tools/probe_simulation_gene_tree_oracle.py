"""Fixed-candidate generating gene-tree control on the retained 70-cell panel."""

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import dendropy
from dendropy.calculate import treecompare
from Bio import SeqIO

from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.reconstruct_reconciliation_trace import apply_logged_constraints
from benchmark_tools.simulation_conditions import canonical_pairs, group_pairs, score_pairs


METHOD = "orthohmm_satellite_v2"
ARMS = ("inferred", "generating_root", "generating_rerooted")
CONDITIONS = {"baseline", "divergent", "turnover", "divergent_turnover",
              "missing20", "uneven_taxa", "taxon_count_control"}
PROTOCOL_SHA = "d6e0dbe11c8394eb51ce937762b68d1cb68c04e79b198c2addb06875db169c61"
ADMISSION_SHA = "17a1be71823b9fb7fb13082d05001d7667bb89c3366da2e3138a124e90907820"
PREPARED_SHA = "0caffab16f73019fabafe1fcfaac8c2de90960c5f1317803bf83308a6d2bf914"
SOURCE_SHA = "216d97608e6dede8960f3b06fe7036bd23da8fae5677df88583ba9369ec3c1bf"


def checked(item, inputs):
    path = Path(item.get("absolute_path") or item["path"]).resolve()
    observed = record(path)
    if any(observed[k] != item[k] for k in ("bytes", "sha256")):
        raise ValueError("Input changed: " + str(path))
    if str(path) in inputs and inputs[str(path)] != observed:
        raise ValueError("Conflicting retained identity")
    inputs[str(path)] = observed
    return path


def pinned(path, digest, inputs):
    item = record(path)
    if item["sha256"] != digest:
        raise ValueError("Pinned source or manifest changed: " + str(path))
    checked(item, inputs)
    return json.loads(path.read_text())


def frozen_module(path, inputs):
    item = record(path)
    if item["sha256"] != SOURCE_SHA:
        raise ValueError("Frozen phylogeny source changed")
    checked(item, inputs)
    spec = importlib.util.spec_from_file_location("_simulation_oracle_phylogeny", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def select_cells(admission):
    rows = [r for r in admission["records"]
            if r["method"] == METHOD and r["variant"] == "generating"]
    expected = {(c, s) for c in CONDITIONS for s in range(20261101, 20261111)}
    if len(rows) != 70 or {(r["condition"], r["seed"]) for r in rows} != expected:
        raise ValueError("Incomplete or duplicate fixed 70-cell inventory")
    return sorted(rows, key=lambda r: (r["condition"], r["seed"]))


def partition(path, universe):
    groups, seen = {}, set()
    for i, line in enumerate(path.read_text().splitlines()):
        genes = line.split()
        if not genes or len(set(genes)) != len(genes) or seen & set(genes):
            raise ValueError("Empty or duplicate candidate membership")
        if not set(genes) <= universe:
            raise ValueError("Unknown candidate gene")
        groups[f"Family{i:07d}"] = set(genes)
        seen.update(genes)
    if seen != universe:
        raise ValueError("Candidates do not partition the complete input universe")
    return groups


def tree_labels(tree):
    labels = [leaf.taxon.label for leaf in tree.leaf_node_iter()]
    if not labels or len(labels) != len(set(labels)):
        raise ValueError("Empty or duplicate gene-tree leaves")
    return set(labels)


def induced_tree(path, ancestor, parent_genes, genes):
    tree = dendropy.Tree.get(path=str(path), schema="newick",
                             preserve_underscores=True, rooting="force-rooted")
    namespace = dendropy.TaxonNamespace()
    for leaf in tree.leaf_node_iter():
        leaf.taxon = namespace.require_taxon(label=f"F{ancestor}__{leaf.taxon.label}")
    tree.taxon_namespace = namespace
    if tree_labels(tree) != set(parent_genes) or not genes <= set(parent_genes):
        raise ValueError("Generating tree differs from validated parent family")
    tree = tree.extract_tree_with_taxa_labels(sorted(genes))
    for node in tree.preorder_node_iter():
        if not node.is_leaf():
            node.label = None
    if tree_labels(tree) != genes or tree.seed_node.num_child_nodes() != 2:
        raise ValueError("Induced generating tree has incorrect leaves or root")
    return tree


def distances(inferred, generating):
    if tree_labels(inferred) != tree_labels(generating):
        raise ValueError("Tree comparison requires identical gene sets")
    result = {}
    for rooted in (True, False):
        namespace = dendropy.TaxonNamespace(sorted(tree_labels(inferred)))
        trees = []
        for original in (inferred, generating):
            tree = original.clone(depth=2)
            tree.migrate_taxon_namespace(namespace)
            if not rooted:
                tree.deroot()
            tree.is_rooted = rooted
            trees.append(tree)
        result["rooted_clade_distance" if rooted else "unrooted_split_distance"] = (
            treecompare.symmetric_difference(*trees))
    return result


def constrained_pairs(value, events, genes, active):
    pairs = set(value.ortholog_pairs)
    if not active:
        if events:
            raise ValueError("Constraints supplied without native constraint policy")
        return pairs
    high = {(a, b) for a, b, confidence in value.ortholog_pair_confidence
            if confidence == "high"}
    groups, _ = apply_logged_constraints(
        {"root_groups": list(map(set, value.root_groups)), "high_confidence_pairs": high},
        events, genes)
    index = {g: i for i, group in enumerate(groups) for g in group}
    if set(index) != genes:
        raise ValueError("Constraint application lost candidate genes")
    return {p for p in pairs if index[p[0]] == index[p[1]]}


def read_native_pairs(path, owners):
    rows = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["gene_a", "species_a", "gene_b", "species_b"]:
            raise ValueError("Unexpected native pair table")
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError("Malformed native pair row")
            if any(owners.get(row["gene_" + s]) != row["species_" + s] for s in ("a", "b")):
                raise ValueError("Native pair ownership mismatch")
            rows.append((row["gene_a"], row["gene_b"]))
    pairs, duplicates = canonical_pairs(rows, owners)
    if duplicates:
        raise ValueError("Duplicate native prediction pairs")
    return pairs


def cell(row, prepared, module, inputs):
    if row["status"] != "admitted":
        return {"label": row["label"], "condition": row["condition"], "seed": row["seed"],
                "status": "native_unavailable", "native_record": row}
    checked(row["native_report"], inputs)
    checked(row["preflight"], inputs)
    execution = json.loads(checked(row["execution"], inputs).read_text())
    status = execution["methods"][METHOD]
    if status["status"] != "process_succeeded" or status["exit_code"] != 0:
        raise ValueError("Native execution did not succeed")
    artifacts = {r["absolute_path"]: r for r in status["outputs"]}
    if len(artifacts) != len(status["outputs"]):
        raise ValueError("Duplicate native artifact path")
    directory = Path(row["native_report"]["path"]).parent / METHOD
    def artifact(relative):
        path = directory / relative
        if str(path) not in artifacts:
            raise ValueError("Artifact not in native execution inventory: " + str(path))
        return checked(artifacts[str(path)], inputs)
    dataset = next(d for d in prepared["datasets"] if
                   (d["condition"], d["seed"]) == (row["condition"], row["seed"]))
    truth_path = checked(row["truth"], inputs)
    truth = json.loads(truth_path.read_text())
    if row["truth"] != dataset["input_evidence"]["truth"]:
        raise ValueError("Dataset truth differs from frozen preparation")
    owners = {}
    for item in dataset["input_evidence"]["inputs"]:
        path = checked(item, inputs)
        for sequence in SeqIO.parse(path, "fasta"):
            if sequence.id in owners:
                raise ValueError("Duplicate FASTA ID")
            owners[sequence.id] = path.stem
    ancestors = {}
    for family, genes in truth["families"].items():
        for gene in genes:
            if gene in ancestors:
                raise ValueError("Repeated ancestral-family membership")
            ancestors[gene] = family
    if set(ancestors) != set(owners) or len(owners) != truth["extant_genes"]:
        raise ValueError("Truth family/input universe mismatch")
    candidates = partition(artifact("orthohmm_working_res/phylogeny_candidate_superfamilies.txt"), set(owners))
    candidate_index = {g: family for family, genes in candidates.items() for g in genes}
    native = read_native_pairs(checked(row["prediction_files"][0], inputs), owners)
    native_by_family = defaultdict(set)
    for pair in native:
        if candidate_index[pair[0]] != candidate_index[pair[1]]:
            raise ValueError("Native pair crosses candidate families")
        native_by_family[candidate_index[pair[0]]].add(pair)
    manifest = json.loads(artifact("orthohmm_phylogeny/provenance_manifest.json").read_text())
    if manifest["root_duplication_rule"] != "species_overlap" or manifest["pair_orthology_rule"] != "positive_paralogy":
        raise ValueError("Unexpected frozen reconciliation rules")
    species_path = artifact("orthohmm_phylogeny/species_tree.rooted.nwk")
    if record(species_path)["sha256"] != manifest["species_tree_sha256"]:
        raise ValueError("Species-tree manifest identity mismatch")
    checked(row["tree"], inputs)
    species = dendropy.Tree.get(path=str(species_path), schema="newick",
                                preserve_underscores=True, rooting="force-rooted")
    supplied = dendropy.Tree.get(path=row["tree"]["path"], schema="newick",
                                 preserve_underscores=True, rooting="force-rooted")
    if tree_labels(species) != set(owners.values()) or distances(species, supplied)["rooted_clade_distance"]:
        raise ValueError("Native species tree differs from generating supplied topology")
    policy = manifest["membership_reconciliation"]
    active = policy is not None
    constraint_file = directory / "orthohmm_working_res/phylogeny_candidate_merges.json"
    events = json.loads(artifact("orthohmm_working_res/phylogeny_candidate_merges.json").read_text()) if str(constraint_file) in artifacts else []
    if active and (policy["policy"] != "high_confidence_pair" or len(events) != policy["constraints"]):
        raise ValueError("Native constraint policy or inventory mismatch")
    constraints = defaultdict(list)
    for i, event in enumerate(events):
        families = {candidate_index[g] for key in ("source_genes", "target_genes") for g in event[key]}
        if len(families) != 1:
            raise ValueError("Constraint crosses candidate boundaries")
        constraints[next(iter(families))].append((i, event))
    parent = next(d for d in prepared["datasets"] if
                  f"{d['condition']}_{d['seed']}" == dataset["parent"])
    parent_truth_path = checked(parent["input_evidence"]["truth"], inputs)
    parent_truth = json.loads(parent_truth_path.read_text())
    tree_inventory = {r["path"]: r for r in parent_truth["inputs"]}
    native_parent = Path(parent["generating_tree"]["path"]).parent.parent
    reference, duplicates = canonical_pairs(truth["ortholog_pairs"], owners)
    if duplicates:
        raise ValueError("Duplicate reference pairs")
    reports, predicted = [], {arm: set() for arm in ARMS}
    for family, genes in candidates.items():
        row_report = {"family": family, "genes": len(genes),
                      "ancestral_families": sorted({ancestors[g] for g in genes})}
        raw_relative = f"orthohmm_phylogeny/gene_trees/{family}.raw.nwk"
        has_tree = str(directory / raw_relative) in artifacts
        inferred_pairs = native_by_family[family]
        choices = {arm: inferred_pairs for arm in ARMS}
        if not has_tree:
            counts = Counter(owners[g] for g in genes)
            if len(genes) >= 3 and len(counts) >= 2 and max(counts.values()) > 1:
                raise ValueError("Ambiguous candidate lacks a native tree")
            if set(group_pairs([genes], owners)) != inferred_pairs:
                raise ValueError("Native bypass predictions do not reproduce")
            row_report["status"] = "unambiguous_bypass"
        else:
            checkpoint = json.loads(artifact(f"orthohmm_phylogeny/checkpoints/{family}.json").read_text())
            raw_path = artifact(raw_relative)
            rooted_path = artifact(f"orthohmm_phylogeny/gene_trees/{family}.rooted.nwk")
            if (checkpoint["status"] != "complete" or checkpoint["family_id"] != family
                    or set(checkpoint["genes"]) != genes
                    or checkpoint["species_tree_sha256"] != manifest["species_tree_sha256"]
                    or checkpoint["raw_tree_sha256"] != record(raw_path)["sha256"]
                    or checkpoint["rooted_tree_sha256"] != record(rooted_path)["sha256"]):
                raise ValueError("Tree checkpoint identity mismatch")
            inferred = module.parse_gene_tree(rooted_path.read_text())
            if tree_labels(inferred) != genes:
                raise ValueError("Inferred tree differs from candidate membership")
            def reconcile(tree):
                value = module.reconcile_gene_tree(tree, species, owners, family_id=family,
                    root_duplication_rule="species_overlap", pair_orthology_rule="positive_paralogy")
                return constrained_pairs(value, constraints[family], genes, active)
            if reconcile(inferred) != inferred_pairs:
                raise ValueError("Frozen inferred-tree baseline does not reproduce native pairs")
            if len(row_report["ancestral_families"]) != 1:
                row_report["status"] = "mixed_ancestry_ineligible"
            else:
                ancestor = row_report["ancestral_families"][0]
                relative = f"G/Gene_trees/{ancestor}_prunedtree.nwk"
                generating_path = checked({**tree_inventory[relative], "absolute_path": str(native_parent / relative)}, inputs)
                generating = induced_tree(generating_path, ancestor, parent_truth["families"][ancestor], genes)
                row_report.update(status="oracle_eligible", **distances(inferred, generating))
                choices["generating_root"] = reconcile(generating)
                rerooted = module.root_gene_tree_min_duplication_loss(generating.clone(depth=2), species, owners)
                choices["generating_rerooted"] = reconcile(rerooted)
        local_owners = {g: owners[g] for g in genes}
        local_truth = {p for p in reference if set(p) <= genes}
        row_report["arms"] = {arm: score_pairs(pairs, local_truth, local_owners) for arm, pairs in choices.items()}
        row_report["changes"] = {arm: {"added": len(choices[arm] - inferred_pairs),
                                        "removed": len(inferred_pairs - choices[arm])}
                                  for arm in ARMS[1:]}
        for arm in ARMS:
            predicted[arm].update(choices[arm])
        reports.append(row_report)
    if predicted["inferred"] != native:
        raise ValueError("Whole-dataset baseline reproduction failed")
    return {"label": row["label"], "condition": row["condition"], "seed": row["seed"],
            "status": "baseline_reproduced_oracle_scored", "candidates": reports,
            "candidate_status_counts": dict(Counter(r["status"] for r in reports)),
            "arms": {arm: score_pairs(pairs, reference, owners) for arm, pairs in predicted.items()},
            "changes": {arm: {"added": len(predicted[arm] - native), "removed": len(native - predicted[arm])}
                        for arm in ARMS[1:]}}


def run(repo):
    inputs = {}
    results = repo / "benchmark_tools/results"
    protocol = record(results / "SIMULATION_GENE_TREE_ORACLE_PROTOCOL_20261004.md")
    if protocol["sha256"] != PROTOCOL_SHA:
        raise ValueError("Diagnostic protocol changed")
    checked(protocol, inputs)
    admission = pinned(repo / "benchmarks/results/simulation_tree_panel_admission_v1/results.json", ADMISSION_SHA, inputs)
    prepared = pinned(results / "simulation_tree_controls_prepared_20260917.json", PREPARED_SHA, inputs)
    module = frozen_module(repo / "benchmarks/work/publication_method_native_v2/orthohmm/phylogeny.py", inputs)
    for name, digest in (("simulation_conditions.py", "befe1c6ebb91e216b3029e4bf8ab2ef522a357b3419950a7edfba4d2f4460964"),
                         ("reconstruct_reconciliation_trace.py", "cd80d2a50e959bf9c530c3f19e1c0c1a99d53251d20791078a74f5a70b739506")):
        path = Path(__file__).with_name(name)
        if record(path)["sha256"] != digest:
            raise ValueError("Scoring or constraint helper changed")
        checked(record(path), inputs)
    checked(record(Path(__file__).with_name("prepare_ob_candidate_neighborhood.py")), inputs)
    checked(record(__file__), inputs)
    rows = []
    for item in select_cells(admission):
        rows.append(cell(item, prepared, module, inputs))
        print(json.dumps({"label": item["label"], "status": rows[-1]["status"]}), flush=True)
    summary = []
    for condition in sorted(CONDITIONS):
        cells = [r for r in rows if r["condition"] == condition and "arms" in r]
        summary.append({"condition": condition, "scored_cells": len(cells),
            "mean_metrics": {arm: {metric: sum(r["arms"][arm][metric] for r in cells) / len(cells)
                                    if cells else None for metric in ("f1", "precision", "recall")} for arm in ARMS},
            "candidate_status_counts": dict(sum((Counter(r["candidate_status_counts"]) for r in cells), Counter()))})
    for item in list(inputs.values()):
        checked(item, inputs)
    return {"status": "fixed_candidate_gene_tree_oracle_complete", "source": record(__file__),
            "protocol": protocol, "inputs": list(inputs.values()), "cells": rows, "summary": summary,
            "runtime": {"python": sys.version, "dendropy": dendropy.__version__},
            "defaults_changed": False, "independent_confirmation": False, "publication_ready": False,
            "limitations": ["Development-exposed finite-panel oracle mechanism diagnostic; no population CI or significance claim.",
                "Mixed-ancestry candidates and unambiguous bypasses retain native predictions.",
                "Generating trees/root/support/lengths contain information unavailable to real inference.",
                "A change does not identify topology-only causality or validate root-HOG ancestral-copy membership.",
                "No new HMM search, end-to-end inference or matched timing; unrelated shared-host analyses are untouched."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    report = run(args.repo.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as handle:
        handle.write(json.dumps(report, indent=2, sort_keys=True) + "\n")
