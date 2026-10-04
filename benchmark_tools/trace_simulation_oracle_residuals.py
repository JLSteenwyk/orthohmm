"""Explain the complete retained generating-root within-candidate error cohort."""

import argparse
from collections import Counter, defaultdict
from dataclasses import asdict
from itertools import combinations
import json
from pathlib import Path

import dendropy
from Bio import SeqIO

from benchmark_tools import probe_simulation_gene_tree_oracle as oracle
from benchmark_tools.readback_simulation_gene_tree_oracle import validate_counts
from benchmark_tools.reconstruct_reconciliation_trace import apply_logged_constraints
from benchmark_tools.zombi_truth import event_graph, xml_graph, ortholog_truth


ORACLE_SHA = "d2299d73238f1a4bb511320add9a62dbb427bc1ea9b4f414e74aefab20e7ad53"
READBACK_SHA = "db016a1b4c68ab39caed061122e2196aadfb3ebdc38a4d4f5ad38d9edd664dec"


def cohort(report):
    expected = {(c, s) for c in oracle.CONDITIONS for s in range(20261101, 20261111)}
    cells = report["cells"]
    if len(cells) != 70 or {(c["condition"], c["seed"]) for c in cells} != expected:
        raise ValueError("Incomplete fixed 70-cell oracle inventory")
    selected = []
    for cell in cells:
        if cell["status"] != "baseline_reproduced_oracle_scored":
            raise ValueError("Oracle cell is not complete")
        seen = set()
        for candidate in cell["candidates"]:
            if candidate["family"] in seen:
                raise ValueError("Repeated oracle candidate")
            seen.add(candidate["family"])
            score = candidate["arms"]["generating_root"]
            validate_counts(score)
            if score["fp"] or score["fn"]:
                if candidate["status"] not in {"oracle_eligible", "unambiguous_bypass"}:
                    raise ValueError("Unsupported error-cohort eligibility")
                if len(candidate["ancestral_families"]) != 1:
                    raise ValueError("Mixed ancestry cannot receive a single-family history")
                selected.append((cell, candidate))
    return selected


def history_index(graph):
    leaves, pairs = ortholog_truth(graph)
    paths, descendants = {}, {}
    def visit(node, path):
        path = (*path, node)
        event, children = graph[node]
        if event == "F":
            paths[node] = path
            members = {node}
        else:
            members = set().union(*(visit(child, path) for child in children))
        descendants[node] = members
        return members
    visit("Root_1", ())
    return leaves, pairs, paths, descendants


def common_ancestor(paths, left, right):
    if left == right or left not in paths or right not in paths:
        raise ValueError("History requires distinct extant genes")
    shared = [a for a, b in zip(paths[left], paths[right]) if a == b]
    if not shared:
        raise ValueError("History pair has no shared ancestor")
    return shared[-1]


def overlap(left, right, owners):
    return sorted({owners[g] for g in left} & {owners[g] for g in right})


def error_class(row, bypass):
    if row["truth"] == row["predicted"]:
        return None
    if row["truth"]:
        if not row["raw_predicted"]:
            return "pair_rule_exclusion"
        if not row["same_root_group"]:
            return "root_partition_filter"
        if not row["same_final_group"]:
            return "unsupported_satellite_constraint"
        raise ValueError("False negative is not explained by retained stages")
    if row["history_event"] != "D":
        raise ValueError("False positive does not have a true duplication ancestor")
    if bypass:
        return "single_copy_bypass_on_true_duplication"
    if row["pair_node"]["species_overlap_count"] or row["pair_node"]["pair_event"] == "duplication":
        raise ValueError("False positive conflicts with positive-paralogy exclusion")
    return "true_duplication_without_retained_species_overlap"


def trace_candidate(module, species, tree, graph, ancestor, genes, owners,
                    parent_owners, truth, events, active, bypass, family):
    native_names = {g: g.removeprefix(f"F{ancestor}__") for g in genes}
    if any(g == native_names[g] for g in genes):
        raise ValueError("Candidate ID does not belong to its declared ancestor")
    leaves, event_pairs, paths, descendants = history_index(graph)
    parent_names = {f"F{ancestor}__{g}" for g in leaves}
    if not parent_names <= parent_owners.keys() or not set(native_names.values()) <= leaves:
        raise ValueError("Candidate or parent ownership differs from history leaves")
    history_owners = {g: parent_owners[f"F{ancestor}__{g}"] for g in leaves}
    event_reference = {tuple(sorted((f"F{ancestor}__{a}", f"F{ancestor}__{b}")))
                       for a, b in event_pairs if f"F{ancestor}__{a}" in genes
                       and f"F{ancestor}__{b}" in genes}
    if event_reference != truth:
        raise ValueError("Original event history differs from local retained truth")
    if oracle.tree_labels(tree) != genes:
        raise ValueError("Induced tree differs from candidate membership")
    value = module.reconcile_gene_tree(tree, species, owners, family_id=family,
        root_duplication_rule="species_overlap", pair_orthology_rule="positive_paralogy")
    node_members = {n.node_id: set(n.genes) for n in value.nodes}
    children = defaultdict(list)
    for node in value.nodes:
        children[node.parent_node_id].append(node.node_id)
    node_by_pair = {}
    for node in value.nodes:
        for left, right in combinations(children[node.node_id], 2):
            for a in node_members[left]:
                for b in node_members[right]:
                    if owners[a] != owners[b]:
                        pair = tuple(sorted((a, b)))
                        if pair in node_by_pair:
                            raise ValueError("Pair has multiple induced ancestors")
                        node_by_pair[pair] = node
    all_pairs = {tuple(sorted((a, b))) for a, b in combinations(genes, 2) if owners[a] != owners[b]}
    if set(node_by_pair) != all_pairs:
        raise ValueError("Induced ancestor inventory misses cross-species pairs")
    raw = all_pairs if bypass else set(value.ortholog_pairs)
    root_index = {g: i for i, group in enumerate(value.root_groups) for g in group}
    evidence = []
    if bypass:
        predicted, final_index = all_pairs, {g: 0 for g in genes}
    elif active:
        high = {(a, b) for a, b, confidence in value.ortholog_pair_confidence if confidence == "high"}
        groups, evidence = apply_logged_constraints(
            {"root_groups": list(map(set, value.root_groups)), "high_confidence_pairs": high}, events, genes)
        final_index = {g: i for i, group in enumerate(groups) for g in group}
        predicted = {p for p in raw if final_index[p[0]] == final_index[p[1]]}
    else:
        if events:
            raise ValueError("Constraint events present without active native policy")
        final_index = {g: 0 for g in genes}
        predicted = raw
    expected = oracle.constrained_pairs(value, events, genes, active) if not bypass else all_pairs
    if predicted != expected:
        raise ValueError("Diagnostic stage filtering differs from oracle implementation")
    detached = {g: e["event_index"] for e in evidence if not e["supported"] for g in e["source_genes"]}
    pair_rows = []
    for a, b in sorted(all_pairs):
        history_lca = common_ancestor(paths, native_names[a], native_names[b])
        event, branches = graph[history_lca]
        if event not in {"S", "D"} or len(branches) != 2:
            raise ValueError("Cross-species common ancestor is not binary S/D")
        node = node_by_pair[(a, b)]
        original_members = {f"F{ancestor}__{g}" for g in descendants[history_lca]} & genes
        if original_members != node_members[node.node_id]:
            raise ValueError("Induced pair ancestor differs from original history clade")
        induced_species = overlap(*(node_members[c] for c in children[node.node_id]), owners)
        retained = [{g for g in descendants[child] if f"F{ancestor}__{g}" in genes} for child in branches]
        retained_overlap = overlap(*retained, history_owners)
        if retained_overlap != induced_species or len(induced_species) != node.species_overlap_count:
            raise ValueError("Induced topology/overlap differs from original pair ancestor")
        row = {"gene_a": a, "gene_b": b, "truth": (a, b) in truth,
               "predicted": (a, b) in predicted, "raw_predicted": (a, b) in raw,
               "history_node": history_lca, "history_event": event,
               "parent_species_overlap": overlap(*(descendants[c] for c in branches), history_owners),
               "candidate_species_overlap": retained_overlap,
               "pair_node": asdict(node), "same_root_group": root_index[a] == root_index[b],
               "same_final_group": final_index[a] == final_index[b],
               "detachment_events": {g: detached[g] for g in (a, b) if g in detached},
               "root_filter_active": active and not bypass}
        if row["truth"] != (event == "S"):
            raise ValueError("Pair truth differs from independently indexed original ancestor")
        row["error_class"] = error_class(row, bypass)
        pair_rows.append(row)
    return {"pairs": pair_rows, "constraint_evidence": evidence,
            "root_groups": [sorted(group) for group in value.root_groups],
            "reconciliation_nodes": [asdict(n) for n in value.nodes],
            "counts": {"tp": len(predicted & truth), "fp": len(predicted - truth), "fn": len(truth - predicted)},
            "error_classes": dict(Counter(r["error_class"] for r in pair_rows if r["error_class"]))}


def run(repo):
    inputs = {}
    results = repo / "benchmark_tools/results"
    detailed = oracle.pinned(results / "simulation_gene_tree_oracle_20261004.json", ORACLE_SHA, inputs)
    oracle.pinned(results / "simulation_gene_tree_oracle_readback_20261004.json", READBACK_SHA, inputs)
    admission = oracle.pinned(repo / "benchmarks/results/simulation_tree_panel_admission_v1/results.json", oracle.ADMISSION_SHA, inputs)
    prepared = oracle.pinned(results / "simulation_tree_controls_prepared_20260917.json", oracle.PREPARED_SHA, inputs)
    module = oracle.frozen_module(repo / "benchmarks/work/publication_method_native_v2/orthohmm/phylogeny.py", inputs)
    retained = {r["path"]: r for r in detailed["inputs"]}
    def retain(path):
        path = Path(path).resolve()
        if str(path) not in retained:
            raise ValueError("Selected input absent from retained oracle inventory: " + str(path))
        return oracle.checked(retained[str(path)], inputs)
    for name in ("probe_simulation_gene_tree_oracle.py",
                 "reconstruct_reconciliation_trace.py", "simulation_conditions.py", "prepare_ob_candidate_neighborhood.py"):
        retain(Path(__file__).with_name(name))
    oracle.checked(oracle.record(Path(__file__).with_name("readback_simulation_gene_tree_oracle.py")), inputs)
    oracle.checked(oracle.record(Path(__file__)), inputs)
    oracle.checked(oracle.record(Path(__file__).with_name("zombi_truth.py")), inputs)
    oracle.checked(oracle.record(results / "SIMULATION_ORACLE_RESIDUAL_PROTOCOL_20261004.md"), inputs)
    selected = cohort(detailed)
    rows = []
    by_label = {r["label"]: r for r in oracle.select_cells(admission)}
    datasets = {(d["condition"], d["seed"]): d for d in prepared["datasets"]}
    cache = {}
    for cell, candidate in selected:
        label = cell["label"]
        dataset = datasets[(cell["condition"], cell["seed"])]
        parent = next(d for d in prepared["datasets"] if f"{d['condition']}_{d['seed']}" == dataset["parent"])
        if label not in cache:
            native_row = by_label[label]
            execution = json.loads(retain(native_row["execution"]["path"]).read_text())
            status = execution["methods"][oracle.METHOD]
            if status["status"] != "process_succeeded" or status["exit_code"] != 0:
                raise ValueError("Selected native run not successful")
            artifact_refs = {r["absolute_path"]: r for r in status["outputs"]}
            directory = Path(native_row["native_report"]["path"]).parent / oracle.METHOD
            def artifact(relative):
                path = directory / relative
                if str(path) not in artifact_refs:
                    raise ValueError("Selected native artifact missing from execution inventory")
                return oracle.checked(artifact_refs[str(path)], inputs)
            def ownership(data):
                result = {}
                for item in data["input_evidence"]["inputs"]:
                    path = oracle.checked(item, inputs)
                    for sequence in SeqIO.parse(path, "fasta"):
                        if sequence.id in result:
                            raise ValueError("Repeated FASTA ID")
                        result[sequence.id] = path.stem
                return result
            owners, parent_owners = ownership(dataset), ownership(parent)
            if any(parent_owners.get(g) != s for g, s in owners.items()):
                raise ValueError("Derived ownership differs from parent")
            groups = oracle.partition(artifact("orthohmm_working_res/phylogeny_candidate_superfamilies.txt"), set(owners))
            truth = json.loads(oracle.checked(dataset["input_evidence"]["truth"], inputs).read_text())
            parent_truth = json.loads(oracle.checked(parent["input_evidence"]["truth"], inputs).read_text())
            species = dendropy.Tree.get(path=str(artifact("orthohmm_phylogeny/species_tree.rooted.nwk")),
                schema="newick", preserve_underscores=True, rooting="force-rooted")
            manifest = json.loads(artifact("orthohmm_phylogeny/provenance_manifest.json").read_text())
            policy = manifest["membership_reconciliation"]
            constraint_path = directory / "orthohmm_working_res/phylogeny_candidate_merges.json"
            events = json.loads(retain(constraint_path).read_text()) if str(constraint_path) in artifact_refs else []
            if policy is not None and (policy["policy"] != "high_confidence_pair" or policy["constraints"] != len(events)):
                raise ValueError("Native constraint policy differs from sidecar")
            cache[label] = (owners, parent_owners, groups, truth, parent_truth, species, policy is not None, events)
        owners, parent_owners, groups, truth, parent_truth, species, active, events = cache[label]
        family, ancestor = candidate["family"], candidate["ancestral_families"][0]
        genes = groups[family]
        if len(genes) != candidate["genes"] or not genes <= set(truth["families"][ancestor]):
            raise ValueError("Selected candidate size or ancestral membership changed")
        native_parent = Path(parent["generating_tree"]["path"]).parent.parent
        inventory = {r["path"]: r for r in parent_truth["inputs"]}
        def history_file(relative):
            return oracle.checked({**inventory[relative], "absolute_path": str(native_parent / relative)}, inputs)
        graph = event_graph(history_file(f"G/Gene_families/{ancestor}_events.tsv"))
        xml = xml_graph(history_file(f"G/Gene_trees/{ancestor}_rec.xml"))
        normalize = lambda g: {n: (e, tuple(sorted(c))) for n, (e, c) in g.items()}
        if normalize(graph) != normalize(xml):
            raise ValueError("Original event table and reconciled XML disagree")
        leaves, event_pairs = ortholog_truth(graph)
        names = {f"F{ancestor}__{g}" for g in leaves}
        if names != set(parent_truth["families"][ancestor]):
            raise ValueError("History leaves differ from original parent truth")
        parent_pairs = {tuple(sorted(p)) for p in parent_truth["ortholog_pairs"] if set(p) <= names}
        if parent_pairs != {tuple(sorted((f"F{ancestor}__{a}", f"F{ancestor}__{b}"))) for a, b in event_pairs}:
            raise ValueError("Original parent truth differs from event history")
        tree = oracle.induced_tree(history_file(f"G/Gene_trees/{ancestor}_prunedtree.nwk"), ancestor, names, genes)
        constraints = []
        for i, event in enumerate(events):
            sides = set(event["source_genes"]) | set(event["target_genes"])
            if sides & genes:
                if not sides <= genes:
                    raise ValueError("Constraint crosses final candidate boundaries")
                constraints.append((i, event))
        local_truth = {tuple(sorted(p)) for p in truth["ortholog_pairs"] if set(p) <= genes}
        traced = trace_candidate(module, species, tree, graph, ancestor, genes, owners,
            parent_owners, local_truth, constraints, active, candidate["status"] == "unambiguous_bypass", family)
        if traced["counts"] != {k: candidate["arms"]["generating_root"][k] for k in ("tp", "fp", "fn")}:
            raise ValueError("New event trace does not reproduce original oracle counts")
        rows.append({"cell": label, "condition": cell["condition"], "seed": cell["seed"],
                     "family": family, "ancestor": ancestor, "status": candidate["status"],
                     "genes": sorted(genes), "native_constraint_policy_active": active, **traced})
    for item in list(inputs.values()):
        oracle.checked(item, inputs)
    return {"schema_version": 1, "status": "complete_oracle_residual_history_trace",
            "selection": "All within-candidate generating-root FP/FN candidates in all 70 retained cells; post hoc diagnostic",
            "screened_cells": 70, "screened_candidates": sum(len(c["candidates"]) for c in detailed["cells"]),
            "candidates": rows, "inputs": list(inputs.values()),
            "summary": {"candidates": len(rows), "status_counts": dict(Counter(r["status"] for r in rows)),
                "pair_rows": sum(len(r["pairs"]) for r in rows),
                "counts": {k: sum(r["counts"][k] for r in rows) for k in ("tp", "fp", "fn")},
                "error_classes": dict(sum((Counter(r["error_classes"]) for r in rows), Counter()))},
            "defaults_changed": False, "independent_confirmation": False, "publication_ready": False,
            "limitations": ["Complete post hoc residual cohort, not an independent validation or population sample.",
                "Cross-candidate errors are covered by the preceding upstream trace, not this within-candidate diagnostic.",
                "Event histories are simulator truth unavailable to real inference; no new default is selected.",
                "Reuses frozen reconciliation; original XML/events independently check truth and pair ancestor correspondence.",
                "No new HMM search, native inference, or comparative runtime measurement."]}


def render(report):
    lines = ["# Generating-Root Residual Error Trace", "", report["selection"], "",
             f"Screened {report['screened_cells']} cells / {report['screened_candidates']:,} candidates.", "",
             "| Cell | Candidate | Status | Genes | TP | FP | FN | Error Classes |",
             "|---|---|---|---:|---:|---:|---:|---|"]
    for row in report["candidates"]:
        counts = row["counts"]
        classes = "; ".join(f"{k}: {v}" for k, v in sorted(row["error_classes"].items()))
        lines.append(f"| {row['cell']} | {row['family']} | {row['status']} | {len(row['genes'])} | "
                     f"{counts['tp']} | {counts['fp']} | {counts['fn']} | {classes} |")
    lines.extend(["", "## Aggregate Error Counts", "", "| Mechanism Class | Pairs |", "|---|---:|"])
    for name, count in sorted(report["summary"]["error_classes"].items()):
        lines.append(f"| {name} | {count} |")
    lines.extend(["", "All cross-species pairs in selected candidates are traced, including correctly classified pairs.",
                  "Counts are descriptive simulator-specific mechanism evidence, not independent accuracy gains or a tuning prescription.", ""])
    return "\n".join(lines)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    report = run(args.repo.resolve())
    args.output.mkdir(parents=True)
    (args.output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    (args.output / "results.md").write_text(render(report))
    print(json.dumps(report["summary"], sort_keys=True))
