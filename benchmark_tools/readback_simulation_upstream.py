"""Independently verify native upstream pair rows and igraph connectivity."""

import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path

import igraph
import numpy as np

from benchmark_tools.readback_simulation_gene_tree_oracle import record, load


ADMISSION_SHA = "17a1be71823b9fb7fb13082d05001d7667bb89c3366da2e3138a124e90907820"
CONDITIONS = {"baseline", "divergent", "turnover", "divergent_turnover",
              "missing20", "uneven_taxa", "taxon_count_control"}
FLAGS = ("hit_forward", "hit_reverse", "graph_direct", "graph_connected", "same_candidate", "native_predicted")
FIELDS = ("label", "condition", "seed", "ancestor", "gene_a", "gene_b", *FLAGS)


def component_labels(path, names):
    ids, pairs = {g: i for i, g in enumerate(names)}, set()
    if len(ids) != len(names):
        raise ValueError("Duplicate graph vertex name")
    with path.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) != 3 or any(g not in ids for g in row[:2]):
                raise ValueError("Invalid graph endpoints")
            pairs.add(tuple(sorted(row[:2])))
    graph = igraph.Graph(n=len(names), edges=[(ids[a], ids[b]) for a, b in pairs], directed=False)
    components = graph.connected_components()
    return pairs, dict(zip(names, components.membership)), len(components)


def row_flags(row, context):
    if any(row.get(k) not in {"0", "1"} for k in FLAGS):
        raise ValueError("Nonbinary TSV flag")
    a, b = row["gene_a"], row["gene_b"]
    if a >= b or (a, b) not in context["truth"]:
        raise ValueError("Noncanonical or nonreference pair")
    if context["ancestors"][a] != row["ancestor"] or context["ancestors"][b] != row["ancestor"]:
        raise ValueError("Ancestral-family mismatch")
    actual = (int((a, b) in context["hits"]), int((b, a) in context["hits"]),
        int((a, b) in context["edges"]), int(context["components"][a] == context["components"][b]),
        int(context["candidates"][a] == context["candidates"][b]), int((a, b) in context["native"]))
    if tuple(int(row[k]) for k in FLAGS) != actual:
        raise ValueError("TSV stage flags disagree with independently read native inputs")
    return actual


def project(counter):
    total = Counter(true_pairs=sum(counter.values()))
    for state, n in counter.items():
        hit_a, hit_b, edge, connected, candidate, native = state
        total["native_tp" if native else "native_fn"] += n
        if not candidate:
            category = "connected_but_separated" if connected else "different_graph_components"
            total["across_candidates"] += n
            total[category] += n
            total[f"{category}_hit_orientations_{hit_a + hit_b}"] += n
            if edge:
                total["direct_graph_edge_but_separated"] += n
        elif not native:
            total["native_fn_within_candidates"] += n
    return total


def equal_counts(a, b):
    return {k: v for k, v in a.items() if v} == {k: v for k, v in b.items() if v}


def context_for(row, native, inputs):
    directory = Path(native["native_report"]["path"]).parent / "orthohmm_satellite_v2"
    names_path = directory / "orthohmm_working_res/high_sensitivity_checkpoint/gene_names.txt"
    edge_path = directory / "orthohmm_working_res/orthohmm_edges.txt"
    if row["graph_path"] != str(edge_path) or row["gene_names_path"] != str(names_path):
        raise ValueError("Native-directory binding mismatch")
    used = [names_path, edge_path]
    names = names_path.read_text().splitlines()
    edges, components, count = component_labels(edge_path, names)
    q_path, t_path = [names_path.parent / name for name in ("hit_queries.npy", "hit_targets.npy")]
    used += [q_path, t_path]
    q, t = [np.load(p, allow_pickle=False) for p in (q_path, t_path)]
    if q.shape != t.shape or q.ndim != 1:
        raise ValueError("Different hit dimensions")
    hits = {(names[int(a)], names[int(b)]) for a, b in zip(q, t)}
    candidate_path = edge_path.parent / "phylogeny_candidate_superfamilies.txt"
    used.append(candidate_path)
    candidates = {}
    for i, line in enumerate(candidate_path.read_text().splitlines()):
        genes = line.split()
        if not genes:
            raise ValueError("Empty candidate line")
        for g in genes:
            if g in candidates:
                raise ValueError("Repeated candidate membership")
            candidates[g] = i
    truth_ref = native["truth"]
    truth_path = Path(truth_ref.get("absolute_path") or truth_ref["path"])
    used.append(truth_path)
    truth = json.loads(truth_path.read_text())
    ancestors = {g: f for f, genes in truth["families"].items() for g in genes}
    references = {tuple(sorted(pair)) for pair in truth["ortholog_pairs"]}
    if len(references) != len(truth["ortholog_pairs"]) or set(ancestors) != set(names) or set(candidates) != set(names):
        raise ValueError("Truth/input/candidate mismatch")
    native_path = Path(native["prediction_files"][0]["path"])
    used.append(native_path)
    with native_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["gene_a", "species_a", "gene_b", "species_b"]:
            raise ValueError("Unexpected native prediction schema")
        predictions = {tuple(sorted((r["gene_a"], r["gene_b"]))) for r in reader}
    if len(edges) != row["graph_edges"] or count != row["graph_components"]:
        raise ValueError("igraph/SciPy graph count mismatch")
    for path in used:
        if str(path) not in inputs:
            raise ValueError("Independent native input not bound in trace: " + str(path))
    return {"truth": references, "ancestors": ancestors, "native": predictions, "hits": hits,
            "edges": edges, "components": components, "candidates": candidates}


def run(repo, report_path, digest):
    report, ref = load(report_path, digest)
    if report["status"] != "complete_retained_native_upstream_trace":
        raise ValueError("Trace not complete")
    inputs = {}
    for item in report["inputs"]:
        if item["path"] in inputs or record(item["path"]) != item:
            raise ValueError("Duplicate or changed trace input")
        inputs[item["path"]] = item
    if record(report["pair_trace"]["path"]) != report["pair_trace"]:
        raise ValueError("Changed all-pair TSV")
    admission, admission_ref = load(repo / "benchmarks/results/simulation_tree_panel_admission_v1/results.json", ADMISSION_SHA)
    native = {r["label"]: r for r in admission["records"] if r["method"] == "orthohmm_satellite_v2" and r["variant"] == "generating"}
    rows = {r["label"]: r for r in report["cells"]}
    expected = {f"{c}_{s}" for c in CONDITIONS for s in range(20261101, 20261111)}
    if len(report["cells"]) != 70 or set(rows) != expected or set(native) != expected:
        raise ValueError("Incomplete fixed native panel")
    contexts, pairs, states, families = {}, defaultdict(set), defaultdict(Counter), defaultdict(lambda: defaultdict(Counter))
    with Path(report["pair_trace"]["path"]).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != FIELDS:
            raise ValueError("Unexpected all-pair TSV schema")
        for row in reader:
            label = row["label"]
            if label not in rows or row["condition"] != rows[label]["condition"] or int(row["seed"]) != rows[label]["seed"]:
                raise ValueError("Unknown cell or changed condition/seed")
            if label not in contexts:
                contexts[label] = context_for(rows[label], native[label], inputs)
            flags = row_flags(row, contexts[label])
            pair = (row["gene_a"], row["gene_b"])
            if pair in pairs[label]:
                raise ValueError("Repeated truth pair in TSV")
            pairs[label].add(pair)
            states[label][flags] += 1
            families[label][row["ancestor"]][flags] += 1
    checked_cells = []
    for label, row in rows.items():
        if label not in contexts or pairs[label] != contexts[label]["truth"]:
            raise ValueError("Missing truth pair or whole dataset in TSV")
        contingency = {tuple(r[k] for k in FLAGS): r["count"] for r in row["contingency"]}
        if len(contingency) != len(row["contingency"]) or dict(states[label]) != contingency:
            raise ValueError("Contingency projection mismatch")
        totals = project(states[label])
        if not equal_counts(totals, row["totals"]):
            raise ValueError("Cell total projection mismatch")
        projected = {}
        for family, counter in families[label].items():
            counts = project(counter)
            keys = {"true_pairs", "across_candidates", "different_graph_components", "connected_but_separated",
                    "direct_graph_edge_but_separated", "native_tp", "native_fn", "native_fn_within_candidates"}
            projected[family] = {k: v for k, v in counts.items() if k in keys and v}
        if projected != row["families"]:
            raise ValueError("Family-count projection mismatch")
        checked_cells.append({"label": label, "condition": row["condition"], "seed": row["seed"], "totals": dict(totals)})
    summary = []
    for condition in sorted(CONDITIONS):
        selected = [r for r in checked_cells if r["condition"] == condition]
        total = sum((Counter(r["totals"]) for r in selected), Counter())
        original = next(r for r in report["summary"] if r["condition"] == condition)
        if original["cells"] != 10 or len(selected) != 10 or not equal_counts(total, original["totals"]):
            raise ValueError("Condition summary projection mismatch")
        summary.append({"condition": condition, "cells": len(selected), "totals": dict(total)})
    return {"status": "independent_native_rows_and_igraph_components_verified", "source": record(__file__),
            "record_helper": record(Path(__file__).with_name("readback_simulation_gene_tree_oracle.py")),
            "report": ref, "native_admission": admission_ref, "pair_trace": report["pair_trace"],
            "input_identities_rechecked": len(inputs), "pair_rows_verified": sum(map(len, pairs.values())),
            "igraph_version": igraph.__version__, "numpy_version": np.__version__,
            "cells": checked_cells, "summary": summary, "all_stage_flags_verified_against_native_inputs": True,
            "scientific_defaults_changed": False, "independent_confirmation": False, "publication_ready": False,
            "limits": "Independent native row/graph/count verification, not validation of all search decisions, native scores or evolutionary truth."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = run(args.repo.resolve(), args.report, args.report_sha256)
    with args.output.open("x") as handle:
        handle.write(json.dumps(result, indent=2, sort_keys=True) + "\n")
