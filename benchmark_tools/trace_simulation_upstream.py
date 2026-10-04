"""Trace all native simulation true pairs through retained search and graph data."""

import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path
import sys

import numpy as np
import scipy
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from Bio import SeqIO

from orthohmm.accuracy import load_accuracy_checkpoint
from benchmark_tools.audit_accuracy_checkpoint import check_arrays
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.probe_simulation_gene_tree_oracle import (
    checked, pinned, select_cells, partition, read_native_pairs, CONDITIONS, METHOD,
    ADMISSION_SHA, PREPARED_SHA)


PROTOCOL_SHA = "5ce5194d5333802e0f9e79d3903da8b0b9b7d209b449d05ef80f25c15ddd88d5"
ORACLE_SHA = "db016a1b4c68ab39caed061122e2196aadfb3ebdc38a4d4f5ad38d9edd664dec"
FIELDS = ("label", "condition", "seed", "ancestor", "gene_a", "gene_b", "hit_forward",
          "hit_reverse", "graph_direct", "graph_connected", "same_candidate", "native_predicted")
FLAGS = FIELDS[6:]


def graph(path, names):
    ids = {g: i for i, g in enumerate(names)}
    edges = set()
    with path.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) != 3 or row[0] not in ids or row[1] not in ids or row[0] == row[1]:
                raise ValueError("Invalid graph edge endpoints or schema")
            score = float(row[2])
            if not math.isfinite(score) or score <= 0:
                raise ValueError("Invalid graph edge weight")
            key = tuple(sorted((row[0], row[1])))
            if key in edges:
                raise ValueError("Duplicate undirected graph edge")
            edges.add(key)
    sources = [ids[a] for a, _ in edges]
    targets = [ids[b] for _, b in edges]
    matrix = coo_matrix((np.ones(len(edges)), (sources, targets)), shape=(len(names), len(names))).tocsr()
    count, labels = connected_components(matrix, directed=False)
    return edges, dict(zip(names, map(int, labels))), int(count)


def seed_sidecar(path, candidates, metrics, events):
    seen_candidates, seen_seeds = set(), set()
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["candidate_family", "seed_families"]:
            raise ValueError("Unexpected candidate seed sidecar")
        for row in reader:
            family = row["candidate_family"]
            seeds = row["seed_families"].split(",")
            if (family not in candidates or family in seen_candidates or not seeds
                    or len(set(seeds)) != len(seeds) or seen_seeds & set(seeds)):
                raise ValueError("Invalid candidate/seed inventory")
            for seed in seeds:
                if not seed.startswith("Seed") or len(seed) != 11 or not seed[4:].isdigit():
                    raise ValueError("Malformed seed ID")
            seen_candidates.add(family)
            seen_seeds.update(seeds)
    expected = {f"Seed{i:07d}" for i in range(metrics["seed_families"])}
    if (seen_candidates != set(candidates) or seen_seeds != expected
            or metrics["candidate_families"] != len(candidates)
            or metrics["seed_families"] - len(candidates) != metrics["merges"]
            or len(events) != metrics["merges"]):
        raise ValueError("Seed/merge count or coverage mismatch")
    index = {g: f for f, genes in candidates.items() for g in genes}
    for event in events:
        sides = [set(event[k]) for k in ("source_genes", "target_genes")]
        if (not all(sides) or sides[0] & sides[1]
                or any(len(event[k]) != len(side) for k, side in zip(("source_genes", "target_genes"), sides))
                or not set.union(*sides) <= set(index)
                or len({index[g] for g in set.union(*sides)}) != 1):
            raise ValueError("Merge evidence crosses final candidate boundary")
    return {"seed_families": len(seen_seeds), "candidate_families": len(candidates),
            "merges": len(events), "gene_to_seed_inside_merged_candidates_retained": False}


def trace_pairs(truth, owners, ancestors, hits, edges, components, candidates, native):
    if set(owners) != set(ancestors) or set(owners) != set(components) or set(owners) != set(candidates):
        raise ValueError("Stage gene universes differ")
    seen = set()
    for raw in sorted(truth):
        pair = tuple(sorted(raw))
        if (len(raw) != 2 or pair in seen or pair[0] == pair[1]
                or not set(pair) <= set(owners) or owners[pair[0]] == owners[pair[1]]
                or ancestors[pair[0]] != ancestors[pair[1]]):
            raise ValueError("Invalid or duplicated cross-species truth pair")
        seen.add(pair)
        a, b = pair
        yield {"ancestor": ancestors[a], "gene_a": a, "gene_b": b,
            "hit_forward": (a, b) in hits, "hit_reverse": (b, a) in hits,
            "graph_direct": pair in edges, "graph_connected": components[a] == components[b],
            "same_candidate": candidates[a] == candidates[b], "native_predicted": pair in native}


def summarize(rows):
    contingency, totals, families = Counter(), Counter(), {}
    for row in rows:
        flags = tuple(int(row[k]) for k in FLAGS)
        if flags[2] and not flags[3]:
            raise ValueError("Direct edge without graph connectivity")
        if flags[5] and not flags[4]:
            raise ValueError("Native pair crosses candidate boundary")
        contingency[flags] += 1
        family = families.setdefault(row["ancestor"], Counter())
        family["true_pairs"] += 1
        totals["true_pairs"] += 1
        if not flags[4]:
            totals["across_candidates"] += 1
            family["across_candidates"] += 1
            category = "connected_but_separated" if flags[3] else "different_graph_components"
            totals[category] += 1
            family[category] += 1
            orientations = flags[0] + flags[1]
            totals[f"{category}_hit_orientations_{orientations}"] += 1
            if flags[2]:
                totals["direct_graph_edge_but_separated"] += 1
                family["direct_graph_edge_but_separated"] += 1
        if not flags[5]:
            totals["native_fn"] += 1
            family["native_fn"] += 1
            if flags[4]:
                totals["native_fn_within_candidates"] += 1
                family["native_fn_within_candidates"] += 1
        else:
            totals["native_tp"] += 1
            family["native_tp"] += 1
    for key in ("across_candidates", "connected_but_separated", "different_graph_components",
                "direct_graph_edge_but_separated", "native_fn", "native_fn_within_candidates", "native_tp"):
        totals.setdefault(key, 0)
    if totals["native_fn"] != totals["across_candidates"] + totals["native_fn_within_candidates"]:
        raise ValueError("False-negative decomposition mismatch")
    return {"totals": dict(totals), "families": {f: dict(v) for f, v in sorted(families.items())},
            "contingency": [{**dict(zip(FLAGS, key)), "count": n} for key, n in sorted(contingency.items())]}


def cell(row, dataset, previous, inputs, writer):
    if row["status"] != "admitted":
        raise ValueError("Native cell unavailable; preserve failed attempt rather than omit")
    native_report = json.loads(checked(row["native_report"], inputs).read_text())
    execution = json.loads(checked(row["execution"], inputs).read_text())
    native = execution["methods"][METHOD]
    if native["status"] != "process_succeeded" or native["exit_code"] != 0:
        raise ValueError("Native cell execution did not succeed")
    artifacts = {r["absolute_path"]: r for r in native["outputs"]}
    if len(artifacts) != len(native["outputs"]):
        raise ValueError("Duplicate execution artifact")
    directory = Path(row["native_report"]["path"]).parent / METHOD
    def artifact(relative):
        path = directory / relative
        return checked(artifacts[str(path)], inputs)
    metric_path = Path(native_report["provenance"]["configured"]["methods"][METHOD]["metrics"])
    metric_ref = native_report["methods"][METHOD]["admission"]["native_validation"]["metrics"]
    metrics = json.loads(checked({**metric_ref, "absolute_path": str(metric_path)}, inputs).read_text())
    metadata = metrics["metadata"]
    expected = {"accuracy_profile": "high_sensitivity", "search_mode": "builtin", "clustering": "leiden",
                "cpm_resolution": .1, "leiden_seed": 4, "substitution_matrix": "BLOSUM62", "evalue_threshold": 1e-4}
    if (metrics["status"] != "complete" or any(metadata.get(k) != v for k, v in expected.items())
            or Path(metadata["output_directory"]).resolve() != directory.resolve()):
        raise ValueError("Unexpected native method or output-directory binding")
    truth_path = checked(row["truth"], inputs)
    if row["truth"] != dataset["input_evidence"]["truth"]:
        raise ValueError("Truth differs from frozen preparation")
    truth = json.loads(truth_path.read_text())
    owners, ancestors = {}, {}
    for item in dataset["input_evidence"]["inputs"]:
        path = checked(item, inputs)
        for sequence in SeqIO.parse(path, "fasta"):
            if sequence.id in owners:
                raise ValueError("Duplicate FASTA gene")
            owners[sequence.id] = path.stem
    for f, genes in truth["families"].items():
        for g in genes:
            if g in ancestors:
                raise ValueError("Repeated truth-family membership")
            ancestors[g] = f
    checkpoint = directory / "orthohmm_working_res/high_sensitivity_checkpoint"
    manifest = json.loads(artifact("orthohmm_working_res/high_sensitivity_checkpoint/manifest.json").read_text())
    expected_files = {"gene_names.txt", "gene_to_species.npy", "hit_queries.npy", "hit_targets.npy", "hit_scores.npy"}
    if set(manifest["files"]) != expected_files or set(p.name for p in checkpoint.iterdir()) != expected_files | {"manifest.json"}:
        raise ValueError("Unexpected numeric checkpoint inventory")
    for name, ref in manifest["files"].items():
        path = artifact("orthohmm_working_res/high_sensitivity_checkpoint/" + name)
        checked({**ref, "absolute_path": str(path)}, inputs)
    names, species, q, t, scores = load_accuracy_checkpoint(checkpoint, verify=False)
    numeric = check_arrays(names, species, q, t, scores)
    if numeric["genes"] != manifest["genes"] or numeric["hits"] != manifest["hits"] or set(names) != set(owners):
        raise ValueError("Checkpoint/input universe or count mismatch")
    codes, species_codes = {}, {}
    for g, code in zip(names, species):
        if (int(code) in codes and codes[int(code)] != owners[g]) or (owners[g] in species_codes and species_codes[owners[g]] != int(code)):
            raise ValueError("Checkpoint species mapping differs from FASTA ownership")
        codes[int(code)], species_codes[owners[g]] = owners[g], int(code)
    hits = {(names[a], names[b]) for a, b in zip(q, t)}
    if len(hits) != len(q):
        raise ValueError("Duplicate directed significant hit")
    edges, components, component_count = graph(artifact("orthohmm_working_res/orthohmm_edges.txt"), names)
    if len(edges) != metrics["counts"]["network_edges"]:
        raise ValueError("Final graph/metrics count mismatch")
    candidates = partition(artifact("orthohmm_working_res/phylogeny_candidate_superfamilies.txt"), set(names))
    candidate_index = {g: f for f, genes in candidates.items() for g in genes}
    events = json.loads(artifact("orthohmm_working_res/phylogeny_candidate_merges.json").read_text())
    seeds = seed_sidecar(artifact("orthohmm_working_res/phylogeny_candidate_seeds.tsv"),
                         candidates, metadata["phylogeny_candidate_profile"], events)
    predictions = read_native_pairs(checked(row["prediction_files"][0], inputs), owners)
    pair_rows = list(trace_pairs(truth["ortholog_pairs"], owners, ancestors, hits, edges, components, candidate_index, predictions))
    summary = summarize(pair_rows)
    totals = summary["totals"]
    old = previous["arms"]["inferred"]
    if (totals["true_pairs"] != old["eligible_true_pairs"] or totals["native_tp"] != old["tp"]
            or totals["native_fn"] != old["fn"] or len(predictions) - totals["native_tp"] != old["fp"]
            or totals["across_candidates"] != previous["true_pairs_across_candidates"]):
        raise ValueError("Trace does not reproduce original truth/native/oracle counts")
    for pair in pair_rows:
        writer.writerow({"label": row["label"], "condition": row["condition"], "seed": row["seed"],
                         **pair, **{k: int(pair[k]) for k in FLAGS}})
    return {"label": row["label"], "condition": row["condition"], "seed": row["seed"],
            "status": "native_counts_and_candidate_loss_reproduced", "graph_path": str(directory / "orthohmm_working_res/orthohmm_edges.txt"),
            "gene_names_path": str(checkpoint / "gene_names.txt"), "numeric": numeric,
            "graph_components": component_count, "graph_edges": len(edges), "seed_inventory": seeds, **summary}


def render(report):
    lines = ["# Retained Simulation Upstream Trace", "", "Exploratory finite-panel counts, not causal effects or independent confirmation.",
        "All seven conditions include all ten fixed seeds. Direct significant hits, graph edges",
        "and graph paths are not native ortholog predictions. Native TP/FN and cross-candidate",
        "counts reproduce the preceding frozen scores/oracle readback.", "",
        "| Condition | True pairs | Across candidates | Different graph components | Connected but separated | Direct edge but separated | Native FN within candidates |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in report["summary"]:
        keys = ("true_pairs", "across_candidates", "different_graph_components", "connected_but_separated", "direct_graph_edge_but_separated", "native_fn_within_candidates")
        lines.append("| " + row["condition"] + " | " + " | ".join(f"{row['totals'].get(k, 0):,}" for k in keys) + " |")
    lines.extend(["", "## Cross-Candidate Hit Evidence", "",
        "Each cell below sums true pairs over ten seeds; hit orientations refer to significant",
        "initial sequence-HMM hits, not every scored prefilter candidate.", "",
        "| Condition | Different components: 0 / 1 / 2 orientations | Connected but separated: 0 / 1 / 2 orientations |",
        "| --- | ---: | ---: |"])
    for row in report["summary"]:
        values = [" / ".join(f"{row['totals'].get(c + '_hit_orientations_' + str(i), 0):,}" for i in range(3))
                  for c in ("different_graph_components", "connected_but_separated")]
        lines.append(f"| {row['condition']} | {values[0]} | {values[1]} |")
    lines.extend(["", "## Evidence Limits", "", *["- " + s for s in report["limitations"]], ""])
    return "\n".join(lines)


def run(repo, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    inputs, results = {}, repo / "benchmark_tools/results"
    protocol = record(results / "SIMULATION_UPSTREAM_TRACE_PROTOCOL_20261004.md")
    if protocol["sha256"] != PROTOCOL_SHA:
        raise ValueError("Trace protocol changed")
    checked(protocol, inputs)
    admission = pinned(repo / "benchmarks/results/simulation_tree_panel_admission_v1/results.json", ADMISSION_SHA, inputs)
    prepared = pinned(results / "simulation_tree_controls_prepared_20260917.json", PREPARED_SHA, inputs)
    oracle = pinned(results / "simulation_gene_tree_oracle_readback_20261004.json", ORACLE_SHA, inputs)
    for relative, digest in (("orthohmm/accuracy.py", "1a35944ab7fea859143f599b1262aa111272787f6f735b9acc37d299427f2ad6"),
        ("benchmarks/work/publication_method_native_v2/orthohmm/orthohmm.py", "2afb89b9dc683e64e58208188f720e07701d4760ff09d1ac3a53c7a8075b84bb"),
        ("benchmark_tools/audit_accuracy_checkpoint.py", "b71fa9cf1d2cb954206e61179535d7e22fc39946c4de97af6fa45533f7eb9643"),
        ("benchmark_tools/probe_simulation_gene_tree_oracle.py", "3fd49552afcce4b09a58ecbbc5b1290429c0019392bc5a8b29422f9a3a55ff8c")):
        item = record(repo / relative)
        if item["sha256"] != digest:
            raise ValueError("Retained semantic/reader source changed")
        checked(item, inputs)
    checked(record(__file__), inputs)
    checked(record(Path(__file__).with_name("prepare_ob_candidate_neighborhood.py")), inputs)
    indexed = {(r["condition"], r["seed"]): r for r in oracle["cells"]}
    cells = select_cells(admission)
    if len(indexed) != 70 or len(oracle["cells"]) != 70:
        raise ValueError("Incomplete preceding oracle inventory")
    output.mkdir(parents=True)
    rows = []
    try:
        with (output / "true_pair_trace.tsv").open("x") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
            writer.writeheader()
            for row in cells:
                dataset = next(d for d in prepared["datasets"] if (d["condition"], d["seed"]) == (row["condition"], row["seed"]))
                rows.append(cell(row, dataset, indexed[row["condition"], row["seed"]], inputs, writer))
                print(json.dumps({"label": row["label"], "across_candidates": rows[-1]["totals"]["across_candidates"]}), flush=True)
        summaries = []
        for condition in sorted(CONDITIONS):
            selected = [r for r in rows if r["condition"] == condition]
            summaries.append({"condition": condition, "cells": len(selected),
                              "totals": dict(sum((Counter(r["totals"]) for r in selected), Counter()))})
        for item in list(inputs.values()):
            checked(item, inputs)
        report = {"status": "complete_retained_native_upstream_trace", "source": record(__file__), "protocol": protocol,
            "inputs": list(inputs.values()), "cells": rows, "summary": summaries, "pair_trace": record(output / "true_pair_trace.tsv"),
            "runtime": {"python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__},
            "defaults_changed": False, "inference_rerun": False, "independent_confirmation": False, "publication_ready": False,
            "limitations": ["Development-exposed descriptive trace; no population CI, significance or causal intervention.",
                "Absent initial significant hits do not distinguish prefilter rejection, score rejection, caps or ranking.",
                "Saved final graph precedes refinement/candidate expansion; paths may traverse non-reference-family genes.",
                "Connected but separated localizes to grouping/refinement/expansion as a combined boundary, not one isolated algorithm.",
                "Gene-to-seed membership inside merged candidates is not retained; no pre-expansion partition is invented.",
                "Native counts reproduce, but graph closure and direct hits are not orthology predictions.",
                "No search/inference/timing rerun, new default or full publication completion."]}
        with (output / "results.json").open("x") as handle:
            handle.write(json.dumps(report, indent=2, sort_keys=True) + "\n")
        with (output / "results.md").open("x") as handle:
            handle.write(render(report))
        return report
    except Exception as error:
        with (output / "failure.json").open("x") as handle:
            handle.write(json.dumps({"status": "failed", "error": str(error), "error_type": type(error).__name__,
                "completed_labels": [r["label"] for r in rows], "fixed_inventory": [r["label"] for r in cells],
                "source": record(__file__), "inputs_checked_so_far": list(inputs.values())}, indent=2, sort_keys=True) + "\n")
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.repo.resolve(), args.output.absolute())
