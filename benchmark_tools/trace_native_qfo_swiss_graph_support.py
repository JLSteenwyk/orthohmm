"""Trace complete native SwissTrees exclusions through observed final graphs."""

import argparse
from collections import Counter, deque
import csv
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.export_native_factorial_progress import load, require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

CELLS = ("p0_c0_r0", "p0_c0_r1")
SEARCH = ("no_direct_hit", "one_direction", "both_directions")
GRAPH = ("direct_edge", "indirect_path", "disconnected")


def candidate_families(path, universe, selected, expected_groups, expected_genes):
    families, seen, count = {}, set(), 0
    with path.open() as stream:
        for index, line in enumerate(stream):
            genes = line.split()
            local = set(genes)
            require(genes and genes == sorted(local) and local <= universe and not local & seen,
                    "Invalid native candidate partition")
            seen.update(local)
            name = f"Family{index:07d}"
            if name in selected:
                families[name] = genes
            count += 1
    require(count == expected_groups and len(seen) == expected_genes == len(universe)
            and set(families) == selected, "Incomplete candidate partition or selected family inventory")
    return families


def scan_graph(path, universe, membership, expected_rows):
    require(type(expected_rows) is int and expected_rows >= 0, "Invalid expected graph size")
    previous, selected, count, nonpositive = None, [], 0, 0
    with path.open() as stream:
        for count, line in enumerate(stream, 1):
            fields = line.rstrip("\r\n").split("\t")
            require(len(fields) == 3, "Invalid graph columns")
            a, b, raw = fields
            weight = float(raw)
            pair = (a, b)
            require(a in universe and b in universe and a < b and math.isfinite(weight)
                    and (previous is None or previous < pair), "Invalid endpoint, weight or canonical graph order")
            previous = pair
            nonpositive += weight <= 0
            family = membership.get(a)
            if family is not None and membership.get(b) == family:
                selected.append(dict(source_family=family, gene_a=a, gene_b=b, weight=weight, row=count))
    require(count == expected_rows, "Final graph row count differs from admitted metrics")
    return selected, nonpositive


def path_witnesses(families, edges, queries):
    membership = {g: f for f, genes in families.items() for g in genes}
    require(len(membership) == sum(map(len, families.values())), "Overlapping selected candidate families")
    adjacency = {g: set() for g in membership}
    for row in edges:
        a, b = row["gene_a"], row["gene_b"]
        require(a in membership and b in membership and a < b and math.isfinite(row["weight"])
                and membership[a] == membership[b] == row["source_family"]
                and b not in adjacency[a], "Invalid induced graph edge")
        adjacency[a].add(b)
        adjacency[b].add(a)
    result, cache = {}, {}
    for a, b in queries:
        require(a in membership and b in membership and a != b and membership[a] == membership[b],
                "Invalid candidate graph query")
        if a not in cache:
            parents, queue = {a: None}, deque([a])
            while queue:
                node = queue.popleft()
                for neighbor in sorted(adjacency[node]):
                    if neighbor not in parents:
                        parents[neighbor] = node
                        queue.append(neighbor)
            cache[a] = parents
        parents = cache[a]
        if b not in parents:
            result[a, b] = None
        else:
            path = [b]
            while path[-1] != a:
                path.append(parents[path[-1]])
            result[a, b] = path[::-1]
    return result


def classify(search_cases, localized, families, graph_views):
    require(len(search_cases) == len(localized) and search_cases, "Incomplete case join")
    queries = [(r["gene_a"], r["gene_b"]) for r in search_cases]
    require(len(set(tuple(sorted(p)) for p in queries)) == len(queries), "Duplicate native graph query")
    membership = {g: f for f, genes in families.items() for g in genes}
    paths = [path_witnesses(families, edges, queries) for edges in graph_views]
    records = [{(r["gene_a"], r["gene_b"]): r for r in edges} for edges in graph_views]
    cases, counts = [], Counter()
    keys = ("family", "protein_a", "protein_b", "before", "after", "gene_a", "gene_b")
    for original, local in zip(search_cases, localized):
        require(all(original[k] == local[k] for k in keys)
                and original["before"] in ("TP", "FP")
                and original["after"] == ("FN" if original["before"] == "TP" else "TN"), "Changed native case identity/truth")
        a, b = original["gene_a"], original["gene_b"]
        family = local["source_family"]
        require(membership.get(a) == membership.get(b) == family, "Different localized candidate family")
        views = {}
        for cell, witnesses, direct in zip(CELLS, paths, records):
            search = original["direct_search"][cell]
            hits = search["gene_a_to_b"] + search["gene_b_to_a"]
            support = ("both_directions" if search["gene_a_to_b"] and search["gene_b_to_a"] else
                       "one_direction" if hits else "no_direct_hit")
            require(search["support"] == support, "Changed direct-search support")
            edge = direct.get(tuple(sorted((a, b))))
            path = witnesses[a, b]
            require(edge is None or (path is not None and len(path) == 2
                    and any(h["score"] == edge["weight"] for h in hits)), "Graph weight lacks stored direct search evidence")
            status = "direct_edge" if edge else "indirect_path" if path else "disconnected"
            counts[original["before"], cell, support, status] += 1
            views[cell] = dict(search_support=support, graph_support=status, direct_edge=edge,
                               shortest_path_edges=None if path is None else len(path) - 1, path=path)
        cases.append({**{k: original[k] for k in keys}, "source_family": family, "views": views})
    summary = [dict(before=label, cell=cell, search_support=support, graph_support=status,
                    pairs=counts[label, cell, support, status]) for label in ("TP", "FP") for cell in CELLS
               for support in SEARCH for status in GRAPH]
    return cases, summary


def run(search_path, search_sha, readback_path, readback_sha, plan_path, plan_sha, output):
    require(not output.exists() and not output.is_symlink(), "Require fresh graph diagnostic output")
    evidence = [record(Path(__file__).parent / "results/NATIVE_QFO_SWISS_GRAPH_SUPPORT_PROTOCOL_20261006.md")]
    search, search_ref = load(search_path, search_sha, evidence)
    verified, verified_ref = load(readback_path, readback_sha, evidence)
    require(search["schema"] == "native_qfo_swiss_direct_search_support_v1"
            and verified["schema"] == "native_qfo_swiss_search_support_code_readback_v1"
            and verified["report"] == search_ref
            and search["source"] == record(Path(__file__).with_name("trace_native_qfo_swiss_search_support.py"))
            and verified["source"] == record(Path(__file__).with_name("readback_native_qfo_swiss_search_support.py"))
            and all(r[k] is False for r in (search, verified) for k in ("new_scoring_or_admission",
                "uncertainty_admitted", "scientific_timings_admitted", "independent_confirmation", "publication_ready")),
            "Wrong direct-search source/readback/scope")
    evidence.extend([search["source"], verified["source"]])
    localized, _ = load(search["localization"]["path"], search["localization"]["sha256"], evidence)
    require(localized["schema"] == "native_qfo_swiss_reconciliation_localization_v1"
            and localized["source"] == record(Path(__file__).with_name("trace_native_qfo_swiss_reconciliation.py"))
            and all(localized[k] is False for k in ("new_scoring_or_admission", "uncertainty_admitted",
                "scientific_timings_admitted", "independent_confirmation", "publication_ready")), "Wrong localization source/scope")
    evidence.extend([localized["source"], localized["pair_ledger"]])
    check(localized["pair_ledger"])
    with open(localized["pair_ledger"]["path"], newline="") as stream:
        local = list(csv.DictReader(stream, delimiter="\t"))
    require(len(local) == search["changed_pairs"] == localized["changed_pairs_traced"] == 2023,
            "Incomplete changed native cohort")
    plan, plan_ref = load(plan_path, plan_sha, evidence)
    baseline, _ = load(plan["baseline"]["path"], plan["baseline"]["sha256"], evidence)
    require(plan["schema"] == "native_factorial_cost_plan_v1" and plan["core_commit"] == baseline["core_commit"],
            "Wrong frozen native plan/baseline")
    core = {}
    for name in ("accuracy.py", "helpers.py", "externals.py", "orthohmm.py"):
        matches = [r for r in baseline["core_sources"] if r["path"] == "orthohmm/" + name]
        require(len(matches) == 1, "Missing frozen graph source identity")
        r = matches[0]
        ref = dict(path=r["absolute_path"], bytes=r["bytes"], sha256=r["sha256"])
        check(ref)
        core[name] = ref
        evidence.append(ref)
    inventories, metrics, graphs = [], [], []
    require([v["cell"] for v in search["checkpoints"]] == list(CELLS), "Incomplete native graph views")
    for index, view in enumerate(search["checkpoints"]):
        require(view["cell"] == CELLS[index], "Changed native graph view")
        ref = view["native_output_validation"]
        inventory, _ = load(ref["path"], ref["sha256"], evidence)
        require(inventory["native_outputs_validated"] is True and inventory["cell"] == CELLS[index],
                "Invalid native graph inventory")
        pins = {r["path"]: r for r in inventory["checked_files"]}
        working = Path(view["files"]["gene_names.txt"]["path"]).parent.parent
        metric_ref = pins[str(working.parent.parent / "metrics.json")]
        metric, _ = load(metric_ref["path"], metric_ref["sha256"], evidence)
        native = metric["metadata"]["native_factorial"]
        require(metric["status"] == "complete" and native["cell"] == CELLS[index]
                and native["profile_expansion"] is native["candidate_expansion"] is False
                and native["frozen_pipeline_sha256"] == core["orthohmm.py"]["sha256"]
                and inventory["input_genes"] == metric["counts"]["genes"] == view["genes"],
                "Changed native graph settings or universe")
        path = working / "orthohmm_edges.txt"
        require(str(path) not in pins, "Graph unexpectedly inventoried; revise original-admission scope")
        graphs.append(dict(cell=CELLS[index], graph=record(path), native_output_validation=ref,
                           metrics=metric_ref, previously_inventoried=False))
        evidence.append(graphs[-1]["graph"])
        inventories.append(pins)
        metrics.append(metric)
    require(len(graphs) == 2, "Incomplete native graph views")
    names_refs = [v["files"]["gene_names.txt"] for v in search["checkpoints"]]
    require(all(r["sha256"] == names_refs[0]["sha256"] and r["bytes"] == names_refs[0]["bytes"]
                for r in names_refs), "Different native gene universe")
    for ref in names_refs:
        check(ref)
        evidence.append(ref)
    names = Path(names_refs[0]["path"]).read_text().splitlines()
    universe = set(names)
    require(names == sorted(universe) and len(names) == search["checkpoints"][0]["genes"] and all(names),
            "Invalid canonical native gene universe")
    candidate_path = Path(names_refs[0]["path"]).parent.parent / "orthohmm_edges_clustered.txt"
    candidate_ref = inventories[0][str(candidate_path)]
    check(candidate_ref)
    evidence.append(candidate_ref)
    require(candidate_ref["sha256"] == localized["reconstruction"]["reconstructed_candidate_sha256"],
            "Candidate partition differs from reconstructed phylogeny input")
    families = candidate_families(candidate_path, universe, {r["source_family"] for r in local},
                                  metrics[0]["counts"]["orthogroups"], len(names))
    membership = {g: f for f, genes in families.items() for g in genes}
    selected = []
    for graph, metric in zip(graphs, metrics):
        edges, nonpositive = scan_graph(Path(graph["graph"]["path"]), universe, membership,
                                       metric["counts"]["network_edges"])
        graph.update(rows_read=metric["counts"]["network_edges"], nonpositive_weights=nonpositive,
                     selected_induced_edges=edges)
        selected.append(edges)
    cases, summary = classify(search["cases"], local, families, selected)
    for ref in evidence:
        check(ref)
    output.mkdir(parents=True)
    table = output / "pairs.tsv"
    with table.open("x", newline="") as stream:
        fields = ["family", "protein_a", "protein_b", "before", "after", "source_family", "cell",
                  "search_support", "graph_support", "shortest_path_edges", "direct_weight", "direct_row", "path"]
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for case in cases:
            for cell in CELLS:
                view, edge = case["views"][cell], case["views"][cell]["direct_edge"]
                writer.writerow({**{k: case[k] for k in fields[:6]}, "cell": cell,
                    **{k: view[k] for k in ("search_support", "graph_support", "shortest_path_edges")},
                    "direct_weight": None if edge is None else edge["weight"], "direct_row": None if edge is None else edge["row"],
                    "path": "" if view["path"] is None else ",".join(view["path"])})
    result = dict(schema="native_qfo_swiss_graph_support_v1", search=search_ref, search_readback=verified_ref,
        plan=plan_ref, source=record(__file__), evidence=evidence, candidate_partition=candidate_ref,
        selected_families=families, graphs=graphs, cases=cases, summary=summary, pair_ledger=record(table),
        changed_pairs=len(cases), selected_genes=len(membership), graph_bytes_identical=all(
            graphs[0]["graph"][k] == graphs[1]["graph"][k] for k in ("bytes", "sha256")),
        selected_induced_graph_records_identical=selected[0] == selected[1],
        graph_original_admission_established=False, raw_search_rescanned=False, new_scoring_or_admission=False,
        uncertainty_admitted=False, scientific_timings_admitted=False, independent_confirmation=False, publication_ready=False,
        limitations=["Graphs newly observed and checked against admitted counts, not originally inventoried or continuously verified.",
            "Within-candidate unweighted paths include same-species/non-reference intermediates; no path is not global disconnection.",
            "Observed graph connectivity is not causal cluster recruitment, biological orthology or search correctness.",
            "Complete changed development-exposed cohort, not selected-default superiority or a total-HMM ablation.",
            "Direct report/readback/partition bindings reused; no new inference, clustering, scoring or uncertainty.",
            "Failed R1 native timing stays ineligible; shared-host diagnostic contention unknown and potentially tool-dependent."])
    with (output / "report.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("search", "readback", "plan"):
        parser.add_argument("--" + name, type=Path, required=True)
        parser.add_argument("--" + name + "-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.search, args.search_sha256, args.readback, args.readback_sha256, args.plan, args.plan_sha256, args.output)
    print(json.dumps({k: result[k] for k in ("changed_pairs", "selected_genes", "graph_bytes_identical",
                                           "selected_induced_graph_records_identical", "summary")}))
