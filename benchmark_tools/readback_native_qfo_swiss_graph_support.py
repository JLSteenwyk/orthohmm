"""Independent CSV scan and Floyd-Warshall check of native graph diagnostics."""

import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.readback_native_qfo_swiss_transitions import record, require

CELLS = ("p0_c0_r0", "p0_c0_r1")
SUPPORT = ("no_direct_hit", "one_direction", "both_directions")
STATUS = ("direct_edge", "indirect_path", "disconnected")
HEADER = ("family", "protein_a", "protein_b", "before", "after", "source_family", "cell", "search_support",
          "graph_support", "shortest_path_edges", "direct_weight", "direct_row", "path")


def distances(genes, edges):
    require(0 < len(genes) <= 2048 and len(set(genes)) == len(genes), "Invalid or excessive distance-matrix size")
    ids = {g: i for i, g in enumerate(genes)}
    matrix = np.full((len(genes), len(genes)), len(genes) + 1, dtype=np.int32)
    np.fill_diagonal(matrix, 0)
    for a, b in edges:
        require(a in ids and b in ids and a != b, "Invalid induced distance endpoint")
        matrix[ids[a], ids[b]] = matrix[ids[b], ids[a]] = 1
    for k in range(len(genes)):
        matrix = np.minimum(matrix, matrix[:, k:k + 1] + matrix[k:k + 1, :])
    return ids, matrix


def csv_edges(path, universe, families, expected):
    membership = {g: f for f, genes in families.items() for g in genes}
    require(len(membership) == sum(map(len, families.values())), "Duplicate selected family members")
    edges, previous, count, nonpositive = [], None, 0, 0
    with path.open(newline="") as stream:
        for count, fields in enumerate(csv.reader(stream, delimiter="\t"), 1):
            require(len(fields) == 3, "Invalid graph row width")
            a, b = fields[:2]
            value = float(fields[2])
            require(a in universe and b in universe and a < b and math.isfinite(value)
                    and (previous is None or (a, b) > previous), "Invalid canonical graph record")
            previous = a, b
            nonpositive += value <= 0
            if a in membership and b in membership and membership[a] == membership[b]:
                edges.append(dict(source_family=membership[a], gene_a=a, gene_b=b, row=count, weight=value))
    require(type(expected) is int and count == expected, "Wrong complete graph size")
    return edges, nonpositive


def validate_view(view, a, b, search, records, ids, matrix):
    edge = records.get(tuple(sorted((a, b))))
    value = int(matrix[ids[a], ids[b]])
    distance = value if value < len(ids) + 1 else None
    status = "direct_edge" if edge is not None else "indirect_path" if distance is not None else "disconnected"
    hits = search["gene_a_to_b"] + search["gene_b_to_a"]
    support = ("both_directions" if search["gene_a_to_b"] and search["gene_b_to_a"] else
               "one_direction" if hits else "no_direct_hit")
    require(view["direct_edge"] == edge and view["shortest_path_edges"] == distance and view["graph_support"] == status
            and view["search_support"] == search["support"] == support
            and (edge is None or any(h["score"] == edge["weight"] for h in hits)), "Wrong native graph case values")
    path = view["path"]
    if distance is None:
        require(path is None, "Disconnected case has a witness")
    else:
        require(type(path) is list and len(path) == distance + 1 and path[0] == a and path[-1] == b
                and len(set(path)) == len(path) and all(g in ids for g in path), "Invalid shortest-path witness")
        require(all(tuple(sorted((x, y))) in records for x, y in zip(path, path[1:])), "Witness uses absent graph edge")
    return support, status


def verify(report_path, report_sha):
    report_ref = record(report_path)
    require(report_ref["sha256"] == report_sha, "Changed graph diagnostic report")
    report = json.loads(Path(report_path).read_text())
    require(report["schema"] == "native_qfo_swiss_graph_support_v1"
            and report["source"] == record(Path(__file__).with_name("trace_native_qfo_swiss_graph_support.py"))
            and all(report[k] is False for k in ("graph_original_admission_established", "raw_search_rescanned",
                "new_scoring_or_admission", "uncertainty_admitted", "scientific_timings_admitted",
                "independent_confirmation", "publication_ready")), "Wrong graph source/scope")
    refs = [report_ref, report["source"], *report["evidence"], report["pair_ledger"]]
    for ref in refs:
        require(record(ref["path"]) == ref, "Changed bound graph evidence")
    search = json.loads(Path(report["search"]["path"]).read_text())
    readback = json.loads(Path(report["search_readback"]["path"]).read_text())
    require(readback["report"] == report["search"] and report["search"] in refs and report["search_readback"] in refs,
            "Wrong original search readback")
    localized_ref = search["localization"]
    require(record(localized_ref["path"]) == localized_ref, "Changed original localization")
    localized = json.loads(Path(localized_ref["path"]).read_text())
    require(record(localized["pair_ledger"]["path"]) == localized["pair_ledger"], "Changed localized pair ledger")
    with open(localized["pair_ledger"]["path"], newline="") as stream:
        local = list(csv.DictReader(stream, delimiter="\t"))
    families = {}
    selected = {r["source_family"] for r in local}
    names_ref = search["checkpoints"][0]["files"]["gene_names.txt"]
    require(record(names_ref["path"]) == names_ref, "Changed original gene universe")
    names = Path(names_ref["path"]).read_text().splitlines()
    universe = set(names)
    require(names == sorted(universe) and all(names) and len(names) == search["checkpoints"][0]["genes"],
            "Invalid complete gene universe")
    count, seen = 0, set()
    candidate_ref = report["candidate_partition"]
    require(candidate_ref in refs and candidate_ref["sha256"] == localized["reconstruction"]["reconstructed_candidate_sha256"],
            "Wrong original candidate fingerprint")
    with open(report["candidate_partition"]["path"]) as stream:
        for count, line in enumerate(stream, 1):
            genes = line.split()
            require(genes == sorted(set(genes)) and genes and set(genes) <= universe and not seen.intersection(genes),
                    "Invalid candidate partition members")
            seen.update(genes)
            family = "Family%07d" % (count - 1)
            if family in selected:
                families[family] = genes
    require(families == report["selected_families"] and set(families) == selected
            and seen == universe and report["selected_genes"] == sum(map(len, families.values())), "Wrong induced family universe")
    require([g["cell"] for g in report["graphs"]] == list(CELLS), "Wrong graph views")
    actual_edges, matrices = [], []
    for graph, original_view in zip(report["graphs"], search["checkpoints"]):
        working = Path(original_view["files"]["gene_names.txt"]["path"]).parent.parent
        require(all(graph[k] in refs for k in ("graph", "metrics", "native_output_validation"))
                and graph["graph"]["path"] == str(working / "orthohmm_edges.txt")
                and graph["native_output_validation"] == original_view["native_output_validation"],
                "Wrong native graph view binding")
        metric = json.loads(Path(graph["metrics"]["path"]).read_text())
        inventory = json.loads(Path(graph["native_output_validation"]["path"]).read_text())
        pins = {r["path"]: r for r in inventory["checked_files"]}
        require(graph["previously_inventoried"] is False and inventory["native_outputs_validated"] is True
                and inventory["cell"] == graph["cell"] and graph["graph"]["path"] not in {r["path"] for r in inventory["checked_files"]}
                and graph["metrics"] == pins[str(working.parent.parent / "metrics.json")]
                and metric["metadata"]["native_factorial"]["cell"] == graph["cell"]
                and count == localized["reconstruction"]["source_families"]
                and len(names) == metric["counts"]["genes"], "Wrong original graph/partition scope")
        if graph["cell"] == CELLS[0]:
            require(candidate_ref == pins[str(working / "orthohmm_edges_clustered.txt")], "Wrong original candidate admission")
        edges, nonpositive = csv_edges(Path(graph["graph"]["path"]), universe, families, metric["counts"]["network_edges"])
        require(edges == graph["selected_induced_edges"] and nonpositive == graph["nonpositive_weights"]
                and graph["rows_read"] == metric["counts"]["network_edges"], "Induced graph evidence differs")
        by_family = {f: [(r["gene_a"], r["gene_b"]) for r in edges if r["source_family"] == f] for f in families}
        matrices.append({f: distances(genes, by_family[f]) for f, genes in families.items()})
        actual_edges.append({(r["gene_a"], r["gene_b"]): r for r in edges})
    cases = report["cases"]
    require(len(cases) == len(local) == len(search["cases"]) == report["changed_pairs"] == 2023,
            "Incomplete changed graph cases")
    summary, table = Counter(), []
    for case, origin, row in zip(cases, search["cases"], local):
        require(all(case[k] == origin[k] == row[k] for k in ("family", "protein_a", "protein_b", "before", "after", "gene_a", "gene_b"))
                and case["source_family"] == row["source_family"]
                and set(case["views"]) == set(CELLS), "Changed graph-case identities")
        for cell, records, view_matrices in zip(CELLS, actual_edges, matrices):
            view = case["views"][cell]
            support, status = validate_view(view, case["gene_a"], case["gene_b"], origin["direct_search"][cell],
                                             records, *view_matrices[case["source_family"]])
            summary[case["before"], cell, support, status] += 1
            edge = view["direct_edge"]
            values = {**{k: case[k] for k in HEADER[:6]}, "cell": cell,
                "search_support": support, "graph_support": status, "shortest_path_edges": view["shortest_path_edges"],
                "direct_weight": None if edge is None else edge["weight"], "direct_row": None if edge is None else edge["row"],
                "path": "" if view["path"] is None else ",".join(view["path"])}
            table.append({k: "" if v is None else str(v) for k, v in values.items()})
    wanted = [dict(before=label, cell=cell, search_support=support, graph_support=status, pairs=summary[label, cell, support, status])
              for label in ("TP", "FP") for cell in CELLS for support in SUPPORT for status in STATUS]
    require(report["summary"] == wanted, "Wrong graph-support cross-tabulation")
    with open(report["pair_ledger"]["path"], newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        require(tuple(reader.fieldnames or ()) == HEADER and list(reader) == table, "Wrong complete graph-case TSV")
    graphs = report["graphs"]
    require(report["graph_bytes_identical"] is all(graphs[0]["graph"][k] == graphs[1]["graph"][k] for k in ("bytes", "sha256"))
            and report["selected_induced_graph_records_identical"] is (graphs[0]["selected_induced_edges"] == graphs[1]["selected_induced_edges"]),
            "Wrong graph agreement scope")
    for ref in refs:
        require(record(ref["path"]) == ref, "Graph evidence changed during independent scan")
    return dict(schema="native_qfo_swiss_graph_support_floyd_readback_v1", report=report_ref, source=record(__file__),
        checked_inputs=refs, numpy_version=np.__version__, changed_pairs_checked=len(cases), table_rows_checked=len(table),
        selected_families_checked=len(families), selected_genes_checked=report["selected_genes"],
        complete_graph_rows_checked=sum(g["rows_read"] for g in graphs),
        induced_edges_checked=sum(len(g["selected_induced_edges"]) for g in graphs), summary=wanted,
        graph_original_admission_established=False, raw_search_rescanned=False, new_scoring_or_admission=False,
        uncertainty_admitted=False, scientific_timings_admitted=False, independent_confirmation=False, publication_ready=False,
        limitations=["Independent CSV scan and Floyd-Warshall distances/witness check, not native/biological replication.",
                     "Graphs remain newly observed current-byte evidence, not original graph admission or continuous integrity.",
                     "Unweighted paths restricted to original candidates; no induced path is not global disconnection.",
                     "No inference/Leiden/scoring/uncertainty replay or repair of failed R1 timing."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), "Require fresh readback output")
    result = verify(args.report, args.report_sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: result[k] for k in ("changed_pairs_checked", "table_rows_checked", "selected_families_checked",
        "selected_genes_checked", "complete_graph_rows_checked", "induced_edges_checked")}))
