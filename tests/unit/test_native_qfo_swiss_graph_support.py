"""Graph/path oracles and source-bound synthetic diagnostics, not admission."""

import copy
import csv
import itertools
import json
from pathlib import Path

import pytest

from benchmark_tools import trace_native_qfo_swiss_graph_support as trace
from benchmark_tools import readback_native_qfo_swiss_graph_support as reader


def edge(a, b, row=1, weight=1.0):
    return dict(source_family="Family0000000", gene_a=a, gene_b=b, row=row, weight=weight)


def test_path_oracle_includes_nonreference_and_same_species_nodes():
    families = {"Family0000000": ["a", "b", "c", "d", "e"]}
    edges = [edge("a", "b"), edge("a", "c", 2), edge("b", "d", 3), edge("c", "d", 4)]
    paths = trace.path_witnesses(families, edges, [("a", "d"), ("d", "a"), ("a", "e"), ("a", "b")])
    assert paths == {("a", "d"): ["a", "b", "d"], ("d", "a"): ["d", "b", "a"], ("a", "e"): None, ("a", "b"): ["a", "b"]}
    reference, species = {"a", "d"}, {"a": "one", "b": "one", "c": "two", "d": "three", "e": "four"}
    assert paths["a", "d"][1] not in reference and species[paths["a", "d"][1]] == species["a"]
    ids, matrix = reader.distances(families["Family0000000"], [(r["gene_a"], r["gene_b"]) for r in edges])
    assert matrix[ids["a"], ids["d"]] == 2 and matrix[ids["a"], ids["e"]] == 6


@pytest.mark.parametrize("fault", ("duplicate", "unknown", "self", "family", "nan", "query", "overlap"))
def test_invalid_induced_graph_refuses(fault):
    families, edges, queries = {"Family0000000": ["a", "b"]}, [edge("a", "b")], [("a", "b")]
    if fault == "duplicate":
        edges.append(copy.deepcopy(edges[0]))
    elif fault == "unknown":
        edges[0]["gene_b"] = "z"
    elif fault == "self":
        edges[0]["gene_b"] = "a"
    elif fault == "family":
        edges[0]["source_family"] = "wrong"
    elif fault == "nan":
        edges[0]["weight"] = float("nan")
    elif fault == "query":
        queries = [("a", "a")]
    else:
        families["Family0000001"] = ["b"]
    with pytest.raises(ValueError):
        trace.path_witnesses(families, edges, queries)


@pytest.mark.parametrize("function", (trace.scan_graph, reader.csv_edges))
@pytest.mark.parametrize("content,count", (("a\tb\t1\n", 2), ("a\tb\tnan\n", 1), ("b\ta\t1\n", 1),
    ("a\ta\t1\n", 1), ("a\tz\t1\n", 1), ("a\tb\t1\t2\n", 1),
    ("a\tb\t1\na\tb\t2\n", 2), ("b\tc\t1\na\tb\t2\n", 2)))
def test_complete_graph_validation(tmp_path, function, content, count):
    path = tmp_path / "edges.tsv"
    path.write_text(content)
    families = {"Family0000000": ["a", "b", "c"]}
    selection = {g: f for f, genes in families.items() for g in genes} if function is trace.scan_graph else families
    with pytest.raises(ValueError):
        function(path, {"a", "b", "c"}, selection, count)


def test_graph_scanners_keep_nonpositive_weights_and_induced_members(tmp_path):
    path = tmp_path / "edges.tsv"
    path.write_text("a\tb\t-1\na\tc\t0\na\td\t2\n")
    family = {"Family0000000": ["a", "b", "c"]}
    membership = {g: f for f, genes in family.items() for g in genes}
    left = trace.scan_graph(path, {"a", "b", "c", "d"}, membership, 3)
    right = reader.csv_edges(path, {"a", "b", "c", "d"}, family, 3)
    assert left == right == ([edge("a", "b", weight=-1.0), edge("a", "c", row=2, weight=0.0)], 2)


@pytest.mark.parametrize("genes,edges", (([], []), (["a", "a"], []), (["a", "b"], [("a", "z")]),
                                         ([str(i) for i in range(2049)], [])))
def test_floyd_matrix_guards_before_allocation(genes, edges):
    with pytest.raises(ValueError):
        reader.distances(genes, edges)


def save(path, value):
    path.write_text(json.dumps(value))
    return trace.record(path)


@pytest.fixture
def synthetic(tmp_path):
    # Entire fixtures are generated; no original inference or admission is run.
    flags = {k: False for k in ("new_scoring_or_admission", "uncertainty_admitted", "scientific_timings_admitted",
                                "independent_confirmation", "publication_ready")}
    names = ["g%03d" % i for i in range(67)]
    pairs = list(itertools.combinations(names, 2))[:2023]
    core_dir = tmp_path / "synthetic_core/orthohmm"
    core_dir.mkdir(parents=True)
    for name in ("accuracy.py", "helpers.py", "externals.py", "orthohmm.py"):
        (core_dir / name).write_text("# Synthetic source identity fixture; never executed.\n")
    core = [dict(trace.record(core_dir / name), absolute_path=str((core_dir / name).resolve()))
            for name in ("accuracy.py", "helpers.py", "externals.py", "orthohmm.py")]
    for r in core:
        r["path"] = "orthohmm/" + Path(r["absolute_path"]).name
    baseline = save(tmp_path / "baseline.json", dict(core_sources=core, core_commit="synthetic-not-native"))
    plan = save(tmp_path / "plan.json", dict(schema="native_factorial_cost_plan_v1", baseline=baseline, core_commit="synthetic-not-native"))
    views, candidate_ref = [], None
    for cell in trace.CELLS:
        working = tmp_path / cell / "native/orthohmm_working_res"
        checkpoint = working / "high_sensitivity_checkpoint"
        checkpoint.mkdir(parents=True)
        name_path = checkpoint / "gene_names.txt"
        name_path.write_text("\n".join(names) + "\n")
        candidate = working / "orthohmm_edges_clustered.txt"
        candidate.write_text(" ".join(names) + "\n")
        if candidate_ref is None:
            candidate_ref = trace.record(candidate)
        graph = working / "orthohmm_edges.txt"
        graph.write_text("".join("g000\t%s\t1.0\n" % gene for gene in names[1:]))
        metadata = dict(native_factorial=dict(cell=cell, profile_expansion=False, candidate_expansion=False,
                         frozen_pipeline_sha256=core[-1]["sha256"]))
        metric = save(working.parent.parent / "metrics.json", dict(status="complete", metadata=metadata,
                      counts=dict(genes=len(names), orthogroups=1, network_edges=len(names) - 1)))
        inventory = save(tmp_path / (cell + "_inventory.json"), dict(cell=cell, input_genes=len(names),
                          native_outputs_validated=True, checked_files=[metric, trace.record(candidate)]))
        views.append(dict(cell=cell, genes=len(names), native_output_validation=inventory,
                          files={"gene_names.txt": trace.record(name_path)}))
    cases, local = [], []
    for i, (a, b) in enumerate(pairs):
        label = "TP" if i % 2 else "FP"
        case = dict(family="fixture", protein_a=a, protein_b=b, before=label,
                    after="FN" if label == "TP" else "TN", gene_a=a, gene_b=b)
        local.append(dict(case, source_family="Family0000000"))
        forward = [dict(row=i, score=1.0)] if a == "g000" else []
        case["direct_search"] = {cell: dict(support="one_direction" if forward else "no_direct_hit",
                                            gene_a_to_b=forward, gene_b_to_a=[]) for cell in trace.CELLS}
        cases.append(case)
    ledger = tmp_path / "localized.tsv"
    with ledger.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(local[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(local)
    localized = save(tmp_path / "localized.json", dict(flags,
        schema="native_qfo_swiss_reconciliation_localization_v1", source=trace.record(Path(trace.__file__).with_name("trace_native_qfo_swiss_reconciliation.py")),
        pair_ledger=trace.record(ledger), changed_pairs_traced=2023,
        reconstruction=dict(reconstructed_candidate_sha256=candidate_ref["sha256"], source_families=1)))
    search = save(tmp_path / "search.json", dict(flags, schema="native_qfo_swiss_direct_search_support_v1",
        source=trace.record(Path(trace.__file__).with_name("trace_native_qfo_swiss_search_support.py")),
        localization=localized, checkpoints=views, cases=cases, changed_pairs=2023))
    readback = save(tmp_path / "readback.json", dict(flags, schema="native_qfo_swiss_search_support_code_readback_v1",
        source=trace.record(Path(trace.__file__).with_name("readback_native_qfo_swiss_search_support.py")), report=search))
    return search, readback, plan, tmp_path / "output"


def run_fixture(refs):
    search, readback, plan, output = refs
    return trace.run(search["path"], search["sha256"], readback["path"], readback["sha256"], plan["path"], plan["sha256"], output)


def test_complete_synthetic_source_bound_export_and_independent_readback(synthetic):
    result = run_fixture(synthetic)
    output = synthetic[-1]
    ref = trace.record(output / "report.json")
    verified = reader.verify(ref["path"], ref["sha256"])
    assert result["changed_pairs"] == 2023 and len(result["summary"]) == 36
    assert result["graph_bytes_identical"] is result["selected_induced_graph_records_identical"] is True
    assert verified["table_rows_checked"] == 4046 and verified["complete_graph_rows_checked"] == 132
    assert verified["selected_genes_checked"] == 67 and verified["selected_families_checked"] == 1
    assert verified["graph_original_admission_established"] is verified["publication_ready"] is False
    with pytest.raises(ValueError, match="fresh"):
        run_fixture(synthetic)


@pytest.mark.parametrize("fault", ("source", "scope", "readback_binding", "cohort", "views", "candidate_hash", "case_truth",
                                  "case_family", "hit_weight"))
def test_primary_source_bound_refusals(synthetic, fault):
    search_ref, verified_ref, plan_ref, output = synthetic
    search_path = Path(search_ref["path"])
    verified_path = Path(verified_ref["path"])
    search = json.loads(search_path.read_text())
    verified = json.loads(verified_path.read_text())
    if fault == "source":
        search["source"] = trace.record(__file__)
    elif fault == "scope":
        verified["scientific_timings_admitted"] = True
    elif fault == "cohort":
        search["changed_pairs"] -= 1
    elif fault == "views":
        search["checkpoints"].pop()
    elif fault in ("candidate_hash", "case_family"):
        loc_path = Path(search["localization"]["path"])
        localized = json.loads(loc_path.read_text())
        if fault == "candidate_hash":
            localized["reconstruction"]["reconstructed_candidate_sha256"] = "0" * 64
        else:
            table_path = Path(localized["pair_ledger"]["path"])
            table_path.write_text(table_path.read_text().replace("Family0000000", "Family9999999", 1))
            localized["pair_ledger"] = trace.record(table_path)
        search["localization"] = save(loc_path, localized)
    elif fault == "case_truth":
        search["cases"][0]["after"] = "TP"
    elif fault == "hit_weight":
        search["cases"][0]["direct_search"][trace.CELLS[0]]["gene_a_to_b"][0]["score"] = 2.0
    search_ref = save(search_path, search)
    verified["report"] = search_ref
    if fault == "readback_binding":
        verified["report"] = trace.record(__file__)
    verified_ref = save(verified_path, verified)
    with pytest.raises(ValueError):
        run_fixture((search_ref, verified_ref, plan_ref, output))


@pytest.mark.parametrize("fault", ("scope", "source", "case", "edge", "path", "distance", "summary", "table", "candidate", "graph_ref"))
def test_independent_readback_refuses_synthetic_corruption(synthetic, fault):
    run_fixture(synthetic)
    path = synthetic[-1] / "report.json"
    report = json.loads(path.read_text())
    view = report["cases"][0]["views"][trace.CELLS[0]]
    if fault == "scope":
        report["graph_original_admission_established"] = True
    elif fault == "source":
        report["source"] = trace.record(__file__)
    elif fault == "case":
        report["cases"][0]["protein_a"] = "wrong"
    elif fault == "edge":
        report["graphs"][0]["selected_induced_edges"][0]["weight"] = 2.0
    elif fault == "path":
        view["path"] = [report["cases"][0]["gene_a"], "g002", report["cases"][0]["gene_b"]]
    elif fault == "distance":
        view["shortest_path_edges"] = 2
    elif fault == "summary":
        report["summary"][0]["pairs"] += 1
    elif fault == "table":
        table = Path(report["pair_ledger"]["path"])
        table.write_text(table.read_text().replace("direct_weight", "wrong", 1))
        report["pair_ledger"] = trace.record(table)
    elif fault == "candidate":
        report["candidate_partition"] = trace.record(__file__)
    else:
        report["graphs"][0]["graph"] = trace.record(__file__)
    save(path, report)
    with pytest.raises(ValueError):
        reader.verify(path, trace.record(path)["sha256"])


def test_actual_retained_result_and_reader_receipt_without_repeating_graph_scan():
    results = Path(trace.__file__).parent / "results"
    report_path = results / "native_qfo_swiss_graph_support_20261006_v1/report.json"
    result = json.loads(report_path.read_text())
    verified = json.loads((results / "native_qfo_swiss_graph_floyd_readback_20261006.json").read_text())
    assert result["source"] == trace.record(trace.__file__)
    assert verified["source"] == trace.record(reader.__file__) and verified["report"] == trace.record(report_path)
    assert verified["changed_pairs_checked"] == result["changed_pairs"] == 2023
    assert verified["selected_families_checked"] == len(result["selected_families"]) == 23
    assert verified["selected_genes_checked"] == result["selected_genes"] == 1139
    assert verified["table_rows_checked"] == 4046 and verified["complete_graph_rows_checked"] == 51260606
    assert verified["induced_edges_checked"] == 46470 and verified["summary"] == result["summary"]
    assert result["graph_bytes_identical"] is result["selected_induced_graph_records_identical"] is True
    for document in (result, verified):
        assert all(document[k] is False for k in ("graph_original_admission_established", "raw_search_rescanned",
                   "new_scoring_or_admission", "uncertainty_admitted", "scientific_timings_admitted", "publication_ready"))
    for cell in trace.CELLS:
        view = [r for r in result["summary"] if r["cell"] == cell]
        assert sum(r["pairs"] for r in view if r["before"] == "TP" and r["graph_support"] == "direct_edge") == 270
        assert sum(r["pairs"] for r in view if r["before"] == "FP" and r["graph_support"] == "direct_edge") == 1137
        assert sum(r["pairs"] for r in view if r["graph_support"] == "disconnected") == 0
        assert sum(r["pairs"] for r in view if r["search_support"] == "no_direct_hit" and r["graph_support"] == "indirect_path") == 539
