"""Read selected fragment cases through retained native pipeline stages only."""

import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path

from Bio import SeqIO
from orthohmm import accuracy
from benchmark_tools import audit_accuracy_checkpoint as numeric
from benchmark_tools import trace_simulation_upstream as upstream
from benchmark_tools import probe_simulation_gene_tree_oracle as partitions
from benchmark_tools import trace_native_qfo_swiss_reconciliation as reconciliation
from benchmark_tools import reconstruct_reconciliation_trace as reconstruction
from benchmark_tools import simulation_method_outputs as native
from benchmark_tools import select_controlled_fragment_trace as selection
from benchmark_tools.export_native_qfo_three_cell_strata import record, require, write_tsv
from benchmark_tools.score_ygob_groups import read_predictions, membership
from benchmark_tools.report_ygob_validation import read_checkpoint


METHODS, ARMS = selection.METHODS, selection.ARMS
FIELDS = ("case_id", "method", "seed", "fragment_endpoints", "category", "truth", "arm", "predicted",
          "search_status", "hit_forward", "hit_reverse", "graph_direct", "graph_connected", "same_seed_group",
          "same_candidate", "same_root_hog", "pair_event", "membership_filter_active", "observed_location")


def index(groups, universe):
    value = membership(groups)
    require(set(value) == set(universe), "Stage groups do not partition the complete gene universe")
    return value


def pair_location(method, row):
    predicted = row["predicted"]
    if method == METHODS[0]:
        require(predicted is row["same_seed_group"], "Group-derived native prediction disagrees")
        return "group_retention" if predicted else "group_separation"
    if method in METHODS[2:]:
        if method == METHODS[3]:
            require(predicted is row["same_seed_group"], "MCL checkpoint prediction disagrees")
        return ("native_pair_retention" if predicted else "within_mcl_native_pair_exclusion"
                if row["same_seed_group"] else "mcl_group_separation")
    if not row["same_candidate"]:
        require(predicted is False, "Native phylogeny pair crosses candidate boundary")
        return "candidate_separation"
    require(row["pair_event"] in ("speciation", "uncertain", "duplication", "unambiguous_bypass"), "Missing pair event")
    raw = row["pair_event"] != "duplication"
    expected = raw and (not row["membership_filter_active"] or row["same_root_hog"])
    require(predicted is expected, "Native prediction differs from retained event/root-membership rule")
    if predicted:
        return "unambiguous_bypass_retention" if row["pair_event"] == "unambiguous_bypass" else "event_rule_retention"
    if not raw:
        return "observed_duplication_exclusion"
    return "unsupported_satellite_separation" if row["detachment_events"] else "root_partition_exclusion"


def artifact_reader(binding, method, evidence):
    execution = json.loads(selection.checked(binding["execution"], evidence).read_text())
    parent = METHODS[2] if method == METHODS[3] else method
    executed = execution["methods"][parent]
    require(executed["status"] == "process_succeeded" and executed["exit_code"] == 0, "Stage execution did not succeed")
    artifacts = {str(Path(r["absolute_path"]).resolve()): r for r in executed["outputs"]}
    require(len(artifacts) == len(executed["outputs"]), "Duplicate native artifact inventory")
    output = Path(binding["output"]).resolve()
    def read(path):
        path = Path(path).resolve()
        require(str(path) in artifacts, "Stage artifact absent from native execution: " + str(path))
        return selection.checked(artifacts[str(path)], evidence)
    return output, read, artifacts


def root_groups(path, candidates, universe):
    groups, sources = {}, {}
    with path.open() as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        require(rows.fieldnames == ["root_hog", "source_family", "genes"], "Changed RootHOG schema")
        for row in rows:
            require(None not in row and all(v is not None for v in row.values()), "Malformed RootHOG row")
            name, family = row["root_hog"], row["source_family"]
            require(name not in groups and family in candidates, "Invalid RootHOG identity")
            genes = row["genes"].split(",")
            require(set(genes) <= candidates[family], "RootHOG crosses candidate boundary")
            groups[name], sources[name] = genes, family
    value = index(groups, universe)
    require(all(set().union(*(set(groups[g]) for g in groups if sources[g] == f)) == members
                for f, members in candidates.items()), "Incomplete candidate RootHOG coverage")
    return groups, sources, value


def hmm_context(binding, method, queries, owners, metrics_ref, evidence):
    output, read, artifacts = artifact_reader(binding, method, evidence)
    metric_path = Path(binding["configured"]["metrics"])
    metrics = json.loads(selection.checked(dict(metrics_ref, absolute_path=str(metric_path)), evidence).read_text())
    require(metrics["status"] == "complete" and Path(metrics["metadata"]["output_directory"]).resolve() == output,
            "Changed completed metric binding")
    working = output / "orthohmm_working_res"
    checkpoint = working / "high_sensitivity_checkpoint"
    manifest = json.loads(read(checkpoint / "manifest.json").read_text())
    expected = {"gene_names.txt", "gene_to_species.npy", "hit_queries.npy", "hit_targets.npy", "hit_scores.npy"}
    require(manifest["complete"] is True and set(manifest["files"]) == expected, "Incomplete numeric search checkpoint")
    for name, ref in manifest["files"].items():
        path = read(checkpoint / name)
        selection.checked(dict(ref, absolute_path=str(path)), evidence)
    names, species, q, t, scores = accuracy.load_accuracy_checkpoint(checkpoint, verify=False)
    summary = numeric.check_arrays(names, species, q, t, scores)
    require(summary["genes"] == manifest["genes"] and summary["hits"] == manifest["hits"]
            and set(names) == set(owners) and len(q) == metrics["counts"]["significant_hits"], "Changed search universe/count")
    species_codes = defaultdict(set)
    for name, code in zip(names, species):
        species_codes[int(code)].add(owners[name])
    require(all(len(s) == 1 for s in species_codes.values())
            and len(species_codes) == len(set(owners.values())), "Changed checkpoint ownership")
    hits = defaultdict(list)
    wanted = {p for a, b in queries for p in ((a, b), (b, a))}
    for i, (a, b, score) in enumerate(zip(q, t, scores)):
        pair = names[int(a)], names[int(b)]
        if pair in wanted:
            hits[pair].append(dict(row=i, score=float(score)))
    edges, components, component_count = upstream.graph(read(working / "orthohmm_edges.txt"), names)
    require(len(edges) == metrics["counts"]["network_edges"], "Changed graph edge count")
    groups = read_predictions(read(output / "orthohmm_orthogroups.txt"), "named_groups")
    seeds = index(groups, owners)
    context = dict(hits=hits, edges=edges, components=components, seeds=seeds,
                   numeric=summary, graph_components=component_count, metrics=metrics)
    if method == METHODS[0]:
        return context
    candidate_path = read(working / "phylogeny_candidate_superfamilies.txt")
    candidates = partitions.partition(candidate_path, set(owners))
    candidate_index = index(candidates, owners)
    events = json.loads(read(working / "phylogeny_candidate_merges.json").read_text())
    sidecar = upstream.seed_sidecar(read(working / "phylogeny_candidate_seeds.tsv"), candidates,
                                    metrics["metadata"]["phylogeny_candidate_profile"], events)
    doc = json.loads(read(output / "orthohmm_phylogeny/provenance_manifest.json").read_text())
    require(doc["input_cluster_sha256"] == record(candidate_path)["sha256"]
            and doc["pair_orthology_rule"] == "positive_paralogy" and doc["root_duplication_rule"] == "species_overlap",
            "Changed candidate input or native event rule")
    active = doc["membership_reconciliation"] is not None
    require(active is bool(events), "Membership filter activation differs from saved constraints")
    roots, sources, root_index = root_groups(read(output / "orthohmm_phylogeny/orthohmm_root_hogs.tsv"), candidates, owners)
    node_path = read(output / "orthohmm_phylogeny/orthohmm_reconciliation_nodes.tsv")
    with node_path.open() as handle:
        node_families = {row["source_family"] for row in csv.DictReader(handle, delimiter="\t")}
    wanted_families = {candidate_index[a] for a, b in queries if candidate_index[a] == candidate_index[b]}
    require(node_families <= set(candidates), "Unknown node source family")
    nodes, node_rows = reconciliation.read_nodes(node_path, {f:candidates[f] for f in wanted_families & node_families})
    reconstructed, detachments = {}, {}
    for family in sorted(wanted_families):
        genes = candidates[family]
        if family in nodes:
            raw_rows = [dict(row, genes=",".join(sorted(row["genes"])), species=",".join(sorted(row["species"])))
                        for row in nodes[family].values()]
            before = reconstruction.reconstruct_nodes(raw_rows, genes, owners)
        else:
            require(str(output / f"orthohmm_phylogeny/gene_trees/{family}.raw.nwk") not in artifacts,
                    "Reconciled family lacks node evidence")
            before = reconstruction.reconstruct_bypass(genes, owners)
        constraints = [(i, event) for i, event in enumerate(events) if set(event["source_genes"]) <= genes]
        final, details = reconstruction.apply_logged_constraints(before, constraints, genes) if active else (before["root_groups"], [])
        observed = {frozenset(roots[g]) for g in roots if sources[g] == family}
        require({frozenset(g) for g in final} == observed, "Recorded RootHOGs differ from node/constraint reconstruction")
        reconstructed[family] = index({str(i):g for i,g in enumerate(before["root_groups"])}, genes)
        detachments[family] = [event for event in details if not event["supported"]]
    context.update(candidates=candidate_index, roots=root_index, nodes=nodes, active=active,
                   pre_constraint_roots=reconstructed, detachments=detachments,
                   candidate_seed_inventory=sidecar, reconciliation_node_rows=node_rows)
    return context


def observe(method, context, pair, predicted):
    a, b = pair
    row = dict(predicted=predicted, search_status="significant_hit_checkpoint" if method in METHODS[:2] else "unavailable_matched_adapter",
        hit_forward=None, hit_reverse=None, graph_direct=None, graph_connected=None,
        same_seed_group=context["seeds"][a] == context["seeds"][b], same_candidate=None, same_root_hog=None,
        pair_event=None, membership_filter_active=None, detachment_events=[])
    if method in METHODS[:2]:
        row.update(hit_forward=bool(context["hits"][a,b]), hit_reverse=bool(context["hits"][b,a]),
            directed_hits=dict(forward=context["hits"][a,b], reverse=context["hits"][b,a]),
            graph_direct=pair in context["edges"], graph_connected=context["components"][a] == context["components"][b])
        require(not row["graph_direct"] or row["graph_connected"], "Direct edge without connectivity")
    if method == METHODS[1]:
        family = context["candidates"][a]
        same = family == context["candidates"][b]
        row.update(same_candidate=same, candidate_families=[context["candidates"][g] for g in pair],
            same_root_hog=context["roots"][a] == context["roots"][b], membership_filter_active=context["active"])
        if same:
            if family in context["nodes"]:
                node = reconciliation.pair_lca(context["nodes"][family], a, b)
                row.update(pair_event=node["pair_event"], pair_node={k:node[k] for k in (
                    "node_id", "event", "pair_event", "species_overlap_count", "mapping_conflict", "event_confidence", "branch_support")})
            else:
                row["pair_event"] = "unambiguous_bypass"
            row.update(same_pre_constraint_root=context["pre_constraint_roots"][family][a] == context["pre_constraint_roots"][family][b],
                detachment_events=[e for e in context["detachments"][family] if a in e["source_genes"] or b in e["source_genes"]])
    row["observed_location"] = pair_location(method, row)
    return row


def metrics_evidence(outcome, dataset, arm, method):
    if arm == "baseline":
        matches = [r for r in dataset["baseline_admission"] if r["method"] == method]
        require(len(matches) == 1, "Missing or duplicate original metric admission")
        admission = matches[0]["admission"]
    else:
        require(arm == "fragment", "Unknown trace arm")
        admission = outcome["admission"]
    require(admission["status"] == "admitted", "Metric outcome was not admitted")
    return admission["native_validation"]["metrics"]


def run(root, selection_path, selection_sha, output):
    root, output = Path(root).resolve(), Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    require(output.parent.resolve() == root / "benchmark_tools/results", "Use established results directory")
    evidence = {}
    selected_ref = record(selection_path)
    require(selected_ref["sha256"] == selection_sha, "Changed selected representatives")
    selected = json.loads(selection.checked(selected_ref, evidence).read_text())
    require(selected["schema"] == "controlled_fragment_trace_selection_v1" and selected["status"] == "selected_untraced"
            and selected["planned_bins"] == 84 and selected["source"] == record(selection.__file__)
            and selected["new_inference_or_scoring"] is False and selected["publication_ready"] is False, "Wrong selected scope")
    require(selected["protocol"]["sha256"] == selection.PINS["protocol"][1], "Changed trace protocol")
    selection.checked(selected["protocol"], evidence)
    refs = {Path(r["path"]).name:r for r in selected["checked_inputs"] if r["path"].endswith("controlled_fragment_results_20261009_v1/report.json")}
    require(len(refs) == 1 and refs["report.json"]["sha256"] == selection.PINS["result"][1], "Mixed original result binding")
    report = json.loads(selection.checked(refs["report.json"], evidence).read_text())
    prepared = json.loads(selection.checked(report["manifest"], evidence).read_text())
    datasets = {d["seed"]:d for d in prepared["datasets"]}
    records = {(r["arm"],r["seed"],r["method"]):r for r in report["records"]}
    bindings = {r["seed"]:r["arms"] for r in selected["bindings"]}
    query_sets = defaultdict(set)
    for case in selected["cases"]:
        for arm in ARMS:
            query_sets[arm,case["seed"],case["method"]].add((case["gene_a"],case["gene_b"]))
    observations = {}
    for (arm, seed, method), queries in sorted(query_sets.items()):
        data, binding = datasets[seed], bindings[seed][arm][method]
        owners = {}
        for ref in data["verified_inputs"]["inputs"]:
            path = selection.checked(ref, evidence)
            for sequence in SeqIO.parse(path, "fasta"):
                require(sequence.id not in owners, "Repeated input gene")
                owners[sequence.id] = path.stem
        if method in METHODS[:2]:
            metric_ref = metrics_evidence(records[arm,seed,method], data, arm, method)
            context = hmm_context(binding, method, queries, owners, metric_ref, evidence)
        else:
            output_dir, read, _ = artifact_reader(binding, method, evidence)
            clusters = native.unique_path(output_dir, "**/clusters_OrthoFinder_I*.txt_id_pairs.txt")
            mapping = native.unique_path(output_dir, "**/SequenceIDs.txt")
            context = dict(seeds=index(read_checkpoint(read(clusters), read(mapping), owners), owners))
        for case in selected["cases"]:
            if case["seed"] == seed and case["method"] == method:
                pair = case["gene_a"],case["gene_b"]
                observations[case["case_id"],arm] = observe(method, context, pair, case["comparator_predictions"][arm][method])
    cases = [dict(case, stages={a:observations[case["case_id"],a] for a in ARMS}) for case in selected["cases"]]
    summary = Counter((c["method"],c["category"],c["stages"][a]["observed_location"],a) for c in cases for a in ARMS)
    for module in (selection, native, numeric, accuracy, upstream, partitions, reconciliation, reconstruction):
        selection.checked(record(module.__file__), evidence)
    selection.checked(record(__file__), evidence)
    for ref in evidence.values():
        require(record(ref["path"]) == ref, "Stage input changed during read-only trace")
    output.mkdir(parents=True)
    rows = [dict(case_id=c["case_id"], method=c["method"], seed=c["seed"],fragment_endpoints=c["fragment_endpoints"],
        category=c["category"],truth=c["truth"],arm=a,**{k:c["stages"][a][k] for k in FIELDS[7:]}) for c in cases for a in ARMS]
    write_tsv(output / "stages.tsv", rows, FIELDS)
    result = dict(schema="controlled_fragment_stage_trace_v1", status="retained_stages_verified", selection=selected_ref,
        source=record(__file__), cases=cases, stage_rows=len(rows), checked_inputs=list(evidence.values()),
        summary=[dict(method=m,category=c,observed_location=l,arm=a,representatives=n) for (m,c,l,a),n in sorted(summary.items())],
        new_inference_or_scoring=False,new_bootstrap_draws=0,publication_ready=False,
        uncertainty_admitted=False,scientific_timings_admitted=False,independent_confirmation=False,
        limitations=["Deterministic retrospective representatives, not prevalence or independent biological evidence.",
            "Missing significant hits do not identify prefilter versus scoring failures; graph/group support is not orthology.",
            "Node calls locate a recorded exclusion or retention, not proof of inferred-tree or biological correctness.",
            "Matched OrthoFinder search/reconciliation adapters unavailable; only its MCL/native pair transition is traced.",
            "No initial-HMM-off control; synthetic fragment, runtime-inventory and shared-host limitations remain."])
    with (output / "report.json").open("x") as stream:
        json.dump(result,stream,indent=2,sort_keys=True,allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("root","selection","output"):
        parser.add_argument("--"+name,required=True,type=Path)
    parser.add_argument("--selection-sha256",required=True)
    args=parser.parse_args()
    result=run(args.root,args.selection,args.selection_sha256,args.output)
    print(json.dumps({k:result[k] for k in ("status","stage_rows")}))
