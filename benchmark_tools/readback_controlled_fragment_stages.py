"""Independently read retained fragment stages; do not import their producer."""

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
from itertools import combinations, product
import json
import math
from pathlib import Path

import numpy as np

from benchmark_tools.orthofinder_mcl_to_orthogroups import iter_mcl_clusters, load_sequence_ids


METHODS = ("orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full", "orthofinder_sequence_only")
ARMS = ("baseline", "fragment")


def need(condition, message):
    if not condition:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve()
    data = path.read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def check(ref):
    path = Path(ref["absolute_path"] if "absolute_path" in ref else ref["path"])
    observed = pin(path)
    need((observed["bytes"], observed["sha256"]) == (ref["bytes"], ref["sha256"]), "Changed input: " + str(path))
    return path


def partition(groups, universe):
    lookup = {}
    for label, genes in groups.items():
        need(bool(genes), "Empty group")
        for gene in genes:
            need(gene in universe and gene not in lookup, "Invalid or repeated group member")
            lookup[gene] = label
    need(set(lookup) == set(universe), "Incomplete group partition")
    return lookup


def graph(path, universe):
    neighbors = {g: set() for g in universe}
    edges = set()
    with path.open() as stream:
        for row in csv.reader(stream, delimiter="\t"):
            need(len(row) == 3 and row[0] in neighbors and row[1] in neighbors, "Invalid graph row")
            value = float(row[2])
            need(math.isfinite(value) and value > 0, "Invalid graph weight")
            pair = tuple(sorted(row[:2]))
            need(pair not in edges, "Duplicate graph edge")
            edges.add(pair)
            neighbors[row[0]].add(row[1])
            neighbors[row[1]].add(row[0])
    components = {}
    for start in neighbors:
        if start in components:
            continue
        pending = [start]
        components[start] = start
        while pending:
            for gene in neighbors[pending.pop()]:
                if gene not in components:
                    components[gene] = start
                    pending.append(gene)
    return edges, components


def ancestor(nodes, a, b):
    def path(gene):
        result = []
        while gene:
            need(gene in nodes and gene not in result, "Invalid ancestor path")
            result.append(gene)
            gene = nodes[gene]["parent_node_id"]
        return result
    left, right = path(a), path(b)
    common = next((n for n in right if n in set(left)), None)
    need(common is not None, "Disconnected pair nodes")
    return nodes[common]


def tree_partition(nodes, genes, owners):
    children = defaultdict(list)
    for name, row in nodes.items():
        children[row["parent_node_id"]].append(name)
    need(len(children[""]) == 1, "Invalid reconciliation root")
    completed, waiting = {}, set(nodes)
    high_pairs = set()
    # Bottom-up readiness, independent of the producer's recursive traversal.
    while waiting:
        ready = [n for n in waiting if all(c in completed for c in children[n])]
        need(bool(ready), "Cyclic or disconnected reconciliation")
        for name in ready:
            row, offspring = nodes[name], children[name]
            members = set(row["genes"].split(","))
            need(members <= genes and set(row["species"].split(",")) == {owners[g] for g in members}, "Invalid node membership/species")
            if not offspring:
                need(members == {name} and row["event"] == row["pair_event"] == "leaf", "Invalid leaf")
                split, groups = False, [members]
            else:
                parts = [completed[c][2] for c in offspring]
                need(len(parts) >= 2 and set.union(*parts) == members and sum(map(len, parts)) == len(members), "Invalid child partition")
                species = [{owners[g] for g in part} for part in parts]
                overlap = set().union(*(x & y for x, y in combinations(species, 2)))
                conflict = any(nodes[c]["species_tree_node"] == row["species_tree_node"] for c in offspring)
                expected = "duplication" if overlap else "uncertain" if conflict else "speciation"
                need(row["pair_event"] == expected and row["event"] == ("duplication" if overlap or conflict else "speciation")
                     and int(row["species_overlap_count"]) == len(overlap) and row["mapping_conflict"] == str(conflict).lower(), "Invalid observed event")
                support = float(row["branch_support"]) if row["branch_support"] else None
                need(support is None or math.isfinite(support) and 0 <= support <= 1, "Invalid branch support")
                confidence = "medium" if expected == "uncertain" or expected == "duplication" and len(overlap) < 2 and (support or 0) < 0.9 else "high"
                need(row["event_confidence"] == confidence, "Invalid node confidence")
                if expected == "speciation":
                    high_pairs.update(tuple(sorted(p)) for x, y in combinations(parts, 2) for p in product(x, y) if owners[p[0]] != owners[p[1]])
                split = bool(overlap and row["species_tree_node"] == "S0000") or any(completed[c][0] for c in offspring)
                groups = [g for c in offspring for g in completed[c][1]] if split else [members]
            completed[name] = split, groups, members
            waiting.remove(name)
    root = children[""][0]
    need(completed[root][2] == genes and all(not p or p in nodes for p in children), "Incomplete tree")
    return completed[root][1], high_pairs


def constrain(groups, high_pairs, events):
    detached, details = {}, []
    for number, event in events:
        source, target = set(event["source_genes"]), set(event["target_genes"])
        universe = set().union(*groups)
        need(source and target and not source & target and source | target <= universe, "Invalid constraint sides")
        supporting = sorted(p for p in high_pairs if (p[0] in source and p[1] in target) or (p[1] in source and p[0] in target))
        details.append(dict(event_index=number, supported=bool(supporting), supporting_pair=list(supporting[0]) if supporting else None,
                            source_genes=sorted(source), target_genes=sorted(target)))
        if not supporting:
            for gene in source:
                need(gene not in detached, "Multiple detached constraints")
                detached[gene] = number
    final = []
    for group in groups:
        buckets = defaultdict(list)
        for gene in group:
            buckets[detached.get(gene, -1)].append(gene)
        final.extend(buckets.values())
    return final, details


def read_context(binding, method, queries):
    execution = json.loads(check(binding["execution"]).read_text())
    parent = METHODS[2] if method == METHODS[3] else method
    run = execution["methods"][parent]
    need(run["status"] == "process_succeeded" and run["exit_code"] == 0, "Unsuccessful execution")
    inventory = {r["absolute_path"]: r for r in run["outputs"]}
    output = Path(binding["output"])
    input_dir = Path(binding["configured"]["argv"][2] if method in METHODS[:2] else binding["configured"]["copy_inputs_from"]
                     if method == METHODS[2] else execution["methods"][parent]["argv"][2])
    owners = {}
    for path in input_dir.glob("*.fasta"):
        with path.open() as stream:
            for line in stream:
                if line.startswith(">"):
                    gene = line[1:].split()[0]
                    need(gene not in owners, "Duplicate input gene")
                    owners[gene] = path.stem
    need(bool(owners), "Missing FASTA ownership")
    def read(path):
        need(str(path) in inventory, "Artifact not bound to execution: " + str(path))
        return check(inventory[str(path)])
    def unique(pattern):
        paths = list(output.glob(pattern))
        need(len(paths) == 1, "Nonunique artifact " + pattern)
        return read(paths[0])
    if method in METHODS[2:]:
        mapping = load_sequence_ids(unique("**/SequenceIDs.txt"))
        need(set(mapping.values()) == set(owners) and len(mapping) == len(owners), "Changed MCL input universe")
        seeds = partition({str(i): [mapping[g] for g in cluster] for i, cluster in enumerate(iter_mcl_clusters(unique("**/clusters_OrthoFinder_I*.txt_id_pairs.txt")))}, owners)
        predictions = set()
        if method == METHODS[2]:
            for path in output.glob("**/Orthologues/*/*.tsv"):
                with read(path).open() as stream:
                    rows = csv.reader(stream, delimiter="\t")
                    header = next(rows)
                    need(len(header) == 3, "Invalid OrthoFinder pair header")
                    for row in rows:
                        need(len(row) == 3, "Invalid OrthoFinder pair row")
                        for a, b in product(row[1].split(", "), row[2].split(", ")):
                            need(owners[a] == header[1] and owners[b] == header[2], "OrthoFinder species mismatch")
                            predictions.add(tuple(sorted((a, b))))
        return dict(seeds=seeds, predictions=predictions)
    working = output / "orthohmm_working_res"
    cp = working / "high_sensitivity_checkpoint"
    manifest = json.loads(read(cp / "manifest.json").read_text())
    for name, ref in manifest["files"].items():
        check(dict(ref, path=str(read(cp / name))))
    names = (cp / "gene_names.txt").read_text().splitlines()
    need(set(names) == set(owners) and len(names) == len(owners), "Changed checkpoint input universe")
    species, q, t, score = [np.load(cp / (n + ".npy"), allow_pickle=False) for n in ("gene_to_species", "hit_queries", "hit_targets", "hit_scores")]
    need(len(names) == manifest["genes"] == len(species) and len(q) == len(t) == len(score) == manifest["hits"], "Changed checkpoint sizes")
    code_owners = defaultdict(set)
    for name, code in zip(names, species):
        code_owners[int(code)].add(owners[name])
    need(all(len(v) == 1 for v in code_owners.values()) and len(code_owners) == len(set(owners.values())), "Changed species codes")
    hits = defaultdict(list)
    for offset, (i, j, value) in enumerate(zip(q, t, score)):
        need(0 <= i < len(names) and 0 <= j < len(names) and math.isfinite(float(value)), "Invalid numeric hit")
        hits[names[int(i)], names[int(j)]].append(dict(row=offset, score=float(value)))
    edges, components = graph(read(working / "orthohmm_edges.txt"), owners)
    seeds = {}
    for line in read(output / "orthohmm_orthogroups.txt").read_text().splitlines():
        label, separator, members = line.partition(": ")
        need(separator and label not in seeds, "Invalid named groups")
        seeds[label] = members.split()
    context = dict(seeds=partition(seeds, owners), hits=hits, edges=edges, components=components)
    if method == METHODS[0]:
        return context
    phy = output / "orthohmm_phylogeny"
    predictions = set()
    with read(phy / "orthohmm_pairwise_orthologs.tsv").open() as stream:
        rows = csv.DictReader(stream, delimiter="\t")
        for row in rows:
            need(owners[row["gene_a"]] == row["species_a"] and owners[row["gene_b"]] == row["species_b"], "Native pair species mismatch")
            predictions.add(tuple(sorted((row["gene_a"], row["gene_b"]))))
    candidates = {f"Family{i:07d}": line.split() for i, line in enumerate(read(working / "phylogeny_candidate_superfamilies.txt").read_text().splitlines())}
    candidate_index = partition(candidates, owners)
    root_groups, root_sources = {}, {}
    with read(phy / "orthohmm_root_hogs.tsv").open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            root_groups[row["root_hog"]] = row["genes"].split(",")
            root_sources[row["root_hog"]] = row["source_family"]
    roots = partition(root_groups, owners)
    provenance = json.loads(read(phy / "provenance_manifest.json").read_text())
    need(provenance["pair_orthology_rule"] == "positive_paralogy" and provenance["root_duplication_rule"] == "species_overlap", "Changed event rules")
    events = json.loads(read(working / "phylogeny_candidate_merges.json").read_text())
    active = provenance["membership_reconciliation"] is not None
    need(active == bool(events), "Changed constraint activation")
    nodes = defaultdict(dict)
    with read(phy / "orthohmm_reconciliation_nodes.tsv").open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            need(row["node_id"] not in nodes[row["source_family"]], "Duplicate node")
            nodes[row["source_family"]][row["node_id"]] = row
    before, detachments = {}, {}
    for family in {candidate_index[a] for a, b in queries if candidate_index[a] == candidate_index[b]}:
        genes = set(candidates[family])
        if nodes[family]:
            groups, high = tree_partition(nodes[family], genes, owners)
        else:
            counts = Counter(owners[g] for g in genes)
            need(len(genes) < 3 or len(counts) < 2 or max(counts.values()) < 2, "Missing ambiguous tree")
            need(str(phy / "gene_trees" / (family + ".raw.nwk")) not in inventory, "Missing node evidence for inferred tree")
            groups = [genes]
            high = {tuple(sorted(p)) for p in combinations(genes, 2) if owners[p[0]] != owners[p[1]]}
        before[family] = partition({str(i): g for i, g in enumerate(groups)}, genes)
        constraints = [(i, e) for i, e in enumerate(events) if set(e["source_genes"]) <= genes]
        final, details = constrain(groups, high, constraints) if active else (groups, [])
        need({frozenset(g) for g in final} == {frozenset(root_groups[k]) for k in root_groups if root_sources[k] == family}, "RootHOG reconstruction mismatch")
        detachments[family] = [d for d in details if not d["supported"]]
    context.update(predictions=predictions, candidates=candidate_index, roots=roots, nodes=nodes, active=active,
                   before=before, detachments=detachments)
    return context


def observation(method, context, pair):
    a, b = pair
    same = context["seeds"][a] == context["seeds"][b]
    predicted = same if method in (METHODS[0], METHODS[3]) else pair in context["predictions"]
    row = dict(predicted=predicted, same_seed_group=same, search_status="unavailable_matched_adapter",
               hit_forward=None, hit_reverse=None, graph_direct=None, graph_connected=None,
               same_candidate=None, same_root_hog=None, pair_event=None, membership_filter_active=None, detachment_events=[])
    if method in METHODS[:2]:
        row.update(search_status="significant_hit_checkpoint", hit_forward=bool(context["hits"][a,b]), hit_reverse=bool(context["hits"][b,a]),
                   directed_hits=dict(forward=context["hits"][a,b], reverse=context["hits"][b,a]), graph_direct=pair in context["edges"],
                   graph_connected=context["components"][a] == context["components"][b])
    if method == METHODS[0]:
        location = "group_retention" if predicted else "group_separation"
    elif method in METHODS[2:]:
        location = "native_pair_retention" if predicted else "within_mcl_native_pair_exclusion" if same else "mcl_group_separation"
    else:
        family = context["candidates"][a]
        together = family == context["candidates"][b]
        row.update(same_candidate=together, same_root_hog=context["roots"][a] == context["roots"][b],
                   candidate_families=[context["candidates"][g] for g in pair], membership_filter_active=context["active"])
        if not together:
            need(not predicted, "Predicted cross-candidate pair")
            location = "candidate_separation"
        else:
            if context["nodes"][family]:
                node = ancestor(context["nodes"][family], a, b)
                row["pair_event"] = node["pair_event"]
                row["pair_node"] = {k: node[k] for k in ("node_id", "event", "pair_event", "mapping_conflict", "event_confidence")}
                row["pair_node"].update(species_overlap_count=int(node["species_overlap_count"]), branch_support=float(node["branch_support"]) if node["branch_support"] else None)
            else:
                row["pair_event"] = "unambiguous_bypass"
            row.update(same_pre_constraint_root=context["before"][family][a] == context["before"][family][b],
                       detachment_events=[d for d in context["detachments"][family] if a in d["source_genes"] or b in d["source_genes"]])
            allowed = row["pair_event"] != "duplication"
            need(predicted == (allowed and (not context["active"] or row["same_root_hog"])), "Native event/root rule mismatch")
            location = ("unambiguous_bypass_retention" if row["pair_event"] == "unambiguous_bypass" else "event_rule_retention") if predicted else (
                "observed_duplication_exclusion" if not allowed else "unsupported_satellite_separation" if row["detachment_events"] else "root_partition_exclusion")
    row["observed_location"] = location
    return row


def verify(report_path, expected_sha):
    ref = pin(report_path)
    need(ref["sha256"] == expected_sha, "Changed stage report")
    report = json.loads(Path(report_path).read_text())
    need(report["status"] == "retained_stages_verified" and report["new_inference_or_scoring"] is False, "Invalid stage scope")
    selected = json.loads(check(report["selection"]).read_text())
    checked = [pin(check(r)) for r in report["checked_inputs"]]
    bindings = {b["seed"]: b["arms"] for b in selected["bindings"]}
    needs = defaultdict(set)
    for case in selected["cases"]:
        for arm in ARMS:
            needs[case["method"], case["seed"], arm].add((case["gene_a"], case["gene_b"]))
    contexts = {key: read_context(bindings[key[1]][key[2]][key[0]], key[0], pairs) for key, pairs in needs.items()}
    expected_rows, locations = [], Counter()
    need(len(report["cases"]) == len(selected["cases"]), "Changed case count")
    for case, original in zip(report["cases"], selected["cases"]):
        need({k:v for k,v in case.items() if k != "stages"} == original, "Changed representative identity")
        for arm in ARMS:
            actual = observation(case["method"], contexts[case["method"], case["seed"], arm], (case["gene_a"], case["gene_b"]))
            need(actual == case["stages"][arm], "Stage mismatch " + case["case_id"] + "/" + arm)
            need(actual["predicted"] is case["comparator_predictions"][arm][case["method"]], "Comparator label mismatch")
            expected_rows.append(dict(case_id=case["case_id"], method=case["method"], seed=case["seed"], fragment_endpoints=case["fragment_endpoints"],
                                      category=case["category"], truth=case["truth"], arm=arm, **actual))
            locations[case["method"], case["category"], actual["observed_location"], arm] += 1
    need(report["summary"] == [dict(method=m, category=c, observed_location=l, arm=a, representatives=n) for (m,c,l,a),n in sorted(locations.items())], "Summary mismatch")
    with Path(report_path).with_name("stages.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    need(len(rows) == len(expected_rows) == report["stage_rows"], "Changed stage TSV row count")
    for row, expected in zip(rows, expected_rows):
        need(row == {k: "" if expected[k] is None else str(expected[k]) for k in row}, "Stage TSV value mismatch")
    for r in checked:
        check(r)
    return dict(status="independent_native_stages_verified", report=ref, stage_rows=len(rows), selected_cases=len(report["cases"]),
                checked_input_files=len(checked), native_contexts=len(contexts), summary=report["summary"],
                source=pin(__file__), new_inference_or_scoring=False, publication_ready=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    need(not args.output.exists(), "Readback output already exists")
    result = verify(args.report, args.sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: result[k] for k in ("status", "stage_rows", "selected_cases", "native_contexts")}))
