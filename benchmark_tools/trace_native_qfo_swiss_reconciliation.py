"""Locate every changed native SwissTrees pair in observed reconciliation records."""

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import readback_native_qfo_swiss_transitions as readback
from benchmark_tools.export_native_factorial_progress import load, require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from qfo_benchmark.og_to_pairwise import _strip_to_uniprot

NODE_HEADER = ("source_family", "node_id", "parent_node_id", "event", "species_tree_node", "species", "genes",
               "pair_event", "event_confidence", "species_overlap_count", "mapping_conflict", "branch_support")


def source_groups(path, targets, expected_hash):
    digest, current, genes = hashlib.sha256(), None, []
    selected, members = {}, {}
    family_count = root_count = gene_count = 0

    def finish(family, values):
        nonlocal family_count, gene_count
        require(family == f"Family{family_count:07d}" and values and len(set(values)) == len(values),
                "Invalid reconstructed source family")
        digest.update((" ".join(sorted(values)) + "\n").encode())
        family_count += 1
        gene_count += len(values)
        if any(_strip_to_uniprot(g) in targets for g in values):
            members[family] = set(values)

    with path.open(newline="") as stream:
        rows = csv.DictReader(stream, delimiter="\t")
        require(rows.fieldnames == ["root_hog", "source_family", "genes"], "Changed RootHOG header")
        for row in rows:
            require(None not in row and all(v is not None for v in row.values())
                    and row["root_hog"] == f"RootHOG{root_count:07d}", "Invalid RootHOG row/order")
            family = row["source_family"]
            if current is not None and family != current:
                finish(current, genes)
                genes = []
            current = family
            group = row["genes"].split(",")
            require(all(group) and group == sorted(set(group)), "Invalid RootHOG genes")
            genes.extend(group)
            for gene in group:
                accession = _strip_to_uniprot(gene)
                if accession in targets:
                    require(accession not in selected, "Duplicate or ambiguous target accession")
                    selected[accession] = dict(gene=gene, source_family=family, root_hog=row["root_hog"])
            root_count += 1
    require(current is not None, "Empty RootHOG output")
    finish(current, genes)
    require(digest.hexdigest() == expected_hash, "Reconstructed candidate groups differ from R0 input bytes")
    require(set(selected) == targets, "Changed pair endpoint missing from native source groups")
    return selected, members, dict(source_families=family_count, root_hogs=root_count,
        input_genes=gene_count, reconstructed_candidate_sha256=digest.hexdigest())


def validate_nodes(nodes, expected_genes):
    children = defaultdict(list)
    roots = []
    for name, node in nodes.items():
        parent = node["parent_node_id"]
        if parent:
            require(parent in nodes and parent != name, "Missing or self-referential node parent")
            children[parent].append(name)
        else:
            roots.append(name)
    require(len(roots) == 1, "Require one reconciliation root")
    visited, descendants, species = set(), {}, {}
    stack = [(roots[0], False)]
    while stack:
        name, ready = stack.pop()
        node = nodes[name]
        if not ready:
            require(name not in visited, "Cycle or repeated node in reconciliation")
            visited.add(name)
            stack.append((name, True))
            stack.extend((child, False) for child in children[name])
            continue
        require(node["mapping_conflict"] in ("true", "false"), "Invalid mapping conflict flag")
        if not children[name]:
            require(name in expected_genes and node["event"] == node["pair_event"] == "leaf"
                    and node["genes"] == {name} and len(node["species"]) == 1
                    and node["mapping_conflict"] == "false" and node["species_overlap_count"] == 0,
                    "Invalid observed leaf annotation")
            descendants[name], species[name] = {name}, node["species"]
        else:
            require(len(children[name]) >= 2, "Unary reconciliation node")
            child_genes = [descendants[c] for c in children[name]]
            child_species = [species[c] for c in children[name]]
            union = set().union(*child_genes)
            require(sum(map(len, child_genes)) == len(union), "Overlapping descendant genes")
            overlap = set().union(*(a & b for i, a in enumerate(child_species) for b in child_species[i + 1:]))
            mapping = any(nodes[c]["species_tree_node"] == node["species_tree_node"] for c in children[name])
            event = "duplication" if overlap or mapping else "speciation"
            pair_event = "duplication" if overlap else "uncertain" if mapping else "speciation"
            require(node["genes"] == union and node["species"] == set().union(*child_species)
                    and node["species_overlap_count"] == len(overlap)
                    and (node["mapping_conflict"] == "true") is mapping
                    and node["event"] == event and node["pair_event"] == pair_event,
                    "Observed node disagrees with descendants or frozen event rule")
            descendants[name], species[name] = union, node["species"]
        require(node["event_confidence"] in ("high", "medium", "not_applicable"), "Invalid event confidence")
        if node["branch_support"] is not None:
            require(math.isfinite(node["branch_support"]) and 0 <= node["branch_support"] <= 1,
                    "Invalid serialized branch support")
    require(visited == set(nodes) and descendants[roots[0]] == expected_genes,
            "Incomplete or disconnected reconciliation tree")


def read_nodes(path, families):
    selected = {f: {} for f in families}
    observed_rows = 0
    with path.open(newline="") as stream:
        rows = csv.DictReader(stream, delimiter="\t")
        require(tuple(rows.fieldnames or ()) == NODE_HEADER, "Changed reconciliation-node header")
        for row in rows:
            observed_rows += 1
            require(None not in row and all(v is not None for v in row.values()), "Invalid node row width")
            family = row["source_family"]
            if family not in selected:
                continue
            name = row["node_id"]
            require(name and name not in selected[family], "Duplicate reconciliation node")
            row["genes"], row["species"] = set(row["genes"].split(",")), set(row["species"].split(","))
            row["species_overlap_count"] = int(row["species_overlap_count"])
            row["branch_support"] = float(row["branch_support"]) if row["branch_support"] else None
            selected[family][name] = row
    for family, nodes in selected.items():
        validate_nodes(nodes, families[family])
    return selected, observed_rows


def pair_lca(nodes, a, b):
    require(a in nodes and b in nodes and nodes[a]["event"] == nodes[b]["event"] == "leaf",
            "Pair endpoints missing from reconciliation leaves")
    ancestors, cursor = set(), a
    while cursor:
        require(cursor not in ancestors, "Cycle in first ancestor path")
        ancestors.add(cursor)
        cursor = nodes[cursor]["parent_node_id"]
    visited, cursor = set(), b
    while cursor not in ancestors:
        require(cursor and cursor not in visited, "Disconnected or cyclic second ancestor path")
        visited.add(cursor)
        cursor = nodes[cursor]["parent_node_id"]
    return nodes[cursor]


def selected_predictions(path, queries):
    found, count = set(), 0
    with path.open(newline="") as stream:
        rows = csv.reader(stream, delimiter="\t")
        require(next(rows) == ["gene_a", "species_a", "gene_b", "species_b"], "Changed native pair header")
        for row in rows:
            require(len(row) == 4, "Invalid native pair width")
            pair = tuple(sorted((row[0], row[2])))
            if pair in queries:
                found.add(pair)
            count += 1
    return found, count


def run(transition_path, transition_sha, output, ledger):
    require(output != ledger and all(not p.exists() and not p.is_symlink() for p in (output, ledger)),
            "Require distinct fresh output paths")
    readback.verify(transition_path, transition_sha)
    evidence = []
    transition, transition_ref = load(transition_path, transition_sha, evidence)
    evidence.extend(transition["checked_inputs"])
    evidence.extend([transition["changed_relations_ledger"], transition["source"], record(readback.__file__)])
    admissions, reviews = [], []
    for index, row in enumerate(transition["cells"]):
        admission, _ = load(row["admission"]["path"], row["admission"]["sha256"], evidence)
        stage = admission["conversion"]
        require(stage["cell"] == row["cell"] and stage["native_index"] == index + 6,
                "Changed native conversion identity")
        ref = stage["terminal_review"] if index == 0 else stage["scientific_recovery"]
        review, _ = load(ref["path"], ref["sha256"], evidence)
        output_ref = review["reviews"]["outputs_or_failure"] if index == 0 else review["outputs"]
        reviewed, _ = load(output_ref["path"], output_ref["sha256"], evidence)
        require(reviewed["native_outputs_validated"] is True and reviewed["cell"] == row["cell"]
                and stage["native_input"] in reviewed["checked_files"], "Unadmitted native prediction output")
        evidence.extend([review["source"], reviewed["source"]])
        admissions.append(admission)
        reviews.append(reviewed)
    r0 = admissions[0]["conversion"]["native_input"]
    r1 = admissions[1]["conversion"]["native_input"]
    directory = Path(r1["path"]).parent
    def admitted(name):
        refs = [r for r in reviews[1]["checked_files"] if r["path"] == str(directory / name)]
        require(len(refs) == 1, "Missing independently admitted " + name)
        evidence.append(refs[0])
        check(refs[0])
        return refs[0]
    root_ref, manifest_ref = admitted("orthohmm_root_hogs.tsv"), admitted("provenance_manifest.json")
    manifest = json.loads(Path(manifest_ref["path"]).read_text())
    require(manifest["input_cluster_sha256"] == r0["sha256"] and manifest["membership_reconciliation"] is None
            and manifest["pair_orthology_rule"] == "positive_paralogy", "Changed candidate input or pair policy")
    with open(transition["changed_relations_ledger"]["path"], newline="") as stream:
        changes = list(csv.DictReader(stream, delimiter="\t"))
    require(len(changes) == transition["comparison"]["changed_relations"] and changes,
            "Empty or incomplete changed-pair ledger")
    targets = {row[key] for row in changes for key in ("protein_a", "protein_b")}
    evidence.extend([r0, r1])
    for ref in evidence:
        check(ref)
    selected, members, reconstruction = source_groups(Path(root_ref["path"]), targets, r0["sha256"])
    require(reconstruction["input_genes"] == reviews[1]["input_genes"]
            and reconstruction["root_hogs"] == reviews[1]["phylogeny"]["root_hogs"], "Reconstructed universe differs")
    # Node annotations were not inventoried by the original output validator.
    nodes_ref = record(directory / "orthohmm_reconciliation_nodes.tsv")
    evidence.append(nodes_ref)
    nodes, annotation_rows = read_nodes(Path(nodes_ref["path"]), members)
    queries = {tuple(sorted((selected[row["protein_a"]]["gene"], selected[row["protein_b"]]["gene"])))
               for row in changes}
    found, native_rows = selected_predictions(Path(r1["path"]), queries)
    require(not found and native_rows == reviews[1]["phylogeny"]["native_pair_rows"],
            "Removed reference pairs disagree with admitted native predictions")
    traced, summary = [], Counter()
    for row in changes:
        a, b = selected[row["protein_a"]], selected[row["protein_b"]]
        require(a["source_family"] == b["source_family"] and row["before"] in ("TP", "FP")
                and row["after"] in ("FN", "TN"), "Changed pair is not removed from a shared candidate")
        node = pair_lca(nodes[a["source_family"]], a["gene"], b["gene"])
        require(node["pair_event"] == "duplication" and node["species_overlap_count"] > 0,
                "Removed pair lacks positive-paralogy exclusion")
        same_root = a["root_hog"] == b["root_hog"]
        summary[row["before"] + ("_same_root" if same_root else "_different_root")] += 1
        summary[row["before"] + "_" + node["event_confidence"]] += 1
        traced.append(dict(row, source_family=a["source_family"], gene_a=a["gene"], gene_b=b["gene"],
            root_hog_a=a["root_hog"], root_hog_b=b["root_hog"], same_root_hog=same_root,
            lca_node=node["node_id"], pair_event=node["pair_event"], event_confidence=node["event_confidence"],
            species_overlap_count=node["species_overlap_count"], mapping_conflict=node["mapping_conflict"],
            branch_support=node["branch_support"], lca_descendant_genes=len(node["genes"])))
    for ref in evidence:
        check(ref)
    with ledger.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(traced[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(traced)
    result = dict(schema="native_qfo_swiss_reconciliation_localization_v1", transition=transition_ref,
        source=record(__file__), evidence=evidence, reconstruction=reconstruction,
        selected_source_families=len(members), annotation_rows_read=annotation_rows,
        selected_annotation_nodes=sum(len(n) for n in nodes.values()), native_pair_rows_checked=native_rows,
        changed_pairs_traced=len(traced), summary=dict(summary), pair_ledger=record(ledger),
        observed_node_annotations=nodes_ref, node_annotations_previously_inventoried=False,
        new_scoring_or_admission=False, uncertainty_admitted=False, scientific_timings_admitted=False,
        independent_confirmation=False, publication_ready=False,
        limitations=[
            "Retrospective trace of every changed P0/C0 SwissTrees pair, not all reference or submitted relations.",
            "Reconstructed candidate hash matches R0 bytes; this does not prove identical search hit histories.",
            "Original RootHOG/manifest/prediction files were admitted; node annotations are newly observed current-byte evidence.",
            "Annotation topology, descendant sets and frozen pair-event logic checked without rerunning tree inference/reconciliation.",
            "Species overlap explains the software exclusion, not true duplication history or inferred-tree correctness.",
            "Serialized confidence/support are model annotations, not calibrated orthology confidence or independent validation.",
            "Initial HMM search remains on and R changes prediction semantics; no total-HMM or selected-default superiority claim.",
            "Failed R1 inference timing remains ineligible; diagnostics run on shared host with unknown contention effects.",
            "Further direct search/candidate evidence, tree-error analysis and broader uncertainty remain separate requirements."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("transitions", "output", "ledger"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--transitions-sha256", required=True)
    args = parser.parse_args()
    result = run(args.transitions, args.transitions_sha256, args.output, args.ledger)
    print(json.dumps({k: result[k] for k in ("reconstruction", "selected_source_families",
        "annotation_rows_read", "selected_annotation_nodes", "native_pair_rows_checked", "changed_pairs_traced", "summary")}))
