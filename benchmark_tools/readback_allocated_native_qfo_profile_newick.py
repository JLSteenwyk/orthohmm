"""Read saved Newick and root membership for every localized profile pair change."""

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import re

import Bio

from benchmark_tools import readback_native_qfo_swiss_reconciliation as saved
from benchmark_tools.readback_native_qfo_swiss_transitions import record, require

CELLS = ("p0_c0_r1", "p1_c0_r1")
EVENT_CODES = {"speciation": "S", "duplication": "D", "uncertain": "U"}
SCOPE_FLAGS = ("new_scoring_or_admission", "uncertainty_admitted",
               "scientific_timings_admitted", "independent_confirmation", "publication_ready")


def selected_roots(path, families):
    members = {family: set() for family in families}
    roots, seen = {}, set()
    with Path(path).open(newline="") as stream:
        rows = csv.reader(stream, delimiter="\t")
        require(next(rows) == ["root_hog", "source_family", "genes"], "Changed root header")
        for row in rows:
            require(len(row) == 3, "Invalid root row width")
            root, family, genes_text = row
            if family not in members:
                continue
            require(re.fullmatch(r"RootHOG[0-9]{7}", root) is not None and root not in seen,
                    "Invalid or duplicate selected root")
            seen.add(root)
            genes = genes_text.split(",")
            require(all(genes) and len(set(genes)) == len(genes)
                    and not any(gene in roots for gene in genes), "Duplicate or empty selected gene")
            members[family].update(genes)
            roots.update({gene: (family, root) for gene in genes})
    require(all(members.values()), "Missing selected root family")
    return members, roots


def check_point(point, trees, members, roots):
    leaves = {}
    for side in ("a", "b"):
        family, gene = point["source_family_" + side], point["gene_" + side]
        require(family in trees and roots.get(gene) == (family, point["root_hog_" + side]),
                "Pair differs from selected root membership")
        genes = sorted(members[family])
        digest = hashlib.sha256(("\n".join(genes) + "\n").encode()).hexdigest()
        require(type(point["candidate_size_" + side]) is int
                and len(genes) == point["candidate_size_" + side]
                and digest == point["candidate_members_" + side + "_sha256"], "Changed candidate membership")
        terminal = {leaf.name: leaf for leaf in trees[family].get_terminals()}
        require(set(terminal) == members[family] and gene in terminal, "Candidate/tree membership differs")
        leaves[side] = terminal[gene]
    same = point["source_family_a"] == point["source_family_b"]
    require(point["same_source_family"] is same and point["gene_a"] != point["gene_b"],
            "Changed pair family identity")
    if not same:
        require(point["lca"] is None and point["predicted"] is False,
                "Separated pair has invented LCA or prediction")
        return None
    lca = point["lca"]
    require(isinstance(lca, dict) and lca["pair_event"] in EVENT_CODES, "Missing or invalid pair LCA")
    node = trees[point["source_family_a"]].common_ancestor(leaves["a"], leaves["b"])
    prefix = lca["node_id"] + "|" + EVENT_CODES[lca["pair_event"]] + "@"
    require(node.name and node.name.startswith(prefix), "Saved Newick LCA differs from localization")
    require(point["predicted"] is (lca["pair_event"] != "duplication"), "Prediction differs from Newick pair event")
    return point["source_family_a"], lca["node_id"]


def verify(report_path, report_sha):
    report_ref = record(report_path)
    require(report_ref["sha256"] == report_sha, "Changed profile localization report")
    report = json.loads(Path(report_path).read_text())
    require(report["schema"] == "allocated_native_qfo_profile_pair_localization_v1"
            and report["source"] == record(Path(__file__).with_name("trace_allocated_native_qfo_profile.py"))
            and report["node_annotations_previously_inventoried"] is False
            and all(report[key] is False for key in SCOPE_FLAGS), "Changed localization source/scope")
    require([state["cell"] for state in report["states"]] == list(CELLS), "Changed profile cells")
    evidence = [report_ref, report["source"], report["readback"], *report["helpers"], *report["evidence"]]
    for ref in evidence:
        require(record(ref["path"]) == ref, "Changed direct localization evidence")
    helpers = [record(Path(__file__).with_name(name)) for name in (
        "trace_native_qfo_swiss_reconciliation.py", "trace_native_qfo_swiss_transitions.py")]
    require(report["helpers"] == helpers, "Changed localization helpers")
    pairs = report["changed_pairs"]
    require(pairs and len(pairs) == report["changed_pairs_traced"], "Incomplete changed pair inventory")
    keys = [(p["family"], p["protein_a"], p["protein_b"]) for p in pairs]
    require(len(set(keys)) == len(keys) and keys == sorted(keys), "Duplicate or unordered changed pairs")
    observed = Counter((p["family"], p["before_label"] + "->" + p["after_label"]) for p in pairs)
    expected = Counter()
    for family in report["comparison"]["families"]:
        for transition, count in family["transitions"].items():
            before, after = transition.split("->")
            if before != after and count:
                expected[(family["family"], transition)] = count
    require(observed == expected, "Pair inventory differs from complete transition comparison")
    state_trees, state_members, state_roots, counts = [], [], [], []
    for key, state in zip(("before", "after"), report["states"]):
        require(state["root_hogs"] in report["evidence"] and state["manifest"] in report["evidence"],
                "Unbound selected roots or manifest")
        directory = Path(state["root_hogs"]["path"]).parent
        require(Path(state["manifest"]["path"]) == directory / "provenance_manifest.json", "Wrong selected manifest")
        manifest = json.loads(Path(state["manifest"]["path"]).read_text())
        require(manifest["species_tree_sha256"] == state["species_tree_sha256"]
                and manifest["pair_orthology_rule"] == "positive_paralogy"
                and manifest["membership_reconciliation"] is None, "Changed species/pair rule")
        families = {p[key]["source_family_" + side] for p in pairs for side in ("a", "b")}
        require(len(families) == state["selected_source_families"], "Changed selected family count")
        members, roots = selected_roots(state["root_hogs"]["path"], families)
        trees = {family: saved.checked_tree(directory, family, state["species_tree_sha256"], evidence)
                 for family in sorted(families)}
        require(all({leaf.name for leaf in trees[f].get_terminals()} == members[f] for f in families),
                "Saved tree differs from root candidate members")
        state_trees.append(trees)
        state_members.append(members)
        state_roots.append(roots)
        counts.append(dict(cell=state["cell"], source_families_checked=len(families),
                           tree_leaves_checked=sum(len(tree.get_terminals()) for tree in trees.values())))
    summary, lcas = Counter(), set()
    for pair in pairs:
        require((pair["before"]["gene_a"], pair["before"]["gene_b"])
                == (pair["after"]["gene_a"], pair["after"]["gene_b"]), "Changed target gene identity")
        for index, key in enumerate(("before", "after")):
            node = check_point(pair[key], state_trees[index], state_members[index], state_roots[index])
            if node is not None:
                lcas.add((CELLS[index], *node))
            require(pair[key]["predicted"] is (pair[key + "_label"] in ("TP", "FP")), "Changed prediction/label")
        for side in ("a", "b"):
            unchanged = pair["before"]["candidate_members_" + side + "_sha256"] == pair["after"]["candidate_members_" + side + "_sha256"]
            require(pair["candidate_members_" + side + "_unchanged"] is unchanged, "Changed membership equality")
        require(pair["before"]["predicted"] is True and pair["after"]["predicted"] is False
                and (pair["before_label"] in ("TP", "FN")) == (pair["after_label"] in ("TP", "FN")),
                "Unexpected change in prediction or reference truth")
        category = "candidate_separation" if not pair["after"]["same_source_family"] else "positive_paralogy_exclusion"
        require(pair["localization"] == category, "Changed localization category")
        summary[category] += 1
    same_species = report["states"][0]["species_tree_sha256"] == report["states"][1]["species_tree_sha256"]
    require(report["species_tree_bytes_identical"] is same_species and dict(summary) == report["summary"],
            "Changed species-tree or localization summary")
    evidence.append(record(saved.__file__))
    for ref in evidence:
        require(record(ref["path"]) == ref, "Selected evidence changed during readback")
    return dict(schema="allocated_native_qfo_profile_newick_readback_v1", source=record(__file__),
        report=report_ref, checked_inputs=evidence, biopython_version=Bio.__version__, cells=list(CELLS),
        states=counts, source_families_checked=sum(c["source_families_checked"] for c in counts),
        tree_leaves_checked=sum(c["tree_leaves_checked"] for c in counts), changed_pairs_checked=len(pairs),
        distinct_lcas_checked=len(lcas), summary=dict(summary), species_tree_bytes_identical=same_species,
        original_node_annotation_admission_established=False, new_scoring_or_admission=False,
        uncertainty_admitted=False, scientific_timings_admitted=False, independent_confirmation=False,
        publication_ready=False, limitations=[
            "Biopython/Newick plus separate streamed root-membership parser, not the primary node-table LCA code.",
            "Saved trees/checkpoints newly observed and digest-consistent, not originally inventoried biological evidence.",
            "Checks software grouping/pair inclusion, not biological histories, edge causality or tree correctness.",
            "Checks annotated pair events, not the node-table confidence/overlap calculation or an independent inference.",
            "Same admitted outputs and scored SwissTrees contrast; no new scoring, bootstrap or timing admission."])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), "Output already exists")
    result = verify(args.report, args.report_sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({key: result[key] for key in ("source_families_checked", "tree_leaves_checked",
        "changed_pairs_checked", "distinct_lcas_checked", "summary", "species_tree_bytes_identical")}))


if __name__ == "__main__":
    main()
