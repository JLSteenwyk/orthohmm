"""Independent saved-Newick topology check of every localized native exclusion."""

import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import re
import sys

from Bio import Phylo
import Bio

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.readback_native_qfo_swiss_transitions import record, require


def checked_tree(directory, family, species_sha, evidence):
    require(re.fullmatch(r"Family[0-9]{7}", family) is not None, "Invalid source family ID")
    checkpoint_ref = record(directory / "checkpoints" / (family + ".json"))
    checkpoint = json.loads(Path(checkpoint_ref["path"]).read_text())
    require(checkpoint["schema_version"] == 2 and checkpoint["status"] == "complete"
            and checkpoint["family_id"] == family and checkpoint["species_tree_sha256"] == species_sha,
            "Changed saved-tree checkpoint identity")
    genes = checkpoint["genes"]
    require(genes and len(set(genes)) == len(genes), "Duplicate or empty checkpoint genes")
    evidence.append(checkpoint_ref)
    trees = {}
    for suffix, key in (("raw", "raw_tree_sha256"), ("rooted", "rooted_tree_sha256"),
                        ("reconciled", "annotated_tree_sha256")):
        ref = record(directory / "gene_trees" / (family + "." + suffix + ".nwk"))
        require(ref["sha256"] == checkpoint[key], "Saved tree differs from checkpoint")
        evidence.append(ref)
        if suffix != "raw":
            tree = Phylo.read(ref["path"], "newick")
            leaves = [leaf.name for leaf in tree.get_terminals()]
            require(len(set(leaves)) == len(leaves) and set(leaves) == set(genes),
                    "Saved tree leaf universe differs")
            trees[suffix] = tree
    def topology(tree):
        return {frozenset(leaf.name for leaf in node.get_terminals()) for node in tree.find_clades()}
    require(topology(trees["rooted"]) == topology(trees["reconciled"]), "Annotation changed rooted topology")
    return trees["reconciled"]


def verify(report_path, report_sha):
    report_ref = record(report_path)
    require(report_ref["sha256"] == report_sha, "Changed localization report")
    report = json.loads(Path(report_path).read_text())
    require(report["schema"] == "native_qfo_swiss_reconciliation_localization_v1"
            and report["source"] == record(Path(__file__).with_name("trace_native_qfo_swiss_reconciliation.py"))
            and report["node_annotations_previously_inventoried"] is False
            and all(report[k] is False for k in ("new_scoring_or_admission", "uncertainty_admitted",
                "scientific_timings_admitted", "independent_confirmation", "publication_ready")),
            "Wrong localization source/scope")
    evidence = [report_ref, report["pair_ledger"], report["transition"], report["source"]]
    for ref in evidence:
        require(record(ref["path"]) == ref, "Changed direct localization evidence")
    transition = json.loads(Path(report["transition"]["path"]).read_text())
    before_ref = transition["changed_relations_ledger"]
    require(record(before_ref["path"]) == before_ref, "Changed original transition ledger")
    evidence.append(before_ref)
    with open(report["pair_ledger"]["path"], newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    with open(before_ref["path"], newline="") as stream:
        original = list(csv.DictReader(stream, delimiter="\t"))
    require([{k: row[k] for k in ("family", "protein_a", "protein_b", "before", "after")} for row in rows] == original
            and len(rows) == report["changed_pairs_traced"], "Incomplete or changed localized-pair ledger")
    directory = Path(report["observed_node_annotations"]["path"]).parent
    manifest_refs = [r for r in report["evidence"] if r["path"] == str(directory / "provenance_manifest.json")]
    require(len(manifest_refs) == 1 and record(manifest_refs[0]["path"]) == manifest_refs[0], "Changed native manifest")
    evidence.extend(manifest_refs)
    manifest = json.loads(Path(manifest_refs[0]["path"]).read_text())
    families = defaultdict(list)
    for row in rows:
        families[row["source_family"]].append(row)
    require(len(families) == report["selected_source_families"], "Changed selected family inventory")
    summary, nodes, total_leaves = Counter(), set(), 0
    for family, selected in sorted(families.items()):
        tree = checked_tree(directory, family, manifest["species_tree_sha256"], evidence)
        total_leaves += len(tree.get_terminals())
        leaves = {leaf.name: leaf for leaf in tree.get_terminals()}
        for row in selected:
            require(row["gene_a"] in leaves and row["gene_b"] in leaves, "Changed pair/tree membership")
            node = tree.common_ancestor(leaves[row["gene_a"]], leaves[row["gene_b"]])
            require(node.name and node.name.startswith(row["lca_node"] + "|D@")
                    and row["pair_event"] == "duplication"
                    and len(node.get_terminals()) == int(row["lca_descendant_genes"]),
                    "Saved Newick LCA differs from observed-node localization")
            same = row["root_hog_a"] == row["root_hog_b"]
            require(row["same_root_hog"] == str(same) and row["before"] in ("TP", "FP")
                    and row["event_confidence"] in ("high", "medium"), "Changed root/status/confidence")
            summary[row["before"] + ("_same_root" if same else "_different_root")] += 1
            summary[row["before"] + "_" + row["event_confidence"]] += 1
            nodes.add((family, row["lca_node"]))
    require(dict(summary) == report["summary"], "Saved-tree readback summary differs")
    for ref in evidence:
        require(record(ref["path"]) == ref, "Saved-tree evidence changed during readback")
    return dict(schema="native_qfo_swiss_reconciliation_newick_readback_v1", report=report_ref,
        source=record(__file__), checked_inputs=evidence, biopython_version=Bio.__version__,
        source_families_checked=len(families), tree_leaves_checked=total_leaves,
        distinct_exclusion_lcas_checked=len(nodes), changed_pairs_checked=len(rows), summary=dict(summary),
        original_node_annotation_admission_established=False, new_scoring_or_admission=False,
        uncertainty_admitted=False, independent_confirmation=False, publication_ready=False,
        limitations=[
            "Independent Biopython/Newick LCA/topology readback; same saved native data, not biological replication.",
            "Checkpoint/raw/rooted/annotated trees are newly observed and digest-consistent, not originally inventoried evidence.",
            "Does not rerun tree rooting/reconciliation or establish true duplication history or tree correctness.",
            "Direct localization/ledger/selected-tree checks, not repeated whole-study or transitive admission validation."])


if __name__ == "__main__":
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
    print(json.dumps({k: result[k] for k in ("source_families_checked", "tree_leaves_checked",
        "distinct_exclusion_lcas_checked", "changed_pairs_checked", "summary")}))
