"""Prepare a deterministic label-blind rooted species-tree perturbation panel."""

import argparse
from copy import deepcopy
import hashlib
import io
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import Phylo
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.validate_factorial_native import validate


def topology(tree):
    taxa = [n.name for n in tree.get_terminals()]
    if len(taxa) < 4 or None in taxa or len(set(taxa)) != len(taxa):
        raise ValueError("Need at least four uniquely named taxa")
    if any(len(n.clades) != 2 for n in tree.get_nonterminals()):
        raise ValueError("Require a fully bifurcating rooted tree")
    if any(n.branch_length is not None and (not math.isfinite(n.branch_length) or n.branch_length < 0)
           for n in tree.find_clades()):
        raise ValueError("Invalid branch length")
    return frozenset(tuple(sorted(t.name for t in node.get_terminals()))
                     for node in tree.get_nonterminals() if node is not tree.root)


def digest(key):
    return hashlib.sha256(json.dumps(sorted(key), separators=(",", ":")).encode()).hexdigest()


def neighbors(tree):
    """Rooted NNI: swap either grandchild with its parent's sibling subtree."""
    topology(tree)
    results = {}
    for parent_index, parent in enumerate(tree.get_nonterminals(order="preorder")):
        for child_index, child in enumerate(parent.clades):
            if child.is_terminal():
                continue
            sibling_index = 1 - child_index
            for grandchild_index in (0, 1):
                new = deepcopy(tree)
                new_parent = new.get_nonterminals(order="preorder")[parent_index]
                new_child = new_parent.clades[child_index]
                new_parent.clades[sibling_index], new_child.clades[grandchild_index] = (
                    new_child.clades[grandchild_index], new_parent.clades[sibling_index])
                key = topology(new)
                results.setdefault(key, new)
    return results


def panel(tree):
    baseline = topology(tree)
    single = neighbors(tree)
    single = {key: value for key, value in single.items() if len(key ^ baseline) == 2}
    selected_single = sorted(single, key=digest)[:3]
    double = {}
    for first_key in selected_single:
        for key, value in neighbors(single[first_key]).items():
            if len(key ^ baseline) == 4:
                double.setdefault(key, (value, first_key))
    selected_double = sorted(double, key=digest)[:3]
    if len(selected_single) != 3 or len(selected_double) != 3:
        raise ValueError("Tree cannot provide the planned three variants at each distance")
    records = [("supplied_control", deepcopy(tree), 0, None)]
    records += [(f"nni1_{i}", single[key], 2, digest(baseline)) for i, key in enumerate(selected_single)]
    records += [(f"nni2_{i}", double[key][0], 4, digest(double[key][1])) for i, key in enumerate(selected_double)]
    return records


def serialize(tree):
    tree = deepcopy(tree)
    for node in tree.get_nonterminals():
        node.name, node.confidence, node.comment = None, None, None
    for node in tree.get_terminals():
        node.comment = None
    handle = io.StringIO()
    Phylo.write(tree, handle, "newick", format_branch_length="%.17g")
    return "[&R] " + handle.getvalue().strip() + "\n"


def prepare(root, output):
    root, output = root.resolve(), output.resolve()
    if output.exists():
        raise FileExistsError(output)
    native = validate(root, 3)
    source = native["species_tree"]
    tree = Phylo.read(source["path"], "newick")
    variants = panel(tree)
    baseline = topology(tree)
    taxa = sorted(n.name for n in tree.get_terminals())
    output.mkdir(parents=True)
    records = []
    for label, candidate, distance, parent in variants:
        path = output / (label + ".nwk")
        with path.open("x") as handle:
            handle.write(serialize(candidate))
        reread = Phylo.read(path, "newick")
        if topology(reread) != topology(candidate) or sorted(n.name for n in reread.get_terminals()) != taxa:
            raise ValueError("Serialization changed tree topology or taxa")
        if len(topology(reread) ^ baseline) != distance:
            raise ValueError("Wrong rooted clade distance")
        records.append({"label": label, "tree": file_provenance(path), "topology_sha256": digest(topology(reread)),
                        "rooted_rf_clade_distance": distance, "normalized_rooted_rf": distance / (2 * len(baseline)),
                        "parent_topology_sha256": parent})
    verify_file(Path(source["path"]), source)
    report = {"status": "trees_prepared_not_reconciled", "accuracy_evaluated": False,
              "source_tree": source, "native_admission": native, "taxa": taxa, "variants": records,
              "generator": file_provenance(Path(__file__)), "selection": "First three SHA256-ranked distinct topologies per level; no labels",
              "limitations": ["Rooted clade perturbations are controlled diagnostics, not sampled biological tree uncertainty.",
                  "Branch lengths travel with moved subtrees; distances need not retain an evolutionary interpretation.",
                  "Internal support labels removed; leaf names and branch values retained with explicit rooted Newick marker.",
                  "Two-swap candidates drawn from the three selected one-swap trees; not an exhaustive radius-two distribution.",
                  "Supplied-tree control must reproduce inferred-tree output before perturbation effects are interpreted."]}
    (output / "manifest.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    prepare(args.root, args.output)


if __name__ == "__main__":
    main()
