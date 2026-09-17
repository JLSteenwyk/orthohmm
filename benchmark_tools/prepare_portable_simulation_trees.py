"""Adapt frozen rooted trees to plain Newick accepted by both native parsers."""

import argparse
from copy import deepcopy
import io
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import Phylo
from benchmark_tools.prepare_species_tree_robustness import topology
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen

SOURCE_SHA = "0caffab16f73019fabafe1fcfaac8c2de90960c5f1317803bf83308a6d2bf914"


def branch_map(tree):
    return {tuple(sorted(t.name for t in node.get_terminals())): node.branch_length
            for node in tree.find_clades()}


def portable_text(tree):
    topology(tree)
    copied = deepcopy(tree)
    for node in copied.find_clades():
        node.comment = None
    handle = io.StringIO()
    Phylo.write(copied, handle, "newick", format_branch_length="%.17g")
    text = handle.getvalue()
    reread = Phylo.read(io.StringIO(text), "newick")
    if topology(reread) != topology(tree) or branch_map(reread) != branch_map(tree):
        raise ValueError("Plain-Newick serialization changed topology or branch lengths")
    return text


def prepare(root, output):
    if output.exists():
        raise FileExistsError(output)
    source = root / "benchmark_tools/results/simulation_tree_controls_prepared_20260917.json"
    frozen = read_frozen(source, SOURCE_SHA)
    output.mkdir(parents=True)
    rows = []
    for dataset in frozen["datasets"]:
        for variant in dataset["variants"]:
            check(variant["tree"])
            original = Phylo.read(variant["tree"]["path"], "newick")
            path = output / f"{dataset['condition']}_{dataset['seed']}_{variant['label']}.nwk"
            path.write_text(portable_text(original))
            check(variant["tree"])
            rows.append({"condition": dataset["condition"], "seed": dataset["seed"],
                         "variant": variant["label"], "taxa": dataset["taxa"],
                         "source_tree": variant["tree"], "tree": record(path)})
    read_frozen(source, SOURCE_SHA)
    if len(rows) != 210:
        raise ValueError("Incomplete tree panel")
    result = {"status": "portable_trees_prepared_pending_native_parser_checks",
              "source": record(__file__), "source_manifest": record(source), "trees": rows,
              "accuracy_evaluated": False,
              "reason": "OrthoFinder 3.1.5 native parser reads the leading [&R] tree as NoName; plain Newick keeps the binary root without that annotation."}
    (output / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.output.resolve())
