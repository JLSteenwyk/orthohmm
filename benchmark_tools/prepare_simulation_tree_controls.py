"""Prepare generating-tree and deterministic topology controls for all frozen simulation datasets."""

import argparse
from copy import deepcopy
import itertools
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import Phylo
from benchmark_tools.prepare_species_tree_robustness import panel, serialize, topology, digest
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen, verify_inputs

GENERATION_SHA = "806aa1e5f6976c323ff2f6641dd88e25264e7417749e7d77565294666aee776b"
CONDITIONS = ("baseline", "divergent", "turnover", "divergent_turnover", "missing20", "uneven_taxa", "taxon_count_control")
METHODS = ("orthohmm_satellite_v2", "orthofinder_3_1_5_full")


def restrict_tree(tree, retained):
    topology(tree)
    original = {tip.name for tip in tree.get_terminals()}
    if len(retained) < 4 or not set(retained) <= original or len(set(retained)) != len(retained):
        raise ValueError("Require at least four unique retained taxa from the generating tree")
    result = deepcopy(tree)
    for name in sorted(original - set(retained)):
        result.prune(name)
    topology(result)
    if {tip.name for tip in result.get_terminals()} != set(retained):
        raise ValueError("Pruning changed intended taxa")
    for a, b in itertools.combinations(sorted(retained), 2):
        if abs(tree.distance(a, b) - result.distance(a, b)) > 1e-10:
            raise ValueError("Pruning changed retained pairwise tree distances")
    return result


def controls(tree):
    variants = {label: (candidate, distance, parent) for label, candidate, distance, parent in panel(tree)}
    return [(name, *variants[key]) for name, key in
            (("generating", "supplied_control"), ("nni1", "nni1_0"), ("nni2", "nni2_0"))]


def prepare(root, output):
    if output.exists():
        raise FileExistsError(output)
    manifest_path = root / "benchmark_tools/results/publication_variable_simulation_manifest_20260916.json"
    generation = read_frozen(manifest_path, GENERATION_SHA)
    expected = {(condition, seed) for condition in CONDITIONS for seed in range(20261101, 20261111)}
    datasets = generation["datasets"]
    if len(datasets) != 70 or {(d["condition"], d["seed"]) for d in datasets} != expected:
        raise ValueError("Changed frozen simulation inventory")
    simulation_panel = root / "benchmarks/work/publication_variable_simulation_panel_v2"
    runs = {r["label"]: r for r in generation["simulation_runs"]}
    output.mkdir(parents=True)
    rows, planned = [], []
    for dataset in sorted(datasets, key=lambda d: (d["seed"], CONDITIONS.index(d["condition"]))):
        evidence = verify_inputs(dataset, generation, simulation_panel, GENERATION_SHA)
        if evidence["status"] != "ready":
            rows.append({"condition": dataset["condition"], "seed": dataset["seed"], "status": "inapplicable", "reason": evidence})
            continue
        source_tree = Path(runs[dataset["parent"]]["native_output"]) / "T/ExtantTree.nwk"
        tree_record = record(source_tree)
        names = sorted(Path(item["absolute_path"]).stem for item in evidence["inputs"])
        tree = restrict_tree(Phylo.read(source_tree, "newick"), names)
        baseline = topology(tree)
        directory = output / f"{dataset['condition']}_{dataset['seed']}"
        directory.mkdir()
        variants = []
        for label, candidate, distance, parent in controls(tree):
            path = directory / (label + ".nwk")
            path.write_text(serialize(candidate))
            reread = Phylo.read(path, "newick")
            if ({tip.name for tip in reread.get_terminals()} != set(names) or topology(reread) != topology(candidate)
                    or len(topology(reread) ^ baseline) != distance):
                raise ValueError("Serialized tree control changed taxa/topology/distance")
            variants.append({"label": label, "tree": record(path), "rooted_clade_distance": distance,
                             "topology_sha256": digest(topology(reread)), "parent_topology_sha256": parent})
            for method in METHODS:
                planned.append({"index": len(planned), "condition": dataset["condition"], "seed": dataset["seed"],
                                "tree_control": label, "method": method, "tree": record(path)})
        check(tree_record)
        rows.append({"condition": dataset["condition"], "seed": dataset["seed"], "status": "trees_prepared_not_run",
                     "parent": dataset["parent"], "generating_tree": tree_record, "taxa": names,
                     "input_evidence": evidence, "variants": variants})
        (output / "progress.json").write_text(json.dumps({"datasets_prepared": len(rows), "planned_runs": len(planned)}) + "\n")
    read_frozen(manifest_path, GENERATION_SHA)
    report = {"status": "simulation_tree_controls_prepared_unrun", "accuracy_evaluated": False,
              "source": record(__file__), "tree_helper": record(Path(__file__).with_name("prepare_species_tree_robustness.py")),
              "generation_manifest": record(manifest_path), "datasets": rows, "planned_runs": planned,
              "selection": "All frozen condition/seed cells; generating tree plus SHA256-first one-NNI and two-NNI controls.",
              "limitations": ["Generating trees are oracle inputs, not information available in ordinary inference.",
                  "Tree topologies are development-exposed stress controls, not posterior samples or independent validation.",
                  "Lengths travel with NNI subtrees; perturbed trees need not be ultrametric or preserve distances.",
                  "Tree preparation is not method execution, accuracy evaluation, or completed robustness evidence."]}
    (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.output.resolve())
