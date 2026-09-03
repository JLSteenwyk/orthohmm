#!/usr/bin/env python3
"""Independent synthetic holdouts for OrthoHMM phylogeny-aware inference."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from itertools import combinations
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys
import time

from orthohmm.phylogeny import PhylogenyConfig
from orthohmm.phylogeny_pipeline import run_phylogeny_stage


AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
SPECIES = ("A", "B", "C", "D", "E", "F")
TRUE_SPECIES_TREE = "[&R] (((A,B),(C,D)),(E,F));"


@dataclass
class SyntheticDataset:
    input_directory: Path
    cluster_text: str
    species_tree: Path
    truth_groups: tuple[tuple[str, ...], ...]
    truth_by_family: dict[str, tuple[tuple[str, ...], ...]]
    scenarios: dict[str, str]


def random_sequence(rng: random.Random, length: int = 180) -> str:
    return "".join(rng.choice(AMINO_ACIDS) for _ in range(length))


def mutate(sequence: str, rate: float, rng: random.Random) -> str:
    residues = []
    for residue in sequence:
        if rng.random() < rate:
            choices = AMINO_ACIDS.replace(residue, "")
            residues.append(rng.choice(choices))
        else:
            residues.append(residue)
    return "".join(residues)


def evolve_species(sequence: str, rng: random.Random, scale: float = 1.0):
    abcd = mutate(sequence, 0.06 * scale, rng)
    clades = {
        "AB": mutate(abcd, 0.06 * scale, rng),
        "CD": mutate(abcd, 0.08 * scale, rng),
        "EF": mutate(sequence, 0.13 * scale, rng),
    }
    return {
        "A": mutate(clades["AB"], 0.04 * scale, rng),
        "B": mutate(clades["AB"], 0.05 * scale, rng),
        "C": mutate(clades["CD"], 0.05 * scale, rng),
        "D": mutate(clades["CD"], 0.06 * scale, rng),
        "E": mutate(clades["EF"], 0.07 * scale, rng),
        "F": mutate(clades["EF"], 0.09 * scale, rng),
    }


def build_dataset(root: Path, seed: int, marker_families: int = 24):
    rng = random.Random(seed)
    input_directory = root / "inputs"
    input_directory.mkdir(parents=True, exist_ok=True)
    records = {species: [] for species in SPECIES}
    clusters = []
    truth_groups = []
    truth_by_family = {}
    scenarios = {}

    def add_family(
        scenario: str,
        lineages: list[dict[str, str]],
        additions: dict[str, list[tuple[str, str]]] | None = None,
    ):
        family_index = len(clusters)
        family_id = f"Family{family_index:07d}"
        candidate = []
        family_truth = []
        for lineage_index, species_sequences in enumerate(lineages):
            group = []
            for species, sequence in sorted(species_sequences.items()):
                gene = f"{family_id}_L{lineage_index}_{species}"
                records[species].append((gene, sequence))
                candidate.append(gene)
                group.append(gene)
            family_truth.append(tuple(sorted(group)))
        for species, extra_records in (additions or {}).items():
            for gene_suffix, sequence in extra_records:
                gene = f"{family_id}_{gene_suffix}_{species}"
                records[species].append((gene, sequence))
                candidate.append(gene)
                family_truth[0] = tuple(sorted((*family_truth[0], gene)))
        clusters.append(tuple(sorted(candidate)))
        truth_groups.extend(family_truth)
        truth_by_family[family_id] = tuple(family_truth)
        scenarios[family_id] = scenario

    for marker_index in range(marker_families):
        leaves = evolve_species(random_sequence(rng), rng, scale=0.7)
        add_family(f"marker_{marker_index:02d}", [leaves])

    ancestral = random_sequence(rng)
    add_family(
        "ancient_duplication",
        [
            evolve_species(ancestral, rng),
            evolve_species(mutate(ancestral, 0.32, rng), rng),
        ],
    )

    ancestral = random_sequence(rng)
    lineage_one = evolve_species(ancestral, rng)
    lineage_two = evolve_species(mutate(ancestral, 0.30, rng), rng)
    add_family(
        "ancient_duplication_with_differential_loss",
        [
            {species: lineage_one[species] for species in ("A", "C", "E", "F")},
            {species: lineage_two[species] for species in ("B", "C", "D", "F")},
        ],
    )

    recent = evolve_species(random_sequence(rng), rng)
    add_family(
        "recent_duplication",
        [recent],
        additions={
            "A": [
                ("recent2", mutate(recent["A"], 0.03, rng)),
                ("recent3", mutate(recent["A"], 0.05, rng)),
            ]
        },
    )

    expanded = evolve_species(random_sequence(rng), rng)
    add_family(
        "lineage_specific_expansion",
        [expanded],
        additions={
            "E": [("expansion2", mutate(expanded["E"], 0.04, rng))],
            "F": [
                ("expansion2", mutate(expanded["F"], 0.04, rng)),
                ("expansion3", mutate(expanded["F"], 0.06, rng)),
            ],
        },
    )

    ancestral = random_sequence(rng, length=210)
    challenging_one = evolve_species(ancestral, rng, scale=1.2)
    challenging_two = evolve_species(mutate(ancestral, 0.35, rng), rng, scale=1.2)
    challenging_one["C"] = challenging_one["C"][45:145]
    challenging_one["B"] += random_sequence(rng, length=70)
    challenging_two["F"] = mutate(challenging_two["F"], 0.42, rng)
    biased = list(challenging_two["E"])
    for index in range(0, len(biased), 3):
        biased[index] = "A" if index % 2 else "G"
    challenging_two["E"] = "".join(biased)
    add_family(
        "ancient_duplication_remote_fragment_multidomain_composition_bias",
        [challenging_one, challenging_two],
    )

    for species in SPECIES:
        text = "".join(
            f">{gene}\n{sequence}\n" for gene, sequence in records[species]
        )
        (input_directory / f"{species}.faa").write_text(text, encoding="utf-8")
    cluster_text = "".join(" ".join(cluster) + "\n" for cluster in clusters)
    species_tree = root / "true_species_tree.nwk"
    species_tree.write_text(TRUE_SPECIES_TREE + "\n", encoding="utf-8")
    return SyntheticDataset(
        input_directory=input_directory,
        cluster_text=cluster_text,
        species_tree=species_tree,
        truth_groups=tuple(truth_groups),
        truth_by_family=truth_by_family,
        scenarios=scenarios,
    )


def pairs_from_groups(groups, cross_species_only=False):
    pairs = set()
    for group in groups:
        for left, right in combinations(sorted(group), 2):
            if cross_species_only and left.rsplit("_", 1)[-1] == right.rsplit("_", 1)[-1]:
                continue
            pairs.add((left, right))
    return pairs


def score_pairs(predicted, expected):
    true_positive = len(predicted & expected)
    precision = true_positive / len(predicted) if predicted else 0.0
    recall = true_positive / len(expected) if expected else 0.0
    f_score = (
        2.0 * precision * recall / (precision + recall)
        if precision + recall
        else 0.0
    )
    return {
        "precision": precision,
        "recall": recall,
        "f_score": f_score,
        "true_positive_pairs": true_positive,
        "predicted_pairs": len(predicted),
        "expected_pairs": len(expected),
    }


def read_root_groups(path: Path):
    groups_by_family = {}
    with path.open(encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            _, family, genes = line.rstrip("\n").split("\t")
            groups_by_family.setdefault(family, []).append(tuple(genes.split(",")))
    return groups_by_family


def read_native_pairs(path: Path):
    pairs = set()
    with path.open(encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            gene_a, _, gene_b, _ = line.split()
            pairs.add(tuple(sorted((gene_a, gene_b))))
    return pairs


def tree_distance(inferred_path: Path):
    import dendropy
    from dendropy.calculate import treecompare

    namespace = dendropy.TaxonNamespace(SPECIES)
    expected = dendropy.Tree.get(
        data=TRUE_SPECIES_TREE,
        schema="newick",
        taxon_namespace=namespace,
        preserve_underscores=True,
        rooting="default-rooted",
    )
    inferred = dendropy.Tree.get(
        path=str(inferred_path),
        schema="newick",
        taxon_namespace=namespace,
        preserve_underscores=True,
        rooting="default-rooted",
    )
    expected.encode_bipartitions()
    inferred.encode_bipartitions()
    distance = treecompare.symmetric_difference(expected, inferred)
    return {
        "rooted_symmetric_difference": int(distance),
        "exact_rooted_topology": distance == 0,
    }


def evaluate_mode(dataset: SyntheticDataset, output_directory: Path):
    root_groups_by_family = read_root_groups(
        output_directory / "orthohmm_phylogeny" / "orthohmm_root_hogs.tsv"
    )
    predicted_groups = tuple(
        group for groups in root_groups_by_family.values() for group in groups
    )
    partition = score_pairs(
        pairs_from_groups(predicted_groups),
        pairs_from_groups(dataset.truth_groups),
    )
    native_pairs = score_pairs(
        read_native_pairs(
            output_directory
            / "orthohmm_phylogeny"
            / "orthohmm_pairwise_orthologs.tsv"
        ),
        pairs_from_groups(dataset.truth_groups, cross_species_only=True),
    )
    predicted_sets = {frozenset(group) for group in predicted_groups}
    partition["exact_groups"] = sum(
        frozenset(group) in predicted_sets for group in dataset.truth_groups
    )
    partition["total_groups"] = len(dataset.truth_groups)
    by_scenario = {}
    for family_id, truth in dataset.truth_by_family.items():
        predicted = root_groups_by_family.get(family_id, [])
        by_scenario[dataset.scenarios[family_id]] = score_pairs(
            pairs_from_groups(predicted), pairs_from_groups(truth)
        )
    return {
        "root_groups": partition,
        "native_ortholog_pairs": native_pairs,
        "by_scenario": by_scenario,
    }


def run_mode(
    dataset,
    root,
    mode,
    aligner,
    tree_builder,
    cpu,
    root_rule,
    pair_rule,
):
    output = root / mode
    working = output / "orthohmm_working_res"
    working.mkdir(parents=True, exist_ok=True)
    (working / "orthohmm_edges_clustered.txt").write_text(
        dataset.cluster_text, encoding="utf-8"
    )
    config = PhylogenyConfig(
        mode="reconcile",
        species_tree_mode=mode,
        species_tree=str(dataset.species_tree) if mode == "supplied" else None,
        aligner=aligner,
        tree_builder=tree_builder,
        root_duplication_rule=root_rule,
        pair_orthology_rule=pair_rule,
    )
    started = time.perf_counter()
    stage = run_phylogeny_stage(
        str(dataset.input_directory),
        str(output),
        [f"{species}.faa" for species in SPECIES],
        config,
        cpu,
    )
    evaluation = evaluate_mode(dataset, output)
    evaluation["stage"] = stage.__dict__
    evaluation["wall_s"] = time.perf_counter() - started
    if mode == "infer":
        evaluation["species_tree"] = tree_distance(
            output / "orthohmm_phylogeny" / "species_tree.rooted.nwk"
        )
    return evaluation


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("output_json", type=Path)
    parser.add_argument("--work-directory", type=Path, required=True)
    parser.add_argument("--aligner", required=True)
    parser.add_argument("--tree-builder", required=True)
    parser.add_argument("--cpu", type=int, default=4)
    parser.add_argument("--seeds", type=int, nargs="+", default=[101, 211, 307])
    parser.add_argument(
        "--root-rule",
        choices=(
            "supported_children", "confidence", "species_overlap", "mapped_event"
        ),
        default="supported_children",
    )
    parser.add_argument(
        "--pair-rule",
        choices=(
            "lca", "positive_paralogy", "supported_paralogy",
            "sparse_overlap", "depth_two_overlap", "depth_two_closure",
        ),
        default="lca",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    args.work_directory.mkdir(parents=True, exist_ok=True)
    results = []
    for seed in args.seeds:
        replicate_root = args.work_directory / f"seed_{seed}"
        if replicate_root.exists():
            shutil.rmtree(replicate_root)
        dataset = build_dataset(replicate_root / "dataset", seed)
        results.append({
            "seed": seed,
            "supplied": run_mode(
                dataset, replicate_root, "supplied",
                args.aligner, args.tree_builder, args.cpu,
                args.root_rule, args.pair_rule,
            ),
            "infer": run_mode(
                dataset, replicate_root, "infer",
                args.aligner, args.tree_builder, args.cpu,
                args.root_rule, args.pair_rule,
            ),
        })
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        check=False,
    )
    payload = {
        "schema_version": 1,
        "description": "Independent synthetic phylogeny holdouts",
        "source_commit": completed.stdout.strip() or None,
        "command": [sys.executable, *sys.argv],
        "root_duplication_rule": args.root_rule,
        "pair_orthology_rule": args.pair_rule,
        "scenarios": [
            "ancient duplication",
            "recent duplication",
            "differential gene loss",
            "remote homolog",
            "fragment",
            "multidomain protein",
            "lineage-specific expansion",
            "composition bias",
        ],
        "tools": {
            "aligner": str(Path(args.aligner).resolve()),
            "tree_builder": str(Path(args.tree_builder).resolve()),
        },
        "results": results,
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(args.output_json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
