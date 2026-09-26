"""Independent structural readback, not a reimplementation of reconciliation."""

import argparse
import csv
import json
from collections import Counter, defaultdict
from itertools import zip_longest
from pathlib import Path

import dendropy

from benchmark_tools.audit_installed_orthobench import read_root_hogs
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_installed_orthobench import fasta_ids


def rows(path, columns):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != columns:
            raise ValueError(f"Unexpected columns: {path}")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Malformed row: {path}")
            yield row


def tree_leaves(path, expected):
    trees = dendropy.TreeList.get(path=str(path), schema="newick", preserve_underscores=True)
    if len(trees) != 1:
        raise ValueError(f"Expected one tree: {path}")
    labels = [node.taxon.label if node.taxon else None for node in trees[0].leaf_node_iter()]
    if len(labels) != len(set(labels)) or set(labels) != set(expected):
        raise ValueError(f"Unexpected tree leaves: {path}")


def pair_readback(directory, gene_species, gene_family, expected_count):
    columns = ["gene_a", "species_a", "gene_b", "species_b"]
    plain = rows(directory / "orthohmm_pairwise_orthologs.tsv", columns)
    annotated = rows(directory / "orthohmm_pairwise_orthologs_confidence.tsv", columns + ["confidence"])
    previous = None
    counts = Counter()
    for left, right in zip_longest(plain, annotated):
        if left is None or right is None or left != {key: right[key] for key in columns}:
            raise ValueError("Pair/confidence tables differ")
        a, b = left["gene_a"], left["gene_b"]
        pair = (a, b)
        if (a not in gene_species or b not in gene_species or not a < b
                or (previous is not None and pair <= previous)
                or left["species_a"] != gene_species[a] or left["species_b"] != gene_species[b]
                or gene_species[a] == gene_species[b] or gene_family[a] != gene_family[b]
                or right["confidence"] not in {"high", "medium", "low"}):
            raise ValueError("Invalid, repeated, unsorted or foreign ortholog pair")
        previous = pair
        counts[right["confidence"]] += 1
    if sum(counts.values()) != expected_count:
        raise ValueError("Pair count differs from summary")
    return dict(counts)


def audit(directory, input_dir):
    watched = []

    def watch(path):
        item = record(path)
        watched.append(item)
        return item

    def read_json(name):
        path = directory / name
        watch(path)
        return json.loads(path.read_text())

    manifest = read_json("provenance_manifest.json")
    summary = read_json("reconciliation_summary.json")
    if manifest["results"] != summary:
        raise ValueError("Manifest and summary differ")
    gene_species, taxa, filenames = {}, set(), set()
    for item in manifest["input_proteomes"]:
        name, taxon = item["filename"], item["taxon"]
        if (Path(name).name != name or taxon != Path(name).stem
                or name in filenames or taxon in taxa):
            raise ValueError("Invalid or duplicate proteome identity")
        path = input_dir / name
        if watch(path)["sha256"] != item["sha256"]:
            raise ValueError("Input proteome hash differs")
        genes = fasta_ids([path])
        if not genes or genes & gene_species.keys():
            raise ValueError("Empty proteome or duplicate gene identity")
        gene_species.update(dict.fromkeys(genes, taxon))
        taxa.add(taxon)
        filenames.add(name)
    if sorted(taxa) != sorted(manifest["species_tree_taxa"]):
        raise ValueError("Species manifest mismatch")
    species_path = directory / "species_tree.rooted.nwk"
    species_record = watch(species_path)
    if species_record["sha256"] != manifest["species_tree_sha256"]:
        raise ValueError("Species tree hash differs")
    tree_leaves(species_path, taxa)
    roots = directory / "orthohmm_root_hogs.tsv"
    watch(roots)
    groups = read_root_hogs(roots, set(gene_species))
    families, gene_family = defaultdict(set), {}
    for row in rows(roots, ["root_hog", "source_family", "genes"]):
        genes = row["genes"].split(",")
        families[row["source_family"]].update(genes)
        gene_family.update(dict.fromkeys(genes, row["source_family"]))
    if len(groups) != summary["root_hogs"] or len(families) != summary["candidate_families"]:
        raise ValueError("Group/family count mismatch")
    checkpoints = sorted((directory / "checkpoints").glob("*.json"))
    expected_tree_names = set()
    for path in checkpoints:
        watch(path)
        checkpoint = json.loads(path.read_text())
        family = path.stem
        genes = checkpoint["genes"]
        if (checkpoint["status"] != "complete" or checkpoint["family_id"] != family
                or family not in families or len(genes) != len(set(genes))
                or set(genes) != families[family]
                or checkpoint["species_tree_sha256"] != species_record["sha256"]):
            raise ValueError("Checkpoint family differs from final family membership")
        for suffix, hash_key, expected in (
            ("raw", "raw_tree_sha256", {f"g{i:08d}" for i in range(len(genes))}),
            ("rooted", "rooted_tree_sha256", genes),
            ("reconciled", "annotated_tree_sha256", genes),
        ):
            tree_path = directory / "gene_trees" / f"{family}.{suffix}.nwk"
            expected_tree_names.add(tree_path.name)
            if watch(tree_path)["sha256"] != checkpoint[hash_key]:
                raise ValueError("Checkpoint tree hash differs")
            tree_leaves(tree_path, expected)
    actual_tree_names = {p.name for p in (directory / "gene_trees").iterdir()}
    if actual_tree_names != expected_tree_names:
        raise ValueError("Unexpected gene-tree directory inventory")
    if (len(checkpoints) != summary["reconciled_families"]
            or len(families) - len(checkpoints) != summary["bypassed_families"]):
        raise ValueError("Reconciled/bypassed family count mismatch")
    for name in ("orthohmm_pairwise_orthologs.tsv", "orthohmm_pairwise_orthologs_confidence.tsv"):
        watch(directory / name)
    confidence = pair_readback(directory, gene_species, gene_family, summary["ortholog_pairs"])
    for item in watched:
        check(item)
    return dict(status="phylogeny_structure_verified", genes=len(gene_species), species=len(taxa),
        families=len(families), root_groups=len(groups), reconciled_families=len(checkpoints),
        parsed_trees=1 + 3 * len(checkpoints), ortholog_pairs=sum(confidence.values()),
        pair_confidence_counts=confidence, checked_records=watched,
        reader=record(__file__), dendropy_version=dendropy.__version__,
        scientific_scores_admitted=False, reconciliation_semantics_recomputed=False,
        limitations=["Structural consistency does not prove biological correctness or reproduce event inference.",
                     "Native pairs are validated separately; they are not equated with root-group pairs.",
                     "Candidate sequences, alignment content, node events and pair completeness remain separate checks.",
                     "Caller must gate native completion, exact input inventory and frozen provenance."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = audit(args.directory.resolve(), args.input.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
