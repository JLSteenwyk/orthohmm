"""Strict native output adapters for the four frozen simulation methods."""

import csv
import itertools
from pathlib import Path

from benchmark_tools.orthofinder_to_pairwise import iter_pairs
from benchmark_tools.report_ygob_validation import read_checkpoint
from benchmark_tools.score_ygob_groups import read_predictions
from benchmark_tools.simulation_conditions import group_pairs


def unique_path(root, pattern):
    paths = sorted(root.glob(pattern))
    if len(paths) != 1:
        raise ValueError(f"Expected exactly one native artifact for {pattern}, found {len(paths)}")
    return paths[0]


def orthohmm_pairs(path, gene_species):
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["gene_a", "species_a", "gene_b", "species_b"]:
            raise ValueError("Unexpected OrthoHMM pair header")
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError("Malformed OrthoHMM pair row")
            for suffix in ("a", "b"):
                if gene_species.get(row["gene_" + suffix]) != row["species_" + suffix]:
                    raise ValueError("OrthoHMM pair species/ID mismatch")
            if row["species_a"] == row["species_b"]:
                raise ValueError("OrthoHMM pair is not cross-species")
            yield row["gene_a"], row["gene_b"]


def orthofinder_pairs(results, gene_species, species):
    """Default full output requires both orientations, including empty tables."""
    expected = set(itertools.permutations(species, 2))
    observed = {}
    files = sorted((results / "Orthologues").glob("Orthologues_*/*__v__*.tsv"))
    for path in files:
        with path.open() as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, None)
            if header is None or len(header) != 3 or header[0] != "Orthogroup":
                raise ValueError("Unexpected OrthoFinder table header")
            a, b = header[1:]
            if (a, b) not in expected or (a, b) in observed:
                raise ValueError("Duplicate or unexpected species-pair table")
            pairs = set()
            for row in reader:
                if len(row) != 3 or not row[0]:
                    raise ValueError("Malformed OrthoFinder relation row")
                left, right = [[gene.strip() for gene in cell.split(",") if gene.strip()] for cell in row[1:]]
                if not left or not right:
                    raise ValueError("Empty relation endpoint in OrthoFinder row")
                if any(gene_species.get(g) != a for g in left) or any(gene_species.get(g) != b for g in right):
                    raise ValueError("OrthoFinder pair species/ID mismatch")
                pairs.update(tuple(sorted(pair)) for pair in itertools.product(left, right))
            observed[(a, b)] = pairs
    if set(observed) != expected:
        raise ValueError("Incomplete default OrthoFinder species-pair table set")
    for (a, b), pairs in observed.items():
        if pairs != observed[(b, a)]:
            raise ValueError("OrthoFinder table orientations disagree")
    # Reuse the audited native converter, now with complete-table and ID checks.
    return list(iter_pairs(results))


def load_predictions(method, output, gene_species, species):
    if set(gene_species.values()) - set(species):
        raise ValueError("Input species mapping is inconsistent")
    if method == "orthohmm_high_sensitivity":
        path = output / "orthohmm_orthogroups.txt"
        groups = read_predictions(path, "named_groups")
        return list(group_pairs(groups.values(), gene_species)), [path]
    if method == "orthohmm_satellite_v2":
        path = output / "orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv"
        return list(orthohmm_pairs(path, gene_species)), [path]
    if method == "orthofinder_full":
        orthologues = unique_path(output, "**/Orthologues")
        pairs = orthofinder_pairs(orthologues.parent, gene_species, species)
        return pairs, sorted(orthologues.glob("Orthologues_*/*__v__*.tsv"))
    if method == "orthofinder_sequence_only":
        clusters = unique_path(output, "**/clusters_OrthoFinder_I*.txt_id_pairs.txt")
        mapping = unique_path(output, "**/SequenceIDs.txt")
        groups = read_checkpoint(clusters, mapping, gene_species)
        return list(group_pairs(groups.values(), gene_species)), [clusters, mapping]
    raise ValueError("Unknown frozen method")
