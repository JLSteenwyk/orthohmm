"""Freeze retained simulation inputs and ancestral-family search-recall denominators."""

import argparse
from collections import Counter
import json
from pathlib import Path

from Bio import SeqIO

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


PINS = {
    "publication_variable_simulation_manifest_20260916.json": "806aa1e5f6976c323ff2f6641dd88e25264e7417749e7d77565294666aee776b",
    "simulation_variable_native_results_20260916.json": "cc99fc31c3433809098212d3c8dc12f4ad829b4cfeb66dabcb13aede4e95258f",
    "ob_sequence_search_prepared_20260916.json": "9edb1bb8232e114e12a6199412dd45a31d740c8a1dcba36fcf05d2c58b23095b",
    "publication_frozen_overlay_install_20260926.json": "467036ee88d5e6f15d322aea4c94397591ef7ba521ac9f8ddde4b7cb783b97d2",
}
CONDITIONS = {"baseline", "divergent", "turnover", "divergent_turnover", "missing20",
              "uneven_taxa", "taxon_count_control"}
EVALUE_GRID = [1e-100, 1e-80, 1e-60, 1e-40, 1e-30, 1e-20, 1e-15, 1e-10,
               1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]


def family_denominators(families, gene_species):
    seen, result = set(), {}
    for family, genes in sorted(families.items()):
        if len(genes) != len(set(genes)) or set(genes) & seen or not set(genes) <= gene_species.keys():
            raise ValueError("Repeated or foreign ancestral-family gene")
        seen.update(genes)
        counts = Counter(gene_species[gene] for gene in genes)
        result[family] = len(genes) ** 2 - sum(count ** 2 for count in counts.values())
    if seen != gene_species.keys():
        raise ValueError("Ancestral families do not cover every input gene")
    return result


def dataset(row, expected_truth_sha):
    truth_path = Path(row["truth"])
    truth_record = record(truth_path)
    if truth_record["sha256"] != expected_truth_sha:
        raise ValueError("Simulation truth differs from admitted historical identity")
    truth = json.loads(truth_path.read_text())
    inputs, gene_species, lengths = [], {}, {}
    listed = truth.get("prepared_inputs", truth["inputs"])
    for item in listed:
        path = (truth_path.parent / item["path"]).resolve()
        if path.parent != Path(row["input"]).resolve():
            raise ValueError("Input lies outside declared dataset")
        observed = record(path)
        if (observed["bytes"], observed["sha256"]) != (item["bytes"], item["sha256"]):
            raise ValueError("Simulation FASTA differs from retained truth manifest")
        inputs.append(observed)
        for sequence in SeqIO.parse(path, "fasta"):
            if not sequence.id or sequence.id in gene_species or not len(sequence.seq):
                raise ValueError("Empty or repeated simulation sequence")
            gene_species[sequence.id] = path.stem
            lengths[sequence.id] = len(sequence.seq)
    paths = [item["path"] for item in inputs]
    if (len(paths) != len(set(paths))
            or {str(p.resolve()) for p in Path(row["input"]).iterdir()} != set(paths)
            or len(gene_species) != truth["extant_genes"]
            or set(gene_species.values()) != set(truth["species"])):
        raise ValueError("Simulation universe/inventory mismatch")
    denominators = family_denominators(truth["families"], gene_species)
    eligible = sum(value > 0 for value in denominators.values())
    if not eligible:
        raise ValueError("No eligible cross-species homology family")
    species_counts = Counter(gene_species.values())
    possible = len(gene_species) ** 2 - sum(count ** 2 for count in species_counts.values())
    for item in [truth_record, *inputs]:
        check(item)
    return dict(condition=row["condition"], seed=row["seed"], input=row["input"],
        split="calibration" if row["seed"] <= 20261105 else "reporting", truth=truth_record,
        inputs=inputs, genes=len(gene_species), species=len(species_counts),
        family_directed_homology_pairs=denominators, eligible_families=eligible,
        directed_homology_pairs=sum(denominators.values()),
        directed_nonhomology_pairs=possible - sum(denominators.values()),
        min_length=min(lengths.values()), max_length=max(lengths.values()))


def prepare(repo, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    reports, checked = {}, []
    for name, sha in PINS.items():
        path = repo / "benchmark_tools/results" / name
        item = record(path)
        if item["sha256"] != sha:
            raise ValueError("Changed source evidence: " + name)
        checked.append(item)
        reports[name] = json.loads(path.read_text())
    manifest = reports["publication_variable_simulation_manifest_20260916.json"]
    keys = {(r["condition"], r["seed"]) for r in manifest["datasets"]}
    if (len(manifest["datasets"]) != 70
            or keys != {(condition, seed) for condition in CONDITIONS for seed in range(20261101, 20261111)}):
        raise ValueError("Expected complete seven-condition ten-seed panel")
    historical = reports["simulation_variable_native_results_20260916.json"]
    truth_pins = {}
    for row in historical["records"]:
        key = row["condition"], row["seed"]
        if "truth_sha256" in row:
            old = truth_pins.setdefault(key, row["truth_sha256"])
            if old != row["truth_sha256"]:
                raise ValueError("Conflicting historical truth hashes")
    if set(truth_pins) != keys:
        raise ValueError("Incomplete historical truth identity")
    datasets = [dataset(row, truth_pins[row["condition"], row["seed"]])
                for row in sorted(manifest["datasets"], key=lambda row: (row["condition"], row["seed"]))]
    diamond = reports["ob_sequence_search_prepared_20260916.json"]["diamond"]
    check(diamond)
    for item in reports["publication_frozen_overlay_install_20260926.json"]["checked_records"]:
        check(item)
    for item in checked:
        check(item)
    result = dict(status="inputs_verified_search_not_executed", source=record(__file__), checked_evidence=checked,
        datasets=datasets, diamond=diamond, evalue_grid=EVALUE_GRID,
        hmm_settings=dict(matrix="BLOSUM62", kmer_k=4, max_candidates_per_query=100,
                          band_width=64, evalue_threshold=1e-4, cpus=4, threads_per_worker=1),
        diamond_settings=dict(sensitivity="very-sensitive", matrix="BLOSUM62", gapopen=11,
            gapextend=1, comp_based_stats=1, masking=1, max_target_seqs=0, max_hsps=1, evalue=1.0, threads=4),
        primary_statistic="Equal-weight mean of family directed cross-species ancestral-homology recall within dataset, then datasets.",
        cutoff_selection="Calibration-only minimum absolute difference from fixed HMM primary recall; exact tie uses smaller E-value.",
        reporting_match_gate=dict(overall_absolute_recall_difference_max=.02, each_condition_absolute_difference_max=.05),
        new_search_results_inspected=False, development_exposed=True, independent_validation=False,
        scientific_defaults_changed=False, publication_ready=False,
        limitations=["Calibration/reporting seeds separate this diagnostic, not prior development exposure.",
                     "Ancestral homology includes paralogs; it is not orthology F1 or calibrated real-data sensitivity.",
                     "Conditions sharing a seed share histories; do not treat them as independent replicates.",
                     "A single threshold can fail the match gate; do not extend the grid after results.",
                     "Post-search filtering does not match computational effort; shared-host timing remains descriptive."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.repo.resolve(), args.output.absolute())
