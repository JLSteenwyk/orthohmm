"""Admit sequence/alignment evidence and freeze outcome-independent RefOG strata."""

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.assemble_orthobench_factorial import load_reference_snapshot
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.prepare_ob_reference_alignments import normalize, validate_alignment, mean_canonical_identity, tool_inventory
from benchmark_tools.prepare_ob_sequence_features import sequence_features, FIELDS
from benchmark_tools.recover_factorial_postflight import PREPARED_HASH, ENVIRONMENT_HASH
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_ygob_validation import require_completed_job

FEATURE_SHA = "4e3a05668a4f0e69c56472e7bf3836a120f4d76b4b07da576418d869838fd93c"
CATEGORIES = {"size": ("small_2_20", "medium_21_50", "large_gt_50"),
              "copy_number": ("single_copy", "multi_copy"),
              "identity": ("lower_identity", "higher_identity", "missing"),
              "relative_length": ("short_relative", "not_short_relative", "missing"),
              "composition": ("concentrated", "not_concentrated", "missing")}


def assign_strata(families):
    values = [row["mean_pairwise_identity"] for row in families.values() if row["mean_pairwise_identity"] is not None]
    if any(isinstance(value, bool) or not math.isfinite(value) or not 0 <= value <= 1 for value in values):
        raise ValueError("Invalid family identity")
    cutoff = statistics.median(values) if values else None
    strata = {dimension + ":" + category: [] for dimension, categories in CATEGORIES.items() for category in categories}
    assignments = {}
    for name, row in sorted(families.items()):
        lengths = row["residue_lengths"]
        flags = row["composition_flags"]
        size = row["genes"]
        if (type(size) is not int or size < 2 or len(lengths) != size or len(flags) != size or not row["species_counts"]
                or sum(row["species_counts"].values()) != size or any(type(n) is not int or n < 1 for n in row["species_counts"].values())
                or any(isinstance(n, bool) or not isinstance(n, int) or n < 0 for n in lengths)
                or any(flag is not None and type(flag) is not bool for flag in flags)):
            raise ValueError("Invalid family feature inventory")
        length_category = "missing" if 0 in lengths else "short_relative" if min(lengths) < .5 * statistics.median(lengths) else "not_short_relative"
        composition = "concentrated" if any(flag is True for flag in flags) else "missing" if any(flag is None for flag in flags) else "not_concentrated"
        identity = row["mean_pairwise_identity"]
        categories = {"size": "small_2_20" if size <= 20 else "medium_21_50" if size <= 50 else "large_gt_50",
                      "copy_number": "multi_copy" if max(row["species_counts"].values()) > 1 else "single_copy",
                      "identity": "missing" if identity is None else "lower_identity" if identity <= cutoff else "higher_identity",
                      "relative_length": length_category, "composition": composition}
        assignments[name] = categories
        for dimension, category in categories.items():
            strata[dimension + ":" + category].append(name)
    selected = {key: min(names, key=lambda name: (hashlib.sha256(name.encode("ascii")).hexdigest(), name)) if names else None
                for key, names in strata.items()}
    return {"identity_median": cutoff, "strata": strata, "assignments": assignments,
            "illustrative_families_by_stratum": selected,
            "illustrative_families_deduplicated": sorted({name for name in selected.values() if name is not None}),
            "planned_comparisons": 2, "planned_metrics": 3, "multiplicity_endpoints": 84}


def check_alignment_panel(panel, references, source):
    preflight = panel["preflight"]
    names = [run["refog"] for run in panel["runs"]]
    if (panel["status"] != "reference_alignments_prepared_unscored" or panel["accuracy_evaluated"] is not False
            or panel["failed_families"] != [] or preflight["job_id"] != "21309"
            or preflight["accuracy_evaluated"] is not False or preflight["status"] != "ready_unscored"
            or preflight["families"] != sorted(references) or sorted(names) != sorted(references)
            or len(set(names)) != len(names) or preflight["source"] != source):
        raise ValueError("Incomplete, failed or mismatched reference alignment panel")
    if any(run["status"] != "alignment_validated" or run["exit_code"] != 0
           or run["accuracy_evaluated"] is not False for run in panel["runs"]):
        raise ValueError("A reference family was not successfully validated")


def verify_executor(root, name, revision):
    directory = root / "benchmarks/work" / name
    expected = subprocess.check_output(["git", "-C", str(root), "rev-parse", revision + "^{commit}"], text=True).strip()
    actual = subprocess.check_output(["git", "-C", str(directory), "rev-parse", "HEAD"], text=True).strip()
    if expected != actual:
        raise ValueError("Prepared feature executor changed")
    subprocess.run(["git", "-C", str(directory), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    return directory


def prepare(root, output):
    root, output = root.resolve(), output.resolve()
    if output.exists():
        raise FileExistsError(output)
    accounting = subprocess.check_output(["sacct", "-j", "21308,21309", "--parsable2", "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = {str(job): require_completed_job(accounting, job) for job in (21308, 21309)}
    results = root / "benchmark_tools/results"
    feature_path = results / "ob_sequence_features_prepared_20260916.json"
    features = read_frozen(feature_path, FEATURE_SHA)
    feature_executor = verify_executor(root, "publication_ob_sequence_features_v1", "8bccb33")
    alignment_executor = verify_executor(root, "publication_ob_reference_alignments_v1", "a753d0d")
    if (features["status"] != "sequence_features_prepared_unscored" or features["job_id"] != "21308"
            or features["accuracy_evaluated"] is not False or features["genes"] != 251378
            or features["source"] != file_provenance(feature_executor / "benchmark_tools/prepare_ob_sequence_features.py")):
        raise ValueError("Unexpected sequence feature provenance")
    prepared_path = results / "orthobench_factorial_prepared_20260916.json"
    env_path = results / "publication_variable_native_methods_20260916.json"
    prepared = read_frozen(prepared_path, PREPARED_HASH)
    environment = read_frozen(env_path, ENVIRONMENT_HASH)
    if features["inputs"] != prepared["fasta_inputs"] or features["input_manifest"] != file_provenance(prepared_path):
        raise ValueError("Sequence features use different FASTAs")
    references, _, _, records = load_reference_snapshot(results / "orthobench_paired_uncertainty_20260916.json")
    alignment_root = root / "benchmarks/results/ob_reference_alignments_v1"
    alignment_path = alignment_root / "manifest.json"
    alignment_record = file_provenance(alignment_path)
    panel = json.loads(alignment_path.read_text())
    preflight_record = file_provenance(alignment_root / "preflight.json")
    if json.loads((alignment_root / "preflight.json").read_text()) != panel["preflight"]:
        raise ValueError("Alignment preflight differs from completed manifest")
    check_alignment_panel(panel, references, file_provenance(alignment_executor / "benchmark_tools/prepare_ob_reference_alignments.py"))
    mafft = Path(environment["tool_entrypoints"]["mafft"]["absolute_path"])
    prefix, tools = tool_inventory(mafft)
    preflight = panel["preflight"]
    if (preflight["tools"] != tools or preflight["version"] != "v7.525 (2024/Mar/13)"
            or preflight["reference_records"] != records or preflight["fasta_inputs"] != prepared["fasta_inputs"]
            or preflight["input_manifest"] != file_provenance(prepared_path) or preflight["environment_manifest"] != file_provenance(env_path)
            or preflight["environment_overrides"] != {"MAFFT_BINARIES": str(prefix), "OMP_NUM_THREADS": "1",
                 "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "LC_ALL": "C"}):
        raise ValueError("Alignment environment or input provenance changed")
    verify_file(mafft, environment["tool_entrypoints"]["mafft"])
    wanted, sequences, gene_species = set().union(*references.values()), {}, {}
    for item in prepared["fasta_inputs"]:
        verify_file(Path(item["path"]), item)
        for protein in SeqIO.parse(item["path"], "fasta"):
            if protein.id in gene_species:
                raise ValueError("Duplicate protein identifier")
            gene_species[protein.id] = Path(item["path"]).name
            if protein.id in wanted:
                sequences[protein.id] = str(protein.seq).upper()
    if set(sequences) != wanted or len(gene_species) != 251378 or preflight["reference_genes"] != len(wanted) or preflight["workers"] != 8:
        raise ValueError("Incomplete protein/reference universe")
    verify_file(Path(features["table"]["path"]), features["table"])
    seen, ref_features = set(), {}
    with Path(features["table"]["path"]).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != list(FIELDS):
            raise ValueError("Sequence feature columns changed")
        for row in reader:
            gene = row["gene"]
            if gene in seen or gene not in gene_species or row["proteome"] != gene_species[gene]:
                raise ValueError("Sequence feature gene/species inventory differs")
            seen.add(gene)
            if gene in wanted:
                value = sequence_features(sequences[gene])
                if any(row[key] != ("NA" if item is None else str(item)) for key, item in value.items()):
                    raise ValueError("Reference protein feature does not match its sequence")
                ref_features[gene] = value
    if seen != set(gene_species):
        raise ValueError("Sequence feature table is incomplete")
    families, tracked = {}, [alignment_record, preflight_record, features["table"], *records, *prepared["fasta_inputs"], *tools]
    for run in panel["runs"]:
        name, genes = run["refog"], references[run["refog"]]
        directory = alignment_root / Path(name).stem
        command = [str(mafft), "--amino", "--auto", "--thread", "1", str(directory / "input.faa")]
        if run["command"] != command or run["input"]["path"] != command[-1] or run["alignment"]["path"] != str(directory / "aligned.faa"):
            raise ValueError("Reference alignment command or paths differ")
        if json.loads((directory / "status.json").read_text()) != run:
            raise ValueError("Per-family status differs from panel manifest")
        tracked.append(file_provenance(directory / "status.json"))
        for key in ("input", "alignment"):
            verify_file(Path(run[key]["path"]), run[key])
            tracked.append(run[key])
        expected = {gene: normalize(sequences[gene]) for gene in genes}
        input_rows = list(SeqIO.parse(command[-1], "fasta"))
        if len(input_rows) != len(genes) or {row.id: str(row.seq) for row in input_rows} != expected:
            raise ValueError("Alignment input differs from normalized reference sequences")
        aligned = validate_alignment(Path(run["alignment"]["path"]), expected)
        identity = mean_canonical_identity(aligned)
        if identity != run["identity"] or run["genes"] != len(genes) or run["columns"] != len(next(iter(aligned.values()))):
            raise ValueError("Recomputed alignment descriptor differs")
        removed = {gene: sequences[gene].count("*") for gene in sorted(genes) if "*" in sequences[gene]}
        if removed != run["removed_stop_symbols"]:
            raise ValueError("Stop-symbol handling differs")
        families[name] = {"genes": len(genes), "gene_ids": sorted(genes),
            "species_counts": dict(Counter(gene_species[gene] for gene in genes)),
            "residue_lengths": [ref_features[gene]["residue_length"] for gene in sorted(genes)],
            "composition_flags": [ref_features[gene]["composition_concentrated"] for gene in sorted(genes)],
            "mean_pairwise_identity": identity["mean_pairwise_identity"], "alignment": run["alignment"]}
    report = assign_strata(families)
    for item in tracked:
        verify_file(Path(item["path"]), item)
    report.update(status="family_strata_prepared_unscored", accuracy_evaluated=False, families=families,
                  scheduler=scheduler, sequence_features=file_provenance(feature_path), reference_alignments=alignment_record,
                  source=file_provenance(Path(__file__)),
                  helper_sources=[file_provenance(Path(__file__).with_name(name)) for name in
                      ("prepare_ob_reference_alignments.py", "prepare_ob_sequence_features.py", "assemble_orthobench_factorial.py")],
                  protocol=file_provenance(results / "ORTHOBENCH_ERROR_ANALYSIS_PROTOCOL_20260916.md"),
                  limitations=["Outcome-independent descriptors on development-exposed reference families, not independent confirmation.",
                      "Copy number, relative length and global composition are not verified duplication history, fragments or domain architecture.",
                      "Overlapping feature dimensions are not independent samples; empty and small bins remain explicit."])
    output.mkdir(parents=True)
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
