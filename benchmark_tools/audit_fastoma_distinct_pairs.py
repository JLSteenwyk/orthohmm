"""Audit disk-backed conversion against both retained historical QfO pair sets."""

import argparse
import filecmp
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.fastoma_distinct_pairs import write_pairs
from benchmark_tools.fastoma_to_pairwise import input_owners
from benchmark_tools.inventory_swiss_sequences import PREPARED_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.qfo_filter_pairs import filter_pairs, load_mapping
from benchmark_tools.run_simulation_methods import read_frozen


def compare_sets(actual, historical, destination):
    """Use independent external sorting, not the converter's SQLite index."""
    with destination.open("x") as stream:
        subprocess.run(["/usr/bin/sort", "--parallel=1", "-S", "64M", "-T",
                        str(destination.parent), "-u", str(historical)],
                       stdout=stream, env={**os.environ, "LC_ALL": "C"}, check=True)
    if not filecmp.cmp(actual, destination, shallow=False):
        raise ValueError(f"Distinct FastOMA pair set differs from {historical}")


def audit(root, work, output):
    if work.exists() or output.exists():
        raise FileExistsError("Require fresh audit work and report paths")
    prepared_path = root / "benchmark_tools/results/qfo_factorial_prepared_20260917.json"
    prepared = read_frozen(prepared_path, PREPARED_SHA)
    historical = root / "qfo_benchmark/results/fastoma"
    native = historical / "output/orthologs.tsv.gz"
    mapping = root / "qfo_benchmark/benchmark-webservice/reference_data/2020/mapping.json.gz"
    paths = [prepared_path, native, historical / "pairs.tsv", historical / "pairs.qfo.tsv", mapping,
             Path(__file__), Path(__file__).with_name("fastoma_distinct_pairs.py"),
             Path(__file__).with_name("fastoma_to_pairwise.py"),
             Path(__file__).with_name("qfo_filter_pairs.py"), Path("/usr/bin/sort")]
    records = [record(path) for path in paths] + prepared["input_fastas"]
    for item in records:
        check(item)
    owners = input_owners([Path(row["path"]) for row in prepared["input_fastas"]])
    if len(owners) != 976504 or len(set(owners.values())) != 78:
        raise ValueError("Not the original QfO input universe")
    work.mkdir(parents=True)
    raw, filtered = work / "distinct.tsv", work / "distinct.qfo.tsv"
    with raw.open("x") as stream:
        stats = write_pairs(native, owners, stream, work)
    compare_sets(raw, historical / "pairs.tsv", work / "historical.sorted.tsv")
    total, retained = filter_pairs(raw, filtered, load_mapping(mapping))
    if total != stats["distinct_pairs"]:
        raise ValueError("Filtered input count differs from distinct native count")
    compare_sets(filtered, historical / "pairs.qfo.tsv", work / "historical.qfo.sorted.tsv")
    for item in records:
        check(item)
    report = {"status": "distinct_fastoma_pairs_match_historical_raw_and_filtered_sets",
              "inputs": records, "input_accessions": len(owners), "input_species": 78,
              "counts": {**stats, "mapping_retained": retained, "mapping_removed": total - retained},
              "outputs": [record(path) for path in sorted(work.glob("*.tsv"))],
              "corrected_release_result": False, "accuracy_evaluated": False,
              "publication_ready": False,
              "limitations": ["Historical original-release conversion validation only.",
                              "No native workflow-completeness or biological correctness claim.",
                              "No historical pairs or scores were changed.",
                              "Runtime and memory were not measured under matched conditions."]}
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit(args.root.resolve(), args.work.resolve(), args.output.resolve())
