"""Prepare a deterministic, prediction-independent nested proteome scaling panel."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen

INPUT_SHA = "5c325f4d77865e0c7571fe4bb4df0be0959977a49e1d22169f192518700f9382"
ORDER_SALT = "orthohmm-publication-scaling-20260916-v1"
SIZES = (4, 8, 12)
METHODS = ("orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_3_1_5_full")


def ordered_inputs(inputs):
    names = [Path(item["path"]).name for item in inputs]
    if len(inputs) != 12 or len(set(names)) != 12 or len({item["path"] for item in inputs}) != 12:
        raise ValueError("Require twelve unique proteome files and basenames")
    return sorted(inputs, key=lambda item: (hashlib.sha256(
        (ORDER_SALT + "\n" + Path(item["path"]).name).encode("utf-8")).hexdigest(), Path(item["path"]).name))


def planned_runs():
    # Rotate method order within each size/repeat block to balance execution order.
    runs = []
    for repeat in range(3):
        for size_index, size in enumerate(SIZES):
            offset = (repeat + size_index) % len(METHODS)
            for method in METHODS[offset:] + METHODS[:offset]:
                runs.append({"index": len(runs), "repeat": repeat, "proteomes": size, "method": method})
    return runs


def prepare(root, output):
    if output.exists():
        raise FileExistsError(output)
    source = root / "benchmark_tools/results/orthobench_factorial_prepared_20260916.json"
    frozen = read_frozen(source, INPUT_SHA)
    inputs = ordered_inputs(frozen["fasta_inputs"])
    for item in inputs:
        check(item)
    seen, proteomes = set(), []
    for item in inputs:
        proteins = residues = 0
        for protein in SeqIO.parse(item["path"], "fasta"):
            if not protein.id or protein.id in seen or not len(protein.seq):
                raise ValueError("Empty sequence or duplicate/empty protein identifier")
            seen.add(protein.id)
            proteins += 1
            residues += len(protein.seq)
        if not proteins:
            raise ValueError("Empty proteome")
        proteomes.append({"input": item, "proteins": proteins, "sequence_characters": residues})
    if len(seen) != 251378:
        raise ValueError("Changed complete OrthoBench universe")
    output.mkdir(parents=True)
    datasets = []
    for size in SIZES:
        directory = output / f"proteomes_{size}" / "input"
        directory.mkdir(parents=True)
        selected = proteomes[:size]
        for row in selected:
            target = Path(row["input"]["path"])
            (directory / target.name).symlink_to(target)
        actual = sorted((record(path) for path in directory.iterdir()), key=lambda item: item["path"])
        expected = sorted((row["input"] for row in selected), key=lambda item: item["path"])
        if actual != expected:
            raise ValueError("Prepared dataset differs from intended full proteomes")
        datasets.append({"proteomes": size, "input_directory": str(directory), "inputs": expected,
                         "proteins": sum(row["proteins"] for row in selected),
                         "sequence_characters": sum(row["sequence_characters"] for row in selected)})
    for item in inputs:
        check(item)
    report = {"status": "scaling_inputs_prepared_unrun", "accuracy_evaluated": False,
              "source": record(__file__), "input_manifest": record(source), "order_salt": ORDER_SALT,
              "ordered_proteomes": proteomes, "datasets": datasets, "planned_runs": planned_runs(),
              "resource_plan": {"cpus": 32, "memory_gib": 128, "wall_limit_hours_per_run": 24,
                                "fresh_output_each_run": True, "shared_search_cache": False},
              "prediction_files_read": [], "accuracy_reference_files_read": [],
              "limitations": ["One nested OrthoBench taxon series; input size and taxon composition co-vary.",
                  "No subsampling of proteins within a selected proteome; no added or rewritten sequences.",
                  "Preparation is not timing evidence. Native runtime, exact commands, output semantics and workload gates are required before execution.",
                  "No universal scaling law or resource advantage is established by this panel."]}
    (output / "manifest.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.output.resolve())
