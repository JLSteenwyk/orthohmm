"""Prepare corrected QfO search-control inputs without launching searches."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_qfo_corrected_primary import STAGING_SHA, INVENTORY_SHA
from benchmark_tools.prepare_sequence_search_control import search_command, write_queries
from benchmark_tools.run_sequence_search_control import MANIFEST_SHA
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def validate_inputs(stage, inventory):
    if (inventory["status"] != "corrected_staged_inventory_verified_pending_execution_freeze"
            or inventory["total_sequences"] != 984137 or inventory["proteomes"] != 78
            or inventory["inputs"][0]["sha256"] != STAGING_SHA
            or [r["file"] for r in inventory["files"]] != stage["input_fastas"]):
        raise ValueError("Wrong corrected QfO input inventory")
    inputs = stage["input_fastas"]
    paths = [Path(r["path"]) for r in inputs]
    if len(paths) != 78 or len(set(paths)) != 78 or len({p.parent for p in paths}) != 1:
        raise ValueError("Require 78 unique files in one corrected input directory")
    return sorted(inputs, key=lambda r: r["path"])


def search_plan(inputs, binary, queries, output):
    result = []
    for index, item in enumerate(inputs):
        directory = output / f"target_{index:02d}"
        database, hits = directory / "target", directory / "hits.tsv"
        result.append({"index": index, "target_fasta": item, "output": str(hits),
            "makedb": [str(binary), "makedb", "--in", item["path"], "--db", str(database), "--threads", "32"],
            "search": search_command(binary, queries, database, hits)})
    return result


def prepare(root, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    results = root / "benchmark_tools/results"
    stage_path = results / "qfo_corrected_staging_manifest_20260918.json"
    inventory_path = results / "qfo_corrected_staged_inventory_20260918.json"
    old_plan_path = root / "benchmarks/work/ob_sequence_search_control_v1/manifest.json"
    stage = read_frozen(stage_path, STAGING_SHA)
    inventory = read_frozen(inventory_path, INVENTORY_SHA)
    old_plan = read_frozen(old_plan_path, MANIFEST_SHA)
    inputs = validate_inputs(stage, inventory)
    binary = Path(old_plan["diamond"]["path"])
    checked = [record(stage_path), record(inventory_path), record(old_plan_path), old_plan["diamond"], *inputs]
    for item in checked:
        check(item)
    version = subprocess.check_output([str(binary), "version"], text=True).strip()
    if version != "diamond version 2.1.11":
        raise ValueError("Unexpected frozen DIAMOND version")
    output.mkdir(parents=True, exist_ok=False)
    queries = output / "queries.fasta"
    owners, lengths = write_queries(inputs, queries)
    if len(owners) != 984137 or len(set(owners.values())) != 78:
        raise ValueError("Corrected sequence universe changed")
    metadata = output / "gene_metadata.json"
    with metadata.open("x") as stream:
        json.dump({g: {"species": owners[g], "length": lengths[g]} for g in sorted(owners)}, stream, sort_keys=True)
        stream.write("\n")
    searches = search_plan(inputs, binary, queries, output)
    for search in searches:
        Path(search["output"]).parent.mkdir()
    for item in checked:
        check(item)
    report = {"status": "corrected_qfo_search_control_prepared_unrun", "source": record(__file__),
        "inputs": inputs, "checked_records": checked, "queries": record(queries), "gene_metadata": record(metadata),
        "diamond": record(binary), "diamond_version": version, "searches": searches,
        "genes": len(owners), "proteomes": 78, "accuracy_evaluated": False, "execution_authorized": False,
        "normalization": old_plan["normalization"], "helpers": [record(Path(__file__).with_name(n)) for n in
            ("prepare_sequence_search_control.py", "prepare_qfo_corrected_primary.py", "run_sequence_search_control.py")],
        "remaining_gates": ["Frozen corrected-only runner with full-panel completion and failure retention.",
            "Numeric conversion, graph replay and independent scoring after all 78 targets complete.",
            "Measure memory/CPU and hit-set diagnostics; equal E-values do not establish matched sensitivity or cost.",
            "All-hit and post-search top100 variants retain frozen graph settings; no endpoint-driven retuning.",
            "Original search manifest supplies binary/settings provenance only; no historical hit reuse."]}
    with (output / "manifest.json").open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.output.absolute())
