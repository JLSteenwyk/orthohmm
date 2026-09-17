"""Prepare three small native launcher smokes, never the scientific timing panel."""

import argparse
from copy import deepcopy
import json
from pathlib import Path
import sys

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.snapshot_orthohmm_input_order import record

ROOT = Path("/home/jlsteenwyk/projects/orthohmm-publication")


def inputs(results):
    original = read_pinned(results / "dgx_orthohmm_pipeline_smoke_20260917.json",
                           "c35b1f4228d6ff1af7d7459d2500a765604bf267f82fed54cdfe92d539eadf94")
    rows, count, chars = [], 0, 0
    directory = ROOT / "smoke_missing20_20261101_input"
    for row in original["inputs"]:
        if record(row["path"]) != row:
            raise ValueError("Changed original smoke input")
        for seq in SeqIO.parse(row["path"], "fasta"):
            count += 1
            chars += len(seq.seq)
        rows.append({**row, "path": str(directory / Path(row["path"]).name)})
    if count != 645 or len(rows) != 8:
        raise ValueError("Wrong smoke fixture")
    return {"datasets": [{"input_directory": str(directory), "inputs": rows,
                           "proteomes": 8, "proteins": count, "sequence_characters": chars}]}


def specification(results, order):
    plan = read_pinned(results / "dgx_scaling_commands_20260917.json",
                       "7096348236f8372ef7f6ad12a3e829eae3dfee812b17c3b0bdd582480538ff1b")
    dataset = inputs(results)["datasets"][0]
    if (len(order["datasets"]) != 1 or order["datasets"][0]["input_directory"] != dataset["input_directory"]
            or sorted(order["datasets"][0]["inputs_in_native_order"], key=lambda r: r["path"]) != dataset["inputs"]):
        raise ValueError("Native enumeration does not match the frozen smoke fixture")
    runs = []
    for original in plan["runs"][:3]:
        mappings = {original["dataset"]["input_directory"]: dataset["input_directory"],
                    str(ROOT / "scaling_native_v1"): str(ROOT / "launcher_smoke_v1")}
        def remap(value):
            if isinstance(value, dict):
                return {k: remap(v) for k, v in value.items()}
            if isinstance(value, list):
                return [remap(v) for v in value]
            if isinstance(value, str):
                for old, new in mappings.items():
                    if value == old or value.startswith(old + "/"):
                        return new + value[len(old):]
            return value
        run = remap(deepcopy(original))
        run.update(dataset=dataset, proteomes=8)
        runs.append(run)
    return {"purpose": "launcher_smoke", "execution_authorized": True, "runs": runs,
            "scientific_timing_runs_authorized": False, "orders": order["datasets"], "enumerator": order["source"],
            "environment_overrides": plan["environment_overrides"], "environment_paths": plan["environment_paths"],
            "unset_environment": plan["unset_environment"], "resource_plan": plan["resource_plan"],
            "native_timeout_s": 900, "collector_directory": str(ROOT / "collector_load_recipe_v2"),
            "runtime_manifests": [
                {"path": str(ROOT / "runtime_inventory_v1/trees.json"),
                 "sha256": "2f38fc57683e51a7b6768b4709db16293590983640c41cfe12a56ed6036750ae"},
                {"path": str(ROOT / "runtime_inventory_v1/system_trees.json"),
                 "sha256": "4083d15c0ffc40013756c4588763aa8af24ca39165622665ee5074662e4b46ce"}],
            "overhead": {"path": str(ROOT / "native_launcher_recipe_v1/benchmark_tools/dgx_verified_overhead_20260917.json"),
                         "sha256": "0184d54be2f55720aff35f27b6db0123de28ca7dfcdabd50cbf1d8d79d60b437"},
            "limitations": ["Engineering fixture with identical scientific flags but different input/output paths.",
                            "Three small smokes are not the27 scientific scaling measurements."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--order", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = specification(args.results, json.loads(args.order.read_text())) if args.order else inputs(args.results)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
