"""Freeze default-mode SonicParanoid for the corrected QfO input release."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_qfo_corrected_primary import verify as verify_primary
from benchmark_tools.prepare_qfo_corrected_proteinortho import PRIMARY_SHA, REGISTRY_SHA
from benchmark_tools.snapshot_runtime_trees import verify as verify_tree

DEFAULT_MODE = ["default", "diamond", "very-sensitive", 0.0]
TOOLS = {"blastp", "makeblastdb", "diamond", "mcl", "mmseqs"}


def command(entrypoint, root):
    root = Path(root)
    if not root.is_absolute():
        raise ValueError("Require absolute output root")
    return [str(entrypoint), "-i", str(root / "input"), "-o", str(root / "output"), "-t", "32"]


def validate_runtime(runtime):
    if (runtime["status"] != "read_only_current_resolution_observed"
            or runtime["version"] != "2.0.9" or runtime["default_mode"] != DEFAULT_MODE
            or set(runtime["tools"]) != TOOLS):
        raise ValueError("Unexpected SonicParanoid version, default mode or dependencies")
    if any(tool["exit_code"] != 0 for tool in runtime["tools"].values()):
        raise ValueError("SonicParanoid dependency probe failed")


def prepare(root, output, destination, runtime_path, runtime_sha):
    if output.exists() or destination.exists():
        raise FileExistsError("Require new native output root and manifest")
    results = root / "benchmark_tools/results"
    primary_path = results / "qfo_corrected_primary_commands_20260918.json"
    primary = read_frozen(primary_path, PRIMARY_SHA)
    _, inputs = verify_primary(primary)
    registry = results / "publication_comparison_orthomcl_complete_20260916.json"
    read_frozen(registry, REGISTRY_SHA)
    runtime = read_frozen(runtime_path, runtime_sha)
    validate_runtime(runtime)
    verify_tree(runtime["package_inventory"])
    entrypoint = Path("/home/bizon/anaconda3/bin/sonicparanoid")
    sources = [record(root / path) for path in (
        "qfo_benchmark/run_tool.slurm", "benchmark_tools/inspect_sonicparanoid_runtime.py",
        "benchmark_tools/sonicparanoid_to_pairwise.py", "benchmark_tools/normalize_three_kingdoms_orthogroups.py",
        "qfo_benchmark/filter_qfo_pairs.py", "benchmark_tools/run_qfo_corrected_primary.py")]
    records = [record(primary_path), record(registry), record(runtime_path), record(entrypoint),
               runtime["python"], *[tool["file"] for tool in runtime["tools"].values()], *sources, *inputs]
    for item in records:
        check(item)
    report = {
        "status": "corrected_sonic_command_frozen_unrun", "execution_authorized": False,
        "accuracy_admitted": False, "source": record(__file__), "checked_records": records,
        "primary_manifest": record(primary_path), "runtime_manifest": record(runtime_path),
        "input_directory": primary["input_directory"], "input_fastas": inputs,
        "output_root": str(output), "copy_inputs_to": str(output / "input"),
        "cwd": str(output), "native_argv": command(entrypoint, output),
        "default_mode": DEFAULT_MODE, "search_reuse": False,
        "resources": {"node": "bizon", "cpus": 32, "memory_gib": 192, "time_limit_hours": 72},
        "pair_semantics": "Native species-to-species ortholog tables; do not clique-expand global groups.",
        "native_settings_required": {"minimum_bitscore": 40, "max_inparalog_length_difference": 0.75,
            "ortholog_merging_threshold": 0.75, "mcl_inflation": 1.5,
            "only_graph_based_orthology": False, "proteomes": 78, "threads": 32},
        "remaining_gates": [
            "Pin full interpreter/transitive/system runtime inventory and explicit execution environment.",
            "Pinned fresh-copy launcher with effective resolver and input checks before/after execution.",
            "Isolate bytecode lookup and disable writes/user-site; reject unreviewed loader injection.",
            "Native completion, all species-pair tables and ID ownership validation before conversion/scoring.",
            "Freeze conversion, corrected participant identity and scorer invocation before evaluation.",
        ],
        "limitations": ["Current dependency resolution is not historical process execution proof.",
                        "Default DIAMOND mode also permits downstream MMseqs profile searches.",
                        "Shared-host execution is not dedicated timing evidence."],
    }
    with destination.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "output-root", "manifest", "runtime"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--runtime-sha256", required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.output_root.resolve(), args.manifest.resolve(),
            args.runtime.resolve(), args.runtime_sha256)
