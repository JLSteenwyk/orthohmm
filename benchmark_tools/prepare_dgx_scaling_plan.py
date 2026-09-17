"""Translate the frozen 27-command panel to DGX paths without authorizing execution."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_scaling_commands import configurations

INPUT_SHA = "fa95fa60494c12d55f3bc277629c842ddbd2646a691c973882059708714e815e"
ORIGINAL_SHA = "25345e7a5d49e7474b09188498dc4760b298512266cd21510fba56ad6401e53d"
BASELINE_SHA = "bf677728cf9312edd9755d9ac1ceefbce3abb8c3848d43413631dc00ca00eb2f"
ROOT = Path("/home/jlsteenwyk/projects/orthohmm-publication")


def remap(value, mappings):
    if not value.startswith("/"):
        return value
    path = Path(value)
    for source, target in sorted(mappings.items(), key=lambda item: len(item[0].parts), reverse=True):
        if path.is_relative_to(source):
            return str(target / path.relative_to(source))
    raise ValueError("Unmapped absolute command path: " + value)


def translated_native(original, mappings, method):
    argv = [remap(value, mappings) for value in original]
    for flag in (("-t", "-a") if method == "orthofinder_full" else ("-c",)):
        if argv.count(flag) != 1 or argv[argv.index(flag) + 1] != "32":
            raise ValueError("Unexpected original CPU flags")
        argv[argv.index(flag) + 1] = "20"
    return argv


def plan(inputs, original, baseline):
    if inputs["status"] != "nested_inputs_materialized_unrun" or inputs["inference_started"]:
        raise ValueError("Unexpected transferred input status")
    core, output = ROOT / "core_arm_v2", ROOT / "scaling_native_v1"
    python = ROOT / "envs/orthohmm/bin/python"
    orthofinder = ROOT / "envs/orthofinder/bin/orthofinder"
    target = {"core_root": str(core), "tool_entrypoints": {
        "orthohmm_python": {"absolute_path": str(python)},
        "orthofinder": {"absolute_path": str(orthofinder)}}}
    runs = configurations(inputs, target, output, cpu_count=20)
    old_runs = original["runs"]
    if len(old_runs) != 27:
        raise ValueError("Changed frozen run count")
    mappings = {Path(baseline["core_root"]): core,
                Path(baseline["tool_entrypoints"]["orthohmm_python"]["absolute_path"]): python,
                Path(baseline["tool_entrypoints"]["orthofinder"]["absolute_path"]): orthofinder,
                Path(old_runs[0]["configuration"]["output"]).parents[1]: output}
    for dataset in inputs["datasets"]:
        size = dataset["proteomes"]
        matching = [r["dataset"] for r in old_runs if r["proteomes"] == size]
        if not matching or any(d != matching[0] for d in matching):
            raise ValueError("Ambiguous original dataset")
        old = matching[0]
        def identities(data):
            rows = [(Path(r["path"]).name, r["bytes"], r["sha256"]) for r in data["inputs"]]
            if len({r[0] for r in rows}) != size:
                raise ValueError("Duplicate or missing proteome basename")
            return sorted(rows)
        if identities(old) != identities(dataset) or any(old[k] != dataset[k] for k in ("proteins", "sequence_characters")):
            raise ValueError("Dataset membership, bytes or counts changed")
        mappings[Path(old["input_directory"])] = Path(dataset["input_directory"])
    common = {k: remap(v, mappings) for k, v in original["environment_overrides"].items()}
    environments = {
        "orthohmm": [ROOT / "envs/orthohmm/bin", ROOT / "external-prefix-v1/mafft-7.525/bin",
                      ROOT / "external-prefix-v1/fasttree-2.2.0"],
        "orthofinder": [ROOT / "envs/orthofinder/bin", ROOT / "diamond-2.0.13-build-v2",
                        ROOT / "mcl-prefix-v3/bin", ROOT / "external-prefix-v1/fasttree-2.1.11",
                        ROOT / "external-prefix-v1/mafft-7.525/bin", ROOT / "famsa-source-v1",
                        ROOT / "fastme-prefix-v1/bin"],
    }
    for before, after in zip(old_runs, runs):
        keys = ("index", "repeat", "proteomes", "method", "native_method")
        if any(before[k] != after[k] for k in keys):
            raise ValueError("Run identity/order changed")
        expected = translated_native(before["native_argv"], mappings, after["native_method"])
        if after["native_argv"] != expected:
            raise ValueError("Scientific command differs beyond declared paths/CPU flags")
        role = "orthofinder" if after["native_method"] == "orthofinder_full" else "orthohmm"
        after["environment_role"] = role
        after["preparation"] = ("fresh original-basename input copies; native outputs must not exist" if role == "orthofinder"
                                else "create fresh empty native output directory before inference")
    return {"status": "dgx_commands_prepared_not_authorized", "inference_started": False,
            "execution_authorized": False, "root": str(ROOT), "core_root": str(core),
            "core_commit": baseline["core_commit"], "output_root": str(output), "runs": runs,
            "environment_overrides": common,
            "environment_paths": {role: [*map(str, paths), "/usr/bin", "/bin"] for role, paths in environments.items()},
            "unset_environment": ["CONDA_PREFIX", "GOMP_CPU_AFFINITY", "OMP_PROC_BIND", "OMP_PLACES", "OMP_DYNAMIC"],
            "resource_plan": {"host": "spark-7ff0", "partition": "spark", "cpus": 20, "memory_gib": 96,
                              "exclusive": True, "concurrency": 1, "wall_limit_hours_per_run": 24,
                              "gpu": False, "fresh_output_each_run": True, "shared_search_cache": False},
            "measurement_plan": {"resource_interval_s": 1, "host_interval_s": 30,
                                 "boundary": "Native startup through exit; preparation, hashing, conversion and scoring separate"},
            "remaining_gates": ["Completed overhead-panel admission, not calibration alone",
                                "Remote interpreter/package/native/core/companion-file identity freeze and before/after verification",
                                "Actual frozen-tool FASTA enumeration snapshot on DGX, checked unchanged before every run",
                                "Remote input hashes and absent output directories before every run",
                                "Tested execution wrapper with separate preparation/native/conversion/scoring boundaries",
                                "GNU time and scheduler companion accounting with validated measurement boundary",
                                "Independent resource/host/native-output admission and preserved failure inventory"],
            "limitations": ["Paths and settings are prepared locally; this does not inspect remote state or authorize execution.",
                            "Manifest listing order differs from native file enumeration: frozen OrthoHMM fetch_fasta_files uses unsorted glob.",
                            "ARM timings must not be pooled with historical x86 timings.",
                            "Nested taxon composition co-varies with dataset size; no general scaling law is implied."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    paths = [args.results / name for name in ("dgx_scaling_inputs_20260917.json",
             "publication_scaling_commands_20260917.json", "publication_variable_native_methods_20260916.json")]
    records = [record(path) for path in paths]
    if [r["sha256"] for r in records] != [INPUT_SHA, ORIGINAL_SHA, BASELINE_SHA]:
        raise ValueError("Changed frozen plan inputs")
    result = plan(*(json.loads(path.read_text()) for path in paths))
    for item in records:
        check(item)
    result.update(inputs=records[0], original_commands=records[1], baseline=records[2], source=record(Path(__file__)),
                  helper_sources=[record(Path(__file__).with_name(name)) for name in
                                  ("prepare_scaling_commands.py", "prepare_simulation_methods.py", "benchmark_production.py")])
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
