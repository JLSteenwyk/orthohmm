"""Freeze native inference boundaries for the prespecified 27-run scaling panel."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.benchmark_production import parse_args
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_scaling_inputs import planned_runs, METHODS
from benchmark_tools.prepare_simulation_methods import commands
from benchmark_tools.run_simulation_methods import read_frozen, verify_environment, execution_environment
from benchmark_tools.run_simulation_tree_mode_control import METHOD_SHA

INPUT_SHA = "1957b050d33dd89e933ff33f96500a4fbaa6b2155d85d067d7d7d956c3d139af"


def native_orthohmm(config):
    args = parse_args(config["argv"][2:])
    argv = [config["argv"][0], "-m", "orthohmm", str(args.input_directory.resolve()),
            "-o", str(args.output_directory.resolve()), "-c", str(args.cpu)]
    fields = (("--threads_per_worker", "threads_per_worker"), ("-x", "matrix"), ("-e", "evalue"),
              ("--clustering", "clustering"), ("--cpm_resolution", "cpm_resolution"),
              ("--refinement_profile", "refinement_profile"), ("--accuracy_profile", "accuracy_profile"))
    for flag, field in fields:
        argv += [flag, str(getattr(args, field))]
    argv += ["--metrics_json", str(args.result_json.resolve())]
    if not args.full_output:
        argv += ["--stop", "infer"]
    if args.phylogeny == "reconcile":
        for field in ("phylogeny", "species_tree_mode", "aligner", "tree_builder", "phylogeny_candidates",
                      "phylogeny_root_rule", "phylogeny_pair_rule", "species_tree_rooting"):
            argv += ["--" + field, str(getattr(args, field))]
        if args.species_tree is not None:
            argv += ["--species_tree", str(args.species_tree.resolve())]
    return argv


def configurations(inputs, baseline, output, cpu_count=32):
    if type(cpu_count) is not int or cpu_count < 1:
        raise ValueError("Require a positive integer CPU allocation")
    if inputs["planned_runs"] != planned_runs() or [d["proteomes"] for d in inputs["datasets"]] != [4, 8, 12]:
        raise ValueError("Changed complete scaling inventory or run order")
    datasets = {d["proteomes"]: d for d in inputs["datasets"]}
    core = Path(baseline["core_root"])
    python = Path(baseline["tool_entrypoints"]["orthohmm_python"]["absolute_path"])
    orthofinder = Path(baseline["tool_entrypoints"]["orthofinder"]["absolute_path"])
    runs = []
    for row in inputs["planned_runs"]:
        dataset = datasets[row["proteomes"]]
        directory = output / ("run_%02d" % row["index"])
        all_configs = commands({"input": dataset["input_directory"]}, directory, core, python, orthofinder)
        method = "orthofinder_full" if row["method"] == METHODS[2] else row["method"]
        config = all_configs[method]
        flags = ("-t", "-a") if method == "orthofinder_full" else ("--cpu",)
        for flag in flags:
            if config["argv"].count(flag) != 1:
                raise ValueError("Ambiguous thread option")
            config["argv"][config["argv"].index(flag) + 1] = str(cpu_count)
        native = config["argv"] if method == "orthofinder_full" else native_orthohmm(config)
        runs.append({**row, "native_method": method, "dataset": dataset, "configuration": config,
                     "native_argv": native, "cwd": str(core), "measurement_directory": str(directory / "measurement"),
                     "boundary": "Native CLI launch through exit; includes native metrics/output writing, excludes harness hashing and external validation"})
    return runs


def prepare(root, output, destination):
    if output.exists() or destination.exists():
        raise FileExistsError("Require unused scaling output and manifest paths")
    results = root / "benchmark_tools/results"
    input_path = results / "publication_scaling_inputs_20260916.json"
    baseline_path = results / "publication_variable_native_methods_20260916.json"
    inputs, baseline = read_frozen(input_path, INPUT_SHA), read_frozen(baseline_path, METHOD_SHA)
    verify_environment(baseline)
    _, resolved = execution_environment(baseline)
    for dataset in inputs["datasets"]:
        for item in dataset["inputs"]:
            check(item)
        actual = sorted((record(p) for p in Path(dataset["input_directory"]).iterdir()), key=lambda r: r["path"])
        if actual != dataset["inputs"]:
            raise ValueError("Scaling input directory differs from its frozen inventory")
    runs = configurations(inputs, baseline, output)
    report = {"status": "native_scaling_commands_prepared_unrun", "accuracy_evaluated": False,
              "source": record(__file__), "inputs": record(input_path), "baseline": record(baseline_path),
              "helper_sources": [record(Path(__file__).with_name(n)) for n in
                                 ("benchmark_production.py", "prepare_simulation_methods.py", "prepare_scaling_inputs.py")],
              "environment_overrides": {**baseline["environment_overrides"], "PYTHONUNBUFFERED": "1"},
              "prepend_path": baseline["prepend_path"], "resolved_executables": resolved,
              "resource_plan": inputs["resource_plan"], "runs": runs,
              "measurement_plan": {"resource_interval_s": 1., "host_interval_s": 30., "monitor_host": True},
              "remaining_gates": ["Dedicated task allocation, controlled host window and collector overhead assessment",
                  "Fresh per-run input copies for OrthoFinder outside inference timer; no search reuse",
                  "Native process/output validation independent of the reporting harness, including original .fa inputs",
                  "Before/after source/runtime/input checks; provenance hashing outside inference timer",
                  "No accuracy, comparative timings, or completed scaling runs are established by this manifest"]}
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.output_root.resolve(), args.manifest.resolve())
