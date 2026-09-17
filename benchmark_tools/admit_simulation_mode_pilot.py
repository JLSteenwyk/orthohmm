"""Independently reread and admit only the fixed supplied-own-tree pilot."""

import argparse
from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import Phylo
from benchmark_tools.run_simulation_tree_mode_control import (
    METHOD_SHA, RESULT_SHA, METHODS, baseline_records, native_tree, artifact_inventory,
    compare_inventory, compare_pairs,
)
from benchmark_tools.assemble_simulation_results import admit_method, input_universe
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_species_tree_robustness import topology
from benchmark_tools.run_simulation_methods import read_frozen, verify_environment, verify_inputs
from benchmark_tools.simulation_method_outputs import load_predictions
from benchmark_tools.simulation_supplied_commands import fresh_supplied_method
from benchmark_tools.verify_ygob_validation import require_completed_job

PILOT_SHA = "446a3c19cebf87fe562cc9e88cf183986e4e644b71116eaa5ba2839c9083f561"
EXECUTOR = "0e797bab52f91e494c1e3be4b728de3dbedd3fe8"


def retained_equivalence(method, before, after):
    comparison = compare_inventory(before, after)
    if method == "orthohmm_satellite_v2":
        return comparison["identical"]
    if method != "orthofinder_full":
        raise ValueError("Unsupported method")
    if comparison["extra"] or set(comparison["missing"]) - {"Alignments_ids/SpeciesTreeAlignment.fa"}:
        return False
    allowed = {"clusters_OrthoFinder_I1.2.txt", "clusters_OrthoFinder_I1.2.txt_id_pairs.txt"}
    if set(comparison["changed"]) - allowed:
        return False
    for key in comparison["changed"]:
        contents = []
        for inventory in (before, after):
            check(inventory[key])
            lines = Path(inventory[key]["path"]).read_text().splitlines()
            comments = [line for line in lines if line.startswith("# cline: ")]
            if len(comments) != 1:
                return False
            contents.append([line for line in lines if not line.startswith("# cline: ")])
        if contents[0] != contents[1]:
            return False
    return True


def admit(root, output):
    if output.exists():
        raise FileExistsError(output)
    if Path.cwd().resolve() != root:
        raise ValueError("Run environment checks from original repository")
    results = root / "benchmark_tools/results"
    path = results / "simulation_mode_control_baseline_seed1_20260917.json"
    report = read_frozen(path, PILOT_SHA)
    provenance = report["provenance"]
    if (report["status"] != "complete_pending_independent_admission" or report["accuracy_evaluated"] is not False
            or provenance["job_id"] != "21332" or provenance["executor_commit"] != EXECUTOR
            or set(report["methods"]) != set(METHODS)):
        raise ValueError("Unexpected fixed pilot identity or scope")
    directory = root / "benchmarks/results/simulation_mode_control_baseline_seed1_v1"
    if json.loads((directory / "results.json").read_text()) != report:
        raise ValueError("Snapshot differs from original report")
    if json.loads((directory / "preflight.json").read_text()) != provenance:
        raise ValueError("Pilot preflight changed")
    executor = root / "benchmarks/work/publication_simulation_mode_control_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR:
        raise ValueError("Pilot executor changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    expected_sources = [record(executor / "benchmark_tools" / name) for name in
        ("run_simulation_tree_mode_control.py", "simulation_supplied_commands.py", "run_simulation_methods.py",
         "assemble_simulation_results.py", "validate_simulation_outputs.py", "simulation_method_outputs.py")]
    if provenance["sources"] != expected_sources or provenance["source"] != expected_sources[0]:
        raise ValueError("Pilot executor source inventory differs")
    accounting = subprocess.check_output(["sacct", "-j", "21332", "--parsable2", "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 21332)
    manifest = read_frozen(results / "publication_variable_native_methods_20260916.json", METHOD_SHA)
    evidence = read_frozen(results / "simulation_variable_native_results_20260916.json", RESULT_SHA)
    for item in (provenance["manifest"], provenance["baseline_results"], report["execution"]):
        check(item)
    verify_environment(manifest)
    dataset = next(d for d in manifest["datasets"] if d["label"] == "baseline_20261101")
    gen = manifest["generation_manifest"]
    verified = verify_inputs(dataset, read_frozen(Path(gen["absolute_path"]), gen["sha256"]),
                             root / "benchmarks/work/publication_variable_simulation_panel_v2", gen["sha256"])
    baseline = baseline_records(dataset, evidence, manifest, verified)
    if baseline != provenance["baseline"]:
        raise ValueError("Baseline changed since pilot")
    # Reconstruct the exact fresh-run commands without invoking a builder that refuses existing outputs.
    configured = deepcopy(dataset)
    configured["methods"] = {}
    for method in METHODS:
        original = deepcopy(dataset["methods"][method])
        destination = directory / method
        argv = original["argv"]
        if method == "orthohmm_satellite_v2":
            metrics = str(directory / (method + ".json"))
            argv[3:5] = [str(destination), metrics]
            argv[argv.index("--species-tree-mode") + 1] = "supplied"
            argv.extend(["--species-tree", baseline[method]["tree"]["path"]])
            original["metrics"] = metrics
        else:
            argv[argv.index("-f") + 1] = str(destination / "input")
            argv.extend(["-s", baseline[method]["tree"]["path"]])
            original["copy_inputs_to"] = str(destination / "input")
        original.update(output=str(destination), supplied_tree=baseline[method]["tree"]["path"],
                        execution_scope="fresh supplied-tree inference; no upstream cache reuse")
        configured["methods"][method] = original
    if configured != provenance["configured"]:
        raise ValueError("Pilot changed non-tree configuration")
    execution = json.loads(Path(report["execution"]["path"]).read_text())
    if (execution["provenance"] != provenance or execution["verified_inputs"] != verified
            or execution["dataset"] != dataset["label"] or execution["failed_methods"]):
        raise ValueError("Pilot execution inputs or provenance differ")
    owners, species = input_universe(verified["inputs"], json.loads(Path(dataset["truth"]).read_text()))
    summaries = {}
    for method in METHODS:
        admission = admit_method(method, configured, execution, verified, manifest)
        before, after = Path(dataset["methods"][method]["output"]), Path(configured["methods"][method]["output"])
        a, _ = load_predictions(method, before, owners, species)
        b, _ = load_predictions(method, after, owners, species)
        pairs = compare_pairs(a, b)
        ta, tb = Phylo.read(native_tree(method, before), "newick"), Phylo.read(native_tree(method, after), "newick")
        same_tree = topology(ta) == topology(tb) and {t.name for t in ta.get_terminals()} == {t.name for t in tb.get_terminals()}
        original, retained = baseline[method]["retained_artifacts"], artifact_inventory(method, after)
        row = report["methods"][method]
        if (admission != row["admission"] or admission["status"] != "admitted" or pairs != row["pairs"]
                or same_tree != row["rooted_tree_identical"] or retained != row["retained_artifacts"]
                or compare_inventory(original, retained) != row["artifact_comparison"]):
            raise ValueError("Independent pilot reread disagrees with report")
        if not pairs["identical"] or not same_tree or not retained_equivalence(method, original, retained):
            raise ValueError("Supplied-tree pilot is not equivalent")
        summaries[method] = {"pairs": len({tuple(sorted(p)) for p in b}), "rooted_tree_identical": True,
                             "retained_artifacts_equivalent": True, "byte_comparison": row["artifact_comparison"]}
    verify_environment(manifest)
    read_frozen(path, PILOT_SHA)
    output.mkdir(parents=True)
    result = {"status": "pilot_equivalence_verified", "accuracy_evaluated": False, "publication_ready": False,
              "source": record(__file__), "pilot_report": record(path), "scheduler": scheduler,
              "accounting": accounting, "methods": summaries, "executor_sources": expected_sources,
              "scope": "Only baseline_20261101 and the retained native artifacts; not all-dataset equivalence"}
    (output / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.output.resolve())
