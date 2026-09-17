"""Independently admit the complete generating/NNI panel before truth scoring."""

import argparse
from collections import Counter
import csv
import io
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import Phylo
from benchmark_tools.admit_simulation_mode_panel import expected_configuration, PANEL_SHA
from benchmark_tools.assemble_simulation_results import TERMINAL, admit_method, input_universe, verify_scoring_dependencies
from benchmark_tools.audit_mode_partitions import ADMISSION_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_portable_simulation_trees import SOURCE_SHA
from benchmark_tools.prepare_species_tree_robustness import topology
from benchmark_tools.run_simulation_methods import read_frozen, verify_environment, verify_inputs, execution_environment
from benchmark_tools.run_simulation_tree_experiment import PORTABLE_SHA, require_mode_equivalence, select_tree
from benchmark_tools.run_simulation_tree_mode_control import METHOD_SHA, METHODS, native_tree, artifact_inventory, parser_checks
from benchmark_tools.simulation_method_outputs import load_predictions
from benchmark_tools.verify_ygob_validation import require_completed_job

EXECUTOR = "f30eb87f107002e51f1a37b2935982a3a5facb04"
PARTITION_SHA = "bc091032e411ae878275a5965533f1edc45ce1087a54f97dddd867765f9a779d"
HELPERS = ("simulation_supplied_commands.py", "run_simulation_methods.py", "run_simulation_tree_mode_control.py",
           "assemble_simulation_results.py", "validate_simulation_outputs.py", "prepare_portable_simulation_trees.py")


def terminal_panel(accounting):
    rows = list(csv.DictReader(io.StringIO(accounting), delimiter="|"))
    result = []
    for index in range(210):
        job = 21405 if index == 0 else 21406
        name = f"{job}_{index}"
        matches = [r for r in rows if r.get("JobID") == name]
        if len(matches) != 1 or matches[0]["State"].split()[0] not in TERMINAL:
            raise ValueError(f"Cell {name} not uniquely terminal; do not admit a partial panel")
        row = matches[0]
        if row["State"] == "COMPLETED" and row["ExitCode"] != "0:0":
            raise ValueError("Contradictory scheduler completion")
        result.append(row)
    return result


def check_identity(provenance, index, task, tree, configured, executor):
    if (provenance["index"] != index or provenance["tree"] != tree
            or provenance["configured"] != configured or provenance["executor_commit"] != EXECUTOR
            or provenance["job_id"] != task["JobIDRaw"]
            or provenance["array_job_id"] != task["JobID"].split("_")[0]
            or provenance["array_task_id"] != str(index)
            or provenance["source"] != record(executor / "benchmark_tools/run_simulation_tree_experiment.py")
            or provenance["helpers"] != [record(executor / "benchmark_tools" / n) for n in HELPERS]):
        raise ValueError("Tree-cell identity, command or executor differs")


def validate_cell(root, index, task, tree, manifest, generation, executor, gates, mode):
    label = f"{tree['condition']}_{tree['seed']}"
    dataset = next(d for d in manifest["datasets"] if d["label"] == label)
    destination = root / "benchmarks/results/simulation_tree_experiments_v1" / f"cell_{index}"
    common = {"index": index, "label": label, "condition": tree["condition"], "seed": tree["seed"],
              "variant": tree["variant"], "scheduler": task, "tree": tree["tree"]}
    failed = lambda reason: [{**common, "method": m, "status": "failed", "reason": reason} for m in METHODS]
    preflight = destination / "preflight.json"
    if not preflight.exists():
        if task["State"] == "COMPLETED":
            raise ValueError("Completed cell missing preflight")
        return failed("Unsuccessful scheduler task without preflight evidence")
    configured = expected_configuration(dataset, {m: {"tree": tree["tree"]} for m in METHODS}, destination, METHODS)
    provenance = json.loads(preflight.read_text())
    check_identity(provenance, index, task, tree, configured, executor)
    if any(provenance[key] != value for key, value in gates.items()):
        raise ValueError("Tree-cell manifest gates differ")
    if provenance["mode_controls"] != [r for r in mode["records"] if r["label"] == label]:
        raise ValueError("Tree-cell mode controls differ")
    common["preflight"] = record(preflight)
    if task["State"] != "COMPLETED":
        return failed("Unsuccessful terminal scheduler task; no successful inference claim")
    report_path = destination / "results.json"
    report = json.loads(report_path.read_text())
    if (report["status"] != "supplied_tree_inference_complete_pending_independent_admission"
            or report["accuracy_evaluated"] is not False or report["provenance"] != provenance
            or set(report["methods"]) != set(METHODS)
            or report["execution"] != record(destination / "execution/status.json")):
        raise ValueError("Tree-cell completion or execution link differs")
    verified = verify_inputs(dataset, generation, root / "benchmarks/work/publication_variable_simulation_panel_v2",
                             manifest["generation_manifest"]["sha256"])
    if verified["status"] != "ready":
        raise ValueError("Tree-cell inputs not ready")
    owners, species = input_universe(verified["inputs"], json.loads(Path(dataset["truth"]).read_text()))
    env, resolved = execution_environment(manifest)
    if (provenance["resolved_tools"] != resolved or provenance["native_parsers"] !=
            parser_checks(manifest, {m: {"tree": tree["tree"]} for m in METHODS}, species, env)):
        raise ValueError("Resolved tools or native parser evidence differs")
    execution = json.loads((destination / "execution/status.json").read_text())
    if (execution["provenance"] != provenance or execution["verified_inputs"] != verified
            or execution["dataset"] != label or set(execution["methods"]) != set(METHODS)
            or execution["status"] != "finished_pending_native_validation"):
        raise ValueError("Tree-cell execution identity differs")
    common.update(native_report=record(report_path), execution=report["execution"], truth=verified["truth"])
    records = []
    for method in METHODS:
        admission = admit_method(method, configured, execution, verified, manifest)
        observed = report["methods"][method]
        if admission != observed["admission"]:
            raise ValueError("Fresh native admission disagrees with runner")
        row = {**common, "method": method, "status": admission["status"], "admission": admission}
        if admission["status"] == "admitted":
            output = Path(configured["methods"][method]["output"])
            supplied = Phylo.read(tree["tree"]["path"], "newick")
            actual = Phylo.read(native_tree(method, output), "newick")
            retained = topology(actual) == topology(supplied) and {t.name for t in actual.get_terminals()} == set(species)
            inventory = artifact_inventory(method, output)
            if retained != observed["supplied_topology_retained"] or inventory != observed["retained_artifacts"]:
                raise ValueError("Fresh tree/artifact validation disagrees with runner")
            pairs, files = load_predictions(method, output, owners, species)
            inventoried = {r["absolute_path"] for r in execution["methods"][method]["outputs"]}
            if not {str(p) for p in files} <= inventoried:
                raise ValueError("Predictions absent from verified native output inventory")
            row.update(status="admitted" if retained else "supplied_tree_not_retained",
                       supplied_topology_retained=retained, retained_artifacts=inventory,
                       prediction_files=[record(p) for p in files], native_pair_count=len(set(tuple(sorted(p)) for p in pairs)))
        records.append(row)
    return records


def admit(root, output):
    if output.exists():
        raise FileExistsError(output)
    if Path.cwd().resolve() != root:
        raise ValueError("Verify from original repository")
    accounting = subprocess.check_output(["sacct", "-j", "21405,21406", "--parsable2",
                                         "--format=JobID,JobIDRaw,State,ExitCode,Elapsed"], text=True)
    tasks = terminal_panel(accounting)
    results = root / "benchmark_tools/results"
    paths = {"mode_admission": root / "benchmarks/results/simulation_mode_panel_admission_v1/results.json",
             "portable_manifest": results / "simulation_portable_trees_prepared_20260917.json",
             "prepared_manifest": results / "simulation_tree_controls_prepared_20260917.json",
             "method_manifest": results / "publication_variable_native_methods_20260916.json"}
    hashes = dict(zip(paths, (ADMISSION_SHA, PORTABLE_SHA, SOURCE_SHA, METHOD_SHA)))
    data = {key: read_frozen(path, hashes[key]) for key, path in paths.items()}
    mode, manifest = data["mode_admission"], data["method_manifest"]
    for item in [mode["source"], *mode["helper_sources"], mode["fresh_pilot_admission"]]:
        check(item)
    mode_accounting = subprocess.check_output(["sacct", "-j", "21367", "--parsable2",
                                               "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    mode_scheduler = require_completed_job(mode_accounting, 21367)
    panel = read_frozen(results / "simulation_mode_panel_prepared_20260917.json", PANEL_SHA)
    require_mode_equivalence(mode, panel)
    partitions = read_frozen(results / "simulation_mode_partitions_verified_20260917.json", PARTITION_SHA)
    check(partitions["source"])
    check(partitions["mode_admission"])
    if not partitions["all_identical"] or len(partitions["checks"]) != 199 or not all(r["identical"] for r in partitions["checks"]):
        raise ValueError("Downstream mode partition gate failed")
    for row in partitions["checks"]:
        check(row["before"])
        check(row["after"])
    executor = root / "benchmarks/work/publication_simulation_tree_experiments_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR:
        raise ValueError("Frozen executor revision differs")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    verify_environment(manifest)
    verify_scoring_dependencies(manifest)
    gen = manifest["generation_manifest"]
    generation = read_frozen(Path(gen["absolute_path"]), gen["sha256"])
    gates = {key: record(path) for key, path in paths.items()}
    cell_gates = {**gates, "mode_admission_scheduler": mode_scheduler}
    output.mkdir(parents=True)
    report = {"status": "validating", "accuracy_evaluated": False, "source": record(__file__),
              "gates": gates, "accounting": accounting, "records": [],
              "helper_sources": [record(Path(__file__).with_name(n)) for n in
                  (*HELPERS, "admit_simulation_mode_panel.py", "audit_mode_partitions.py", "run_simulation_tree_experiment.py",
                   "simulation_method_outputs.py", "prepare_species_tree_robustness.py")]}
    try:
        for index, task in enumerate(tasks):
            tree = select_tree(data["portable_manifest"], data["prepared_manifest"], index)
            report["records"].extend(validate_cell(root, index, task, tree, manifest, generation, executor, cell_gates, mode))
            (output / "progress.json").write_text(json.dumps({"cells_validated": index + 1}) + "\n")
        expected = {(i, m) for i in range(210) for m in METHODS}
        if len(report["records"]) != 420 or {(r["index"], r["method"]) for r in report["records"]} != expected:
            raise ValueError("Incomplete 420-method outcome inventory")
        verify_environment(manifest)
        for item in [report["source"], *report["helper_sources"], *gates.values()]:
            check(item)
        report.update(status="tree_panel_verified_unscored", counts=dict(Counter(r["status"] for r in report["records"])),
                      limitations=["No truth scoring or isolated runtime claim.",
                                   "Cross-arm upstream equivalence still requires comparison before tree-only causal claims."])
    except Exception as error:
        report.update(status="validation_failed", error=str(error))
        raise
    finally:
        (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.output.resolve())
