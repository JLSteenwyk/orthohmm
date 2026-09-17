"""Execute a fixed generating/NNI tree cell only after mode-control admission."""

import argparse
from copy import deepcopy
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import Phylo
from benchmark_tools.admit_simulation_mode_panel import PANEL_SHA
from benchmark_tools.assemble_simulation_results import admit_method, input_universe
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_portable_simulation_trees import branch_map, SOURCE_SHA
from benchmark_tools.prepare_species_tree_robustness import topology
from benchmark_tools.run_simulation_methods import read_frozen, verify_environment, verify_inputs, execution_environment, execute
from benchmark_tools.run_simulation_tree_mode_control import METHOD_SHA, METHODS, native_tree, parser_checks, artifact_inventory
from benchmark_tools.simulation_supplied_commands import fresh_supplied_method
from benchmark_tools.verify_ygob_validation import require_completed_job

PORTABLE_SHA = "b9ed4fb8dc27da28dd56c674d1ece2edbb3a04697dec14ea3d4538bf6d9dbc0b"
VARIANTS = ("generating", "nni1", "nni2")


def require_mode_equivalence(admission, panel):
    if (admission["status"] != "mode_panel_verified_unscored" or admission["accuracy_evaluated"] is not False
            or admission["panel"]["sha256"] != PANEL_SHA):
        raise ValueError("Mode panel is not independently admitted")
    expected = {(row["label"], method): "unavailable" if method in row["unavailable"] else "equivalent"
                for row in panel["rows"] for method in METHODS}
    records = admission["records"]
    observed = {(row["label"], row["method"]): row["status"] for row in records}
    if len(records) != 140 or len(expected) != 140 or observed != expected:
        raise ValueError("Review mode failures or non-equivalence before authorizing tree experiments")
    for row in records:
        if row["status"] == "unavailable":
            original = next(p for p in panel["rows"] if p["label"] == row["label"])
            if row["baseline_failure"] != original["unavailable"][row["method"]]:
                raise ValueError("Unavailable baseline evidence changed")


def select_tree(portable, prepared, index):
    rows = portable["trees"]
    if isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < 210 or len(rows) != 210:
        raise ValueError("Invalid fixed tree-cell index or inventory")
    expected = {(d["condition"], d["seed"], v["label"]) for d in prepared["datasets"] for v in d["variants"]}
    observed = {(r["condition"], r["seed"], r["variant"]) for r in rows}
    if len(expected) != 210 or expected != observed:
        raise ValueError("Incomplete or duplicated portable tree inventory")
    row = rows[index]
    if row["variant"] not in VARIANTS:
        raise ValueError("Unknown fixed tree variant")
    dataset = next(d for d in prepared["datasets"] if (d["condition"], d["seed"]) == (row["condition"], row["seed"]))
    variant = next(v for v in dataset["variants"] if v["label"] == row["variant"])
    if row["source_tree"] != variant["tree"] or row["taxa"] != dataset["taxa"]:
        raise ValueError("Portable tree source or taxa changed")
    for item in (row["source_tree"], row["tree"]):
        check(item)
    original, actual = [Phylo.read(row[key]["path"], "newick") for key in ("source_tree", "tree")]
    if (topology(original) != topology(actual) or branch_map(original) != branch_map(actual)
            or sorted(t.name for t in actual.get_terminals()) != row["taxa"]):
        raise ValueError("Portable supplied tree differs from frozen topology, lengths or taxa")
    return row


def run(root, index, admission_sha, output):
    if output.exists():
        raise FileExistsError(output)
    if Path.cwd().resolve() != root:
        raise ValueError("Run verification from original repository")
    results = root / "benchmark_tools/results"
    admission_path = root / "benchmarks/results/simulation_mode_panel_admission_v1/results.json"
    admission = read_frozen(admission_path, admission_sha)
    mode_panel = read_frozen(results / "simulation_mode_panel_prepared_20260917.json", PANEL_SHA)
    require_mode_equivalence(admission, mode_panel)
    accounting = subprocess.check_output(["sacct", "-j", "21367", "--parsable2", "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 21367)
    auditor = root / "benchmarks/work/publication_simulation_mode_admission_v1"
    if admission["source"] != record(auditor / "benchmark_tools/admit_simulation_mode_panel.py"):
        raise ValueError("Unexpected mode-panel auditor")
    for item in [admission["source"], *admission["helper_sources"], admission["fresh_pilot_admission"]]:
        check(item)
    portable_path = results / "simulation_portable_trees_prepared_20260917.json"
    prepared_path = results / "simulation_tree_controls_prepared_20260917.json"
    portable, prepared = read_frozen(portable_path, PORTABLE_SHA), read_frozen(prepared_path, SOURCE_SHA)
    tree = select_tree(portable, prepared, index)
    label = f"{tree['condition']}_{tree['seed']}"
    manifest_path = results / "publication_variable_native_methods_20260916.json"
    manifest = read_frozen(manifest_path, METHOD_SHA)
    dataset = next(d for d in manifest["datasets"] if d["label"] == label)
    gen = manifest["generation_manifest"]
    generation = read_frozen(Path(gen["absolute_path"]), gen["sha256"])
    verify_environment(manifest)
    verified = verify_inputs(dataset, generation, root / "benchmarks/work/publication_variable_simulation_panel_v2", gen["sha256"])
    if verified["status"] != "ready" or sorted(Path(r["absolute_path"]).stem for r in verified["inputs"]) != tree["taxa"]:
        raise ValueError("Supplied-tree input universe differs")
    owners, species = input_universe(verified["inputs"], json.loads(Path(dataset["truth"]).read_text()))
    controls = [r for r in admission["records"] if r["label"] == label]
    for control in controls:
        for key in ("task_evidence", "native_report"):
            if key in control:
                check(control[key])
        for item in control.get("prediction_files", []):
            check(item)
    configured = deepcopy(dataset)
    configured["methods"] = {m: fresh_supplied_method(m, dataset["methods"][m], tree["tree"]["path"], output / m) for m in METHODS}
    env, resolved = execution_environment(manifest)
    parsers = parser_checks(manifest, {m: {"tree": tree["tree"]} for m in METHODS}, species, env)
    provenance = {"source": record(__file__), "index": index, "mode_admission": record(admission_path),
                  "mode_admission_scheduler": scheduler, "mode_controls": controls, "tree": tree,
                  "portable_manifest": record(portable_path), "prepared_manifest": record(prepared_path),
                  "method_manifest": record(manifest_path), "configured": configured, "resolved_tools": resolved,
                  "native_parsers": parsers, "job_id": os.environ.get("SLURM_JOB_ID"),
                  "array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"), "array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
                  "executor_commit": subprocess.check_output(["git", "-C", str(Path(__file__).resolve().parent.parent),
                                                               "rev-parse", "HEAD"], text=True).strip(),
                  "helpers": [record(Path(__file__).with_name(name)) for name in
                     ("simulation_supplied_commands.py", "run_simulation_methods.py", "run_simulation_tree_mode_control.py",
                      "assemble_simulation_results.py", "validate_simulation_outputs.py", "prepare_portable_simulation_trees.py")],
                  "scope": "Supplied generating/NNI tree, fresh inference; both tools retained even without original inferred baseline"}
    output.mkdir(parents=True)
    (output / "preflight.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    report = {"status": "running", "provenance": provenance, "accuracy_evaluated": False, "methods": {}}
    try:
        execution = execute(configured, list(METHODS), env, output / "execution", verified, provenance, manifest)
        report["execution"] = record(output / "execution/status.json")
        for method in METHODS:
            native = admit_method(method, configured, execution, verified, manifest)
            row = {"admission": native}
            report["methods"][method] = row
            if native["status"] != "admitted":
                continue
            inferred = Phylo.read(native_tree(method, output / method), "newick")
            supplied = Phylo.read(tree["tree"]["path"], "newick")
            row["supplied_topology_retained"] = (topology(inferred) == topology(supplied)
                and {t.name for t in inferred.get_terminals()} == set(species))
            row["retained_artifacts"] = artifact_inventory(method, output / method)
        verify_environment(manifest)
        for item in [provenance["source"], *provenance["helpers"], tree["tree"], tree["source_tree"]]:
            check(item)
        for path, digest in ((admission_path, admission_sha), (portable_path, PORTABLE_SHA),
                             (prepared_path, SOURCE_SHA), (manifest_path, METHOD_SHA)):
            read_frozen(path, digest)
        if verify_inputs(dataset, generation, root / "benchmarks/work/publication_variable_simulation_panel_v2", gen["sha256"]) != verified:
            raise ValueError("Inputs changed during supplied-tree inference")
        report["status"] = "supplied_tree_inference_complete_pending_independent_admission"
    except Exception as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        (output / "results.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--mode-admission-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.root.resolve(), args.index, args.mode_admission_sha256, args.output.resolve())
