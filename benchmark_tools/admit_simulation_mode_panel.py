"""Validate the entire terminal mode-control array without scoring evolutionary truth."""

import argparse
from collections import Counter
from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import Phylo
from benchmark_tools.admit_simulation_mode_pilot import retained_equivalence, admit as admit_pilot
from benchmark_tools.assemble_simulation_results import terminal_tasks, admit_method, input_universe
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_simulation_mode_panel import panel_rows
from benchmark_tools.prepare_species_tree_robustness import topology
from benchmark_tools.run_simulation_methods import read_frozen, verify_inputs, verify_environment
from benchmark_tools.run_simulation_tree_mode_control import (
    METHOD_SHA, RESULT_SHA, METHODS, baseline_records, native_tree, artifact_inventory, compare_inventory, compare_pairs,
)
from benchmark_tools.simulation_method_outputs import load_predictions

PANEL_SHA = "03f0d115d0f2a1f5354ca6a8a369e419dbd9e10945e450dfd2534d61212bb747"
EXECUTOR = "2d9536966c905868b5b3e168d17b9b9f46d5eae5"
JOB = 21334


def verify_fresh_pilot(fresh, prior):
    if ({k: v for k, v in fresh.items() if k != "source"} != {k: v for k, v in prior.items() if k != "source"}
            or any(fresh["source"][key] != prior["source"][key] for key in ("bytes", "sha256"))):
        raise ValueError("Fresh pilot admission differs from the frozen admission")
    check(fresh["source"])


def expected_configuration(dataset, baseline, destination, methods):
    result = deepcopy(dataset)
    result["methods"] = {}
    for method in methods:
        original = deepcopy(dataset["methods"][method])
        output = destination / method
        tree = baseline[method]["tree"]["path"]
        argv = original["argv"]
        if method == "orthohmm_satellite_v2":
            metrics = str(destination / (method + ".json"))
            argv[3:5] = [str(output), metrics]
            argv[argv.index("--species-tree-mode") + 1] = "supplied"
            argv.extend(["--species-tree", tree])
            original["metrics"] = metrics
        elif method == "orthofinder_full":
            argv[argv.index("-f") + 1] = str(output / "input")
            argv.extend(["-s", tree])
            original["copy_inputs_to"] = str(output / "input")
        else:
            raise ValueError("Unknown mode-control method")
        original.update(output=str(output), supplied_tree=tree,
                        execution_scope="fresh supplied-tree inference; no upstream cache reuse")
        result["methods"][method] = original
    return result


def task_identity(status, row, scheduler, panel_record, executor):
    if (status["row"] != row or status["panel"] != panel_record or status["accuracy_evaluated"] is not False
            or status["job_id"] != scheduler["JobIDRaw"] or status["array_job_id"] != str(JOB)
            or status["array_task_id"] != str(row["index"])
            or status["source"] != record(executor / "benchmark_tools/run_simulation_mode_panel.py")):
        raise ValueError("Mode-control task identity differs")


def compare_method(method, dataset, configured, execution, verified, manifest, baseline, observed, owners, species):
    admission = admit_method(method, configured, execution, verified, manifest)
    if admission != observed["admission"]:
        raise ValueError("Independent native admission disagrees with runner")
    if admission["status"] != "admitted":
        return {"status": "failed", "admission": admission}
    before, after = Path(dataset["methods"][method]["output"]), Path(configured["methods"][method]["output"])
    old, old_files = load_predictions(method, before, owners, species)
    new, new_files = load_predictions(method, after, owners, species)
    pairs = compare_pairs(old, new)
    ta, tb = Phylo.read(native_tree(method, before), "newick"), Phylo.read(native_tree(method, after), "newick")
    tree_equal = topology(ta) == topology(tb) and {t.name for t in ta.get_terminals()} == {t.name for t in tb.get_terminals()}
    prior, current = baseline[method]["retained_artifacts"], artifact_inventory(method, after)
    comparison = compare_inventory(prior, current)
    if (pairs != observed["pairs"] or tree_equal != observed["rooted_tree_identical"]
            or current != observed["retained_artifacts"] or comparison != observed["artifact_comparison"]):
        raise ValueError("Independent mode comparison disagrees with runner")
    artifacts_equal = retained_equivalence(method, prior, current)
    equivalent = pairs["identical"] and tree_equal and artifacts_equal
    return {"status": "equivalent" if equivalent else "not_equivalent", "admission": admission,
            "pairs": pairs, "rooted_tree_identical": tree_equal, "retained_artifacts_equivalent": artifacts_equal,
            "artifact_comparison": comparison, "prediction_files": [record(p) for p in [*old_files, *new_files]]}


def validate_dataset(root, row, scheduler, panel, panel_record, manifest, evidence, generation, executor):
    directory = root / "benchmarks/results/simulation_mode_panel_v1"
    task_path = directory / "tasks" / row["label"] / "status.json"
    common = {"label": row["label"], "condition": row["condition"], "seed": row["seed"], "scheduler": scheduler}
    summaries = [{**common, "method": method, "status": "unavailable", "baseline_failure": reason}
                 for method, reason in row["unavailable"].items()]
    failed = lambda reason: [{**common, "method": method, "status": "failed", "reason": reason}
                             for method in row["methods"]]
    if not task_path.is_file():
        if scheduler["State"] == "COMPLETED":
            raise ValueError("Completed task has no status evidence")
        return summaries + failed("Unsuccessful terminal scheduler task without retained status")
    status = json.loads(task_path.read_text())
    task_identity(status, row, scheduler, panel_record, executor)
    common["task_evidence"] = record(task_path)
    if status["status"] in {"failed", "running"} or scheduler["State"] != "COMPLETED":
        if scheduler["State"] == "COMPLETED":
            raise ValueError("Scheduler success contradicts failed or running task status")
        return summaries + failed(status.get("error", "Interrupted mode-control execution"))
    if row["reuse_pilot"]:
        if status["status"] != "reused_admitted_pilot" or status["admission"] != panel["pilot_admission"]:
            raise ValueError("Pilot reuse identity differs")
        pilot = json.loads(Path(panel["pilot_admission"]["path"]).read_text())
        check(pilot["pilot_report"])
        for method in row["methods"]:
            summaries.append({**common, "method": method, "status": "equivalent", "reused_pilot": True,
                              "admission": panel["pilot_admission"], "pilot_summary": pilot["methods"][method]})
        return summaries
    if not row["methods"]:
        if status["status"] != "no_available_inferred_baseline":
            raise ValueError("Unexpected unavailable-only task state")
        return summaries
    if status["status"] != "executed_pending_independent_admission":
        raise ValueError("Unexpected completed mode-control task state")
    destination = directory / row["label"]
    if status["native_result"] != record(destination / "results.json"):
        raise ValueError("Native report link differs")
    report = json.loads((destination / "results.json").read_text())
    provenance = report["provenance"]
    expected_sources = [record(executor / "benchmark_tools" / name) for name in
        ("run_simulation_tree_mode_control.py", "simulation_supplied_commands.py", "run_simulation_methods.py",
         "assemble_simulation_results.py", "validate_simulation_outputs.py", "simulation_method_outputs.py")]
    if (report["status"] != "complete_pending_independent_admission" or report["status"] != status["native_status"]
            or report["accuracy_evaluated"] is not False or set(report["methods"]) != set(row["methods"])
            or provenance["executor_commit"] != EXECUTOR or provenance["job_id"] != scheduler["JobIDRaw"]
            or provenance["sources"] != expected_sources or provenance["source"] != expected_sources[0]
            or provenance["manifest"] != panel["method_manifest"] or provenance["baseline_results"] != panel["baseline_results"]
            or json.loads((destination / "preflight.json").read_text()) != provenance):
        raise ValueError("Native mode-control provenance differs")
    dataset = next(d for d in manifest["datasets"] if d["label"] == row["label"])
    verified = verify_inputs(dataset, generation, root / "benchmarks/work/publication_variable_simulation_panel_v2",
                             manifest["generation_manifest"]["sha256"])
    baseline = baseline_records(dataset, evidence, manifest, verified, row["methods"])
    configured = expected_configuration(dataset, baseline, destination, row["methods"])
    if baseline != provenance["baseline"] or configured != provenance["configured"]:
        raise ValueError("Mode-control baseline or non-tree configuration changed")
    if report["execution"] != record(destination / "execution/status.json"):
        raise ValueError("Execution status link differs")
    execution = json.loads((destination / "execution/status.json").read_text())
    if (execution["provenance"] != provenance or execution["verified_inputs"] != verified
            or execution["dataset"] != row["label"] or set(execution["methods"]) != set(row["methods"])
            or execution["status"] != "finished_pending_native_validation"):
        raise ValueError("Execution identity or completion state differs")
    owners, species = input_universe(verified["inputs"], json.loads(Path(dataset["truth"]).read_text()))
    for method in row["methods"]:
        summaries.append({**common, "method": method, "native_report": record(destination / "results.json"),
            **compare_method(method, dataset, configured, execution, verified, manifest, baseline,
                             report["methods"][method], owners, species)})
    return summaries


def admit(root, output):
    if output.exists():
        raise FileExistsError(output)
    if Path.cwd().resolve() != root:
        raise ValueError("Verify environment from original repository directory")
    accounting = subprocess.check_output(["sacct", "-j", str(JOB), "--parsable2", "--format=JobID,JobIDRaw,State,ExitCode,Elapsed"], text=True)
    tasks = terminal_tasks(accounting, JOB, 70)
    panel_path = root / "benchmark_tools/results/simulation_mode_panel_prepared_20260917.json"
    panel = read_frozen(panel_path, PANEL_SHA)
    panel_record = record(panel_path)
    for key in ("source", "method_manifest", "baseline_results", "pilot_admission"):
        check(panel[key])
    pilot = json.loads(Path(panel["pilot_admission"]["path"]).read_text())
    if pilot["status"] != "pilot_equivalence_verified" or pilot["accuracy_evaluated"] is not False:
        raise ValueError("Pilot not admitted")
    for item in [pilot["source"], pilot["pilot_report"], *pilot["executor_sources"]]:
        check(item)
    executor = root / "benchmarks/work/publication_simulation_mode_panel_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR:
        raise ValueError("Executor revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    manifest = read_frozen(Path(panel["method_manifest"]["path"]), METHOD_SHA)
    evidence = read_frozen(Path(panel["baseline_results"]["path"]), RESULT_SHA)
    rows = panel_rows(manifest, evidence)
    if rows != panel["rows"]:
        raise ValueError("Frozen mode-control inventory differs")
    verify_environment(manifest)
    gen = manifest["generation_manifest"]
    generation = read_frozen(Path(gen["absolute_path"]), gen["sha256"])
    output.mkdir(parents=True)
    summaries = []
    report = {"status": "validating", "accuracy_evaluated": False, "publication_ready": False,
              "source": record(__file__), "panel": panel_record, "accounting": accounting, "records": summaries,
              "helper_sources": [record(Path(__file__).with_name(name)) for name in
                  ("admit_simulation_mode_pilot.py", "assemble_simulation_results.py", "prepare_simulation_mode_panel.py",
                   "run_simulation_methods.py", "run_simulation_tree_mode_control.py", "simulation_method_outputs.py",
                   "validate_simulation_outputs.py", "prepare_species_tree_robustness.py")]}
    try:
        verify_fresh_pilot(admit_pilot(root, output / "pilot_recheck"), pilot)
        report["fresh_pilot_admission"] = record(output / "pilot_recheck/results.json")
        for row, scheduler in zip(rows, tasks):
            summaries.extend(validate_dataset(root, row, scheduler, panel, panel_record, manifest, evidence, generation, executor))
            (output / "progress.json").write_text(json.dumps({"datasets_validated": row["index"] + 1,
                                                             "method_records": len(summaries)}) + "\n")
        expected = {(r["label"], m) for r in rows for m in METHODS}
        if len(summaries) != 140 or {(r["label"], r["method"]) for r in summaries} != expected:
            raise ValueError("Incomplete all-method outcome inventory")
        verify_environment(manifest)
        read_frozen(panel_path, PANEL_SHA)
        for item in [report["source"], *report["helper_sources"]]:
            check(item)
        report.update(status="mode_panel_verified_unscored", counts=dict(Counter(r["status"] for r in summaries)),
                      limitations=["No evolutionary-truth score or matched timing claim.",
                                   "Non-equivalence and unavailable controls remain explicit; no blanket tree-only interpretation."])
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
