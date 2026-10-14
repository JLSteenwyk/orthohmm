"""Evaluate the shared-hit intervention against admitted canonical controls."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores, replay_native
from benchmark_tools.run_full_search_hillclimb import native_command, runtime, validate_native, verify_input_inventory
from benchmark_tools.run_membership_hillclimb import (
    DATASETS, score_orthobench, score_qfo, score_three_kingdoms, write_json,
)
from benchmark_tools.run_search_recovery import host_activity
from benchmark_tools.run_three_kingdoms_membership_supplement import sources, timed_stage
from benchmark_tools.validate_profile_runtime import require_profile_runtime

PLAN = "hillclimb_shared_hit_pipeline_20261014.json"


def selected_control(report, dataset):
    rows = [row for row in report["arms"] if row["dataset"] == dataset and row["label"] == "canonical"]
    if len(rows) != 1 or rows[0]["status"] != "independent_content_readback_verified":
        raise ValueError("Canonical control unavailable or not independently verified")
    endpoint_vector(dataset, rows[0]["scores"])
    return rows[0]


def prepare(repo, source, output):
    if output.exists() or source == repo:
        raise ValueError("Require new output and private frozen source")
    control_root = repo / "benchmarks/results/canonical_downstream_r1_v1_20261014"
    control_path = control_root / "execution_plan.json"
    control = json.loads(control_path.read_text())
    if runtime() != control["runtime"]:
        raise ValueError("Runtime differs from canonical control")
    plan_path = source / "benchmark_tools" / PLAN
    study = json.loads(plan_path.read_text())
    if tuple(study["datasets"]) != DATASETS:
        raise ValueError("Wrong prespecified dataset inventory")
    evidence = repo / study["selection_evidence"]
    inputs = [record(plan_path), record(evidence), record(control_path),
              *control["inputs"], *control["sources"], *control["runtime"]["files"]]
    for name in ("hmm_viterbi.so", "kmer_prefilter.so", "pair_align.so"):
        suffix = "orthohmm/search/csrc/" + name
        pins = [pin for pin in control["sources"] if pin["path"].endswith(suffix)]
        if len(pins) != 1 or record(source / suffix)["sha256"] != pins[0]["sha256"]:
            raise ValueError("Variant CPU libraries differ from canonical control")
    source_pins = sources(source) + [record(source / "benchmark_tools/run_shared_hit_pipeline.slurm")]
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip()
    for pin in source_pins:
        relative = Path(pin["path"]).relative_to(source)
        if relative.suffix == ".so":
            continue
        blob = subprocess.check_output(["git", "show", revision + ":" + str(relative)], cwd=repo)
        if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
            raise ValueError("Frozen source differs from committed revision")
    cases = []
    for old in control["cases"]:
        folder = output / old["dataset"] / "shared_hit1"
        command = native_command(source, Path(old["fasta_directory"]), folder,
                                 [old["aligner"], old["tree_builder"]], 100, 16)
        command.extend(["--search-result-order", "canonical", "--search-min-total-hits", "1"])
        cases.append({"dataset": old["dataset"], "label": "shared_hit1", "output": str(folder),
                      "fasta": old["fasta_directory"], "inputs": old["fastas"], "golden": old["golden"],
                      "search_max_candidates_per_query": 100, "search_min_total_hits": 1, "command": command})
    for pin in (*inputs, *source_pins):
        check(pin)
    result = {"status": "prepared_waiting_for_control_task", "repo": str(repo), "source_root": str(source),
              "source_revision": revision, "output": str(output), "control_root": str(control_root),
              "inputs": inputs, "sources": source_pins, "runtime": runtime(), "cpu": 16,
              "cases": cases, "tool_path": os.environ["PATH"], "profile_runtime": require_profile_runtime(source),
              "control_dependency_job": "24240", "scope": study["timing_scope"], "defaults_promoted": False}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", result)
    return result


def run(plan_path):
    plan = json.loads(plan_path.read_text())
    source, repo, output = (Path(plan[key]) for key in ("source_root", "repo", "output"))
    if Path(__file__).resolve() != source / "benchmark_tools/run_shared_hit_pipeline.py" or (output / "results.json").exists():
        raise ValueError("Require frozen runner and preserve every attempt")
    accounting = subprocess.check_output(["sacct", "-j", plan["control_dependency_job"], "-X", "-n", "-P",
                                         "-o", "JobID,State%30,ExitCode"], text=True)
    rows = [line.split("|") for line in accounting.splitlines() if line]
    terminal = {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED", "BOOT_FAIL"}
    if len(rows) != 1 or rows[0][0] != plan["control_dependency_job"] or rows[0][1].split()[0] not in terminal:
        raise ValueError("Control task must be authoritatively terminal before execution")
    control_path = Path(plan["control_root"]) / "results.json"
    control_report = json.loads(control_path.read_text())
    checked = [record(plan_path), record(control_path), *plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]
    for pin in checked:
        check(pin)
    if runtime() != plan["runtime"] or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Changed runtime or CUDA visibility")
    result = {"status": "running", "execution_plan": record(plan_path), "control_report": record(control_path),
              "control_accounting": accounting, "source_revision": plan["source_revision"],
              "job_id": os.environ.get("SLURM_JOB_ID"), "arms": [], "scope": plan["scope"],
              "host_before": host_activity(), "defaults_promoted": False, "goal_completed": False}
    write_json(output / "results.json", result)
    env = {**os.environ, "PATH": plan["tool_path"], "PYTHONPATH": str(source), "PYTHONHASHSEED": "0",
           "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "16"}
    nonce = hashlib.sha256(str(output).encode()).hexdigest()[:16]
    control_nonce = hashlib.sha256(plan["control_root"].encode()).hexdigest()[:16]
    for case in plan["cases"]:
        dataset = case["dataset"]
        row = {"dataset": dataset, "label": "shared_hit1", "status": "validating_control", "scores": {}}
        result["arms"].append(row)
        write_json(output / "results.json", result)
        folder = Path(case["output"])
        try:
            control = selected_control(control_report, dataset)
            control_folder = Path(plan["control_root"]) / dataset / "canonical"
            preflight = json.loads((control_folder / "preflight.json").read_text())
            native_control = replay_native(control, {"command": preflight["native"]})
            scores, details, files = readback_scores(repo, dataset, control, native_control, control_folder,
                                                    f"ohmm_hillclimb_r1_canon_{control_nonce}_canonical")
            row["control"] = {"scores": scores, "details": details, "native_summary": native_control["summary"],
                              "checked": [*native_control["checked"], *files, record(control_folder / "preflight.json")],
                              "source_revision": control_report["source_revision"], "initial_min_total_hits": 4,
                              "initial_threshold_provenance": "Frozen control initial-search engine uses its unchanged default 4; not a retroactively recorded metadata field."}
        except Exception as error:
            row.update(status="not_run_unverified_control", error_type=type(error).__name__, error=str(error))
            write_json(output / "results.json", result)
            continue
        folder.mkdir(parents=True, exist_ok=False)
        try:
            for pin in (*checked, *row["control"]["checked"]):
                check(pin)
            verify_input_inventory(case)
            probe = subprocess.check_output([sys.executable, "-c",
                "from orthohmm.search import engine; from orthohmm.search.prefilter import is_c_prefilter_available; "
                "assert not engine.is_cuda_available() and engine.is_c_available() and is_c_prefilter_available()"],
                cwd=source, env=env, text=True)
            row["cpu_backend_probe"] = {"passed": True, "stdout": probe}
            row["status"] = "fresh_native_running"
            write_json(output / "results.json", result)
            row["execution"] = timed_stage(case["command"], folder, "native", source, env)
            if row["execution"]["returncode"]:
                raise RuntimeError("Fresh shared-hit pipeline failed")
            native = validate_native(case, 16)
            metadata = native["metrics"]["metadata"]
            if (metadata.get("search_result_order") != "canonical" or metadata.get("search_min_total_hits") != 1
                    or native["genes"] != native_control["genes"] or native["species"] != native_control["species"]):
                raise ValueError("Native mechanism or complete universe differs from control")
            row["native"] = native
            row["status"] = "native_verified_scoring_pending"
            write_json(output / "results.json", result)
            prediction = folder / "output/orthohmm_working_res/orthohmm_edges_clustered.txt"
            phylo = folder / "output/orthohmm_phylogeny"
            if dataset == "OrthoBench":
                row.update(score_orthobench(repo, prediction))
            elif dataset == "ThreeKingdoms":
                row.update(score_three_kingdoms(repo, prediction, folder))
            else:
                row.update(score_qfo(repo, phylo, folder, f"sharedhit_{nonce}"))
            readback = {**native, "prediction": str(prediction), "phylogeny": str(phylo)}
            actual, details, files = readback_scores(repo, dataset, row, readback, folder,
                                                    f"ohmm_hillclimb_r1_sharedhit_{nonce}")
            row.update(details=details, scoring_readback=files, metric_vector=endpoint_vector(dataset, actual),
                       deltas={name: actual[name] - row["control"]["scores"][name] for name in actual})
            for pin in (*checked, *row["control"]["checked"], *native["products"].values(), *files):
                check(pin)
            if runtime() != plan["runtime"]:
                raise ValueError("Runtime changed during evaluation")
            row["status"] = "independent_content_and_score_readback_verified"
        except Exception as error:
            row.update(status="failed", error_type=type(error).__name__, error=str(error))
        write_json(output / "results.json", result)
        print(dataset, row["status"], flush=True)
    result["status"] = ("finished_all_datasets_verified" if all(
        row["status"] == "independent_content_and_score_readback_verified" for row in result["arms"])
        else "finished_with_explicit_gaps")
    result["host_after"] = host_activity()
    write_json(output / "results.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    preparation = sub.add_parser("prepare")
    for name in ("repo", "source", "output"):
        preparation.add_argument("--" + name, type=Path, required=True)
    worker = sub.add_parser("run")
    worker.add_argument("--execution-plan", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.repo.resolve(), args.source.resolve(), args.output.resolve()) if args.command == "prepare" else run(args.execution_plan.resolve())
    print(result["status"])
    raise SystemExit(0 if result["status"] in {"prepared_waiting_for_control_task", "finished_all_datasets_verified"} else 1)
