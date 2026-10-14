"""Matched native OrthoBench pilot for opt-in singleton tie deferral."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_orthobench_factorial import prepare_partition
from benchmark_tools.replay_phylogeny import load_membership_constraints
from benchmark_tools.report_hillclimb import replay_native, readback_scores
from benchmark_tools.run_canonical_downstream import commands
from benchmark_tools.run_full_search_hillclimb import runtime
from benchmark_tools.run_membership_hillclimb import PRODUCTS, compare_golden, score_orthobench, validate_partition, write_json
from benchmark_tools.run_search_recovery import host_activity
from benchmark_tools.run_three_kingdoms_membership_supplement import sources, timed_stage
from orthohmm.accuracy import load_accuracy_checkpoint
from orthohmm.orthohmm import _expand_phylogeny_candidates

PROTOCOL = "hillclimb_singleton_tie_orthobench_20261014.json"
ARMS = ("first_hit", "defer_ambiguous")


def selected_control(report):
    rows = [r for r in report["arms"] if (r["dataset"], r["label"]) == ("OrthoBench", "canonical")]
    if len(rows) != 1 or rows[0]["status"] != "independent_content_readback_verified":
        raise ValueError("Require the unique verified canonical OrthoBench control")
    return rows[0]


def prepare(repo, source, output):
    root = repo / "benchmarks/results/canonical_downstream_r1_v1_20261014"
    parent_path, report_path = root / "execution_plan.json", root / "results.json"
    parent, report = json.loads(parent_path.read_text()), json.loads(report_path.read_text())
    control = selected_control(report)
    if output.exists() or source == repo or parent["cpu"] != 16:
        raise ValueError("Require fresh output and a private frozen source")
    case = next(c for c in parent["cases"] if c["dataset"] == "OrthoBench")
    control_folder = root / "OrthoBench/canonical"
    preflight = json.loads((control_folder / "preflight.json").read_text())
    native = replay_native(control, {"command": preflight["native"]})
    case = {**case, "baseline_output": str(control_folder / "native")}
    inputs = [record(parent_path), record(report_path), record(source / "benchmark_tools" / PROTOCOL),
              *parent["sources"], *parent["runtime"]["files"], *native["checked"],
              *case["fastas"], *case["scorers"], *control["outputs"],
              control["candidate"]["candidate_partition"], control["candidate"]["membership_constraints"]]
    inputs.extend(record(p) for p in sorted(Path(control["ordering"]["checkpoint"]).iterdir()) if p.is_file())
    phylo = control_folder / "native/orthohmm_phylogeny"
    inputs.extend(record(p) for name in ("alignments", "gene_trees", "checkpoints", "species_tree_inference")
                  for p in sorted((phylo / name).rglob("*")) if p.is_file())
    inputs = list({p["path"]: p for p in inputs}.values())
    source_pins = [*sources(source), record(source / "benchmark_tools/run_singleton_tie_orthobench.slurm")]
    revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    for pin in source_pins:
        relative = Path(pin["path"]).relative_to(source)
        if relative.suffix != ".so":
            blob = subprocess.check_output(["git", "-C", str(source), "show", revision + ":" + str(relative)])
            if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
                raise ValueError("Frozen pilot source differs from committed revision")
    for pin in (*inputs, *source_pins):
        check(pin)
    if runtime() != parent["runtime"]:
        raise ValueError("Scientific runtime differs from control")
    result = {"status": "prepared_singleton_tie_orthobench_pilot", "repo": str(repo), "source_root": str(source),
              "source_revision": revision, "output": str(output), "cpu": 16, "inputs": inputs,
              "sources": source_pins, "runtime": parent["runtime"], "tool_path": parent["tool_path"],
              "case": case, "control": control, "checkpoint": control["ordering"]["checkpoint"],
              "checkpoint_sha": control["ordering"]["manifest"]["sha256"],
              "golden": {name: record(phylo / name) for name in (*PRODUCTS, "species_tree.rooted.nwk")},
              "scope": "Matched cached downstream graph-policy pilot, not fresh search or isolated timing.",
              "defaults_promoted": False, "goal_completed": False}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", result)
    return result


def run(path):
    plan = json.loads(path.read_text())
    source, repo, output = (Path(plan[k]) for k in ("source_root", "repo", "output"))
    if Path(__file__).resolve() != source / "benchmark_tools/run_singleton_tie_orthobench.py" or (output / "results.json").exists():
        raise ValueError("Require frozen runner and preserve all attempts")
    checked = [record(path), *plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]
    for pin in checked:
        check(pin)
    if runtime() != plan["runtime"] or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Runtime or CUDA visibility changed")
    names, species, q, t, scores = load_accuracy_checkpoint(plan["checkpoint"], verify=True)
    if len(names) != 251378 or len(set(species.tolist())) != 12:
        raise ValueError("Wrong full OrthoBench checkpoint universe")
    result = {"status": "running", "execution_plan": record(path), "job_id": os.environ.get("SLURM_JOB_ID"),
              "source_revision": plan["source_revision"], "arms": [], "host_before": host_activity(),
              "scope": plan["scope"], "defaults_promoted": False, "goal_completed": False}
    write_json(output / "results.json", result)
    env = {**os.environ, "PATH": plan["tool_path"], "PYTHONPATH": str(source), "PYTHONHASHSEED": "0",
           "OMP_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "16", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    control_ok = False
    for policy in ARMS:
        row = {"dataset": "OrthoBench", "label": policy, "status": "running", "scores": {}}
        result["arms"].append(row)
        if policy != ARMS[0] and not control_ok:
            row["status"] = "not_run_failed_control_gate"
            write_json(output / "results.json", result)
            continue
        folder = output / policy
        folder.mkdir(exist_ok=False)
        try:
            for pin in checked:
                check(pin)
            if runtime() != plan["runtime"]:
                raise ValueError("Runtime changed before arm execution")
            graph, native_command = commands(source, plan, plan["case"], folder, Path(plan["checkpoint"]), plan["checkpoint_sha"])
            graph.extend(["--singleton-tie-policy", policy])
            write_json(folder / "preflight.json", {"graph": graph, "native": native_command})
            row["graph_execution"] = timed_stage(graph, folder, "graph", source, env)
            if row["graph_execution"]["returncode"]:
                raise RuntimeError("Policy graph replay failed")
            graph_report = json.loads((folder / "graph_metrics.json").read_text())
            if graph_report["parameters"].get("singleton_tie_policy") != policy:
                raise ValueError("Actual graph policy differs from planned policy")
            row["candidate"] = prepare_partition(folder / "graph/orthogroups_profiles_refined.txt", folder / "candidate",
                names, species, (q, t, scores), True, _expand_phylogeny_candidates, load_membership_constraints)
            if policy == ARMS[0]:
                for key in ("candidate_partition", "membership_constraints"):
                    old, new = plan["control"]["candidate"][key], row["candidate"][key]
                    if (old["bytes"], old["sha256"]) != (new["bytes"], new["sha256"]):
                        raise ValueError("Canonical candidate control does not reproduce")
            row["native_execution"] = timed_stage(native_command, folder, "native", source, env)
            if row["native_execution"]["returncode"]:
                raise RuntimeError("Policy native phylogeny replay failed")
            phylo = folder / "native/orthohmm_phylogeny"
            prediction = folder / "native/orthohmm_working_res/orthohmm_edges_clustered.txt"
            row["coverage"] = validate_partition(Path(row["candidate"]["candidate_partition"]["path"]), prediction)
            row["outputs"] = [record(phylo / name) for name in PRODUCTS]
            native = replay_native(row, {"command": native_command})
            if (native["genes"], native["species"]) != (251378, 12):
                raise ValueError("Native universe differs")
            row.update(score_orthobench(repo, prediction))
            actual, details, score_pins = readback_scores(repo, "OrthoBench", row, native, folder, None)
            if policy == ARMS[0]:
                row["golden_matches"] = compare_golden(phylo, plan["golden"])
                if actual != plan["control"]["scores"]:
                    raise ValueError("Canonical control scores differ")
            for pin in (*checked, *native["checked"], *score_pins):
                check(pin)
            if runtime() != plan["runtime"]:
                raise ValueError("Runtime changed during pilot")
            row.update(status="native_and_score_independently_verified", details=details, native_summary=native["summary"],
                       scoring_readback=score_pins)
            if policy == ARMS[0]:
                control_ok = True
        except Exception as error:
            row.update(status="failed", error_type=type(error).__name__, error=str(error))
        write_json(output / "results.json", result)
    result["status"] = ("both_singleton_policies_scored" if all(r["status"] == "native_and_score_independently_verified"
        for r in result["arms"]) else "finished_with_explicit_gaps")
    if result["status"] == "both_singleton_policies_scored":
        delta = result["arms"][1]["scores"]["OrthoBench F1"] - result["arms"][0]["scores"]["OrthoBench F1"]
        result.update(orthobench_f1_delta=delta, passes_development_decline_screen=delta >= -.005)
    result["host_after"] = host_activity()
    write_json(output / "results.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare")
    for name in ("repo", "source", "output"):
        prep.add_argument("--" + name, type=Path, required=True)
    worker = sub.add_parser("run")
    worker.add_argument("--execution-plan", type=Path, required=True)
    args = parser.parse_args()
    result = (prepare(args.repo.resolve(), args.source.resolve(), args.output.resolve()) if args.command == "prepare"
              else run(args.execution_plan.resolve()))
    print(result["status"])
    raise SystemExit(0 if result["status"] in {"prepared_singleton_tie_orthobench_pilot", "both_singleton_policies_scored"} else 1)
