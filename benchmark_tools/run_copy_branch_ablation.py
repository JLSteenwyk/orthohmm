"""Evaluate prespecified copy-split branches using admitted HMM graphs and own trees."""

from __future__ import annotations

import argparse
import ast
from collections import Counter
import csv
import hashlib
import importlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from Bio import SeqIO
from orthohmm.files import fetch_fasta_files
from benchmark_tools.audit_accuracy_checkpoint import audit
from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_orthobench_factorial import prepare_partition
from benchmark_tools.probe_leiden_boundary import saved_fingerprint
from benchmark_tools.replay_high_sensitivity import production_refinement_hits
from benchmark_tools.report_deferred_cells import collect_deferred_cells
from benchmark_tools.report_hillclimb import EXPECTED_GENES, readback_scores, replay_native
from benchmark_tools.run_deferred_copy_split import regenerate_refinement, refinement_changes
from benchmark_tools.run_deferred_inferred_supplement import load_native_edges
from benchmark_tools.run_full_search_hillclimb import runtime as base_runtime
from benchmark_tools.run_held_tree_hillclimb import source_records
from benchmark_tools.run_membership_hillclimb import PRODUCTS, argument, score_qfo, validate_partition, write_json
from benchmark_tools.run_qfo_recovered_assessment import require_darwin_work_path
from benchmark_tools.run_search_recovery import host_activity
from benchmark_tools.run_three_kingdoms_membership_supplement import candidate_control_matches, timed_stage, validate_ownership

PROTOCOL = "hillclimb_copy_branch_ablation_20261013.json"
ARMS = (("large_only500", "phylogeny_deferred_large"), ("medium_only500", "phylogeny_deferred_medium"))
PARENT_LABELS = ("control", "deferred500")


def runtime():
    result = base_runtime()
    dendropy = importlib.import_module("dendropy")
    result["versions"]["dendropy"] = dendropy.__version__
    result["files"].extend(record(p) for p in sorted(Path(dendropy.__file__).parent.rglob("*.py")))
    return result


def terminal_job(accounting, job):
    matches = [r for r in csv.DictReader(io.StringIO(accounting), delimiter="|") if r.get("JobIDRaw") == str(job)]
    if len(matches) != 1:
        raise ValueError("Parent job is not uniquely present in accounting")
    row = matches[0]
    if ((row.get("State"), row.get("ExitCode")) not in {("COMPLETED", "0:0"), ("FAILED", "1:0")}
            or row.get("NodeList") != "bizon" or row.get("AllocCPUS") != "16"):
        raise ValueError("Require a terminal 16-CPU local parent with its explicit exit status")
    return row


def definitions(path):
    return {node.name: ast.dump(node, include_attributes=False) for node in ast.parse(path.read_text()).body
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))}


def refinement_constants(path):
    result = {}
    for node in ast.parse(path.read_text()).body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and (target.id.startswith("DEFAULT_")
                        or target.id == "PHYLOGENY_DEFERRED_COPY_SPLIT_MAX_SIZE"):
                    result[target.id] = ast.dump(node.value, include_attributes=False)
    return result


def replay_equivalence(previous, source):
    excluded = {"refinement.py", "args_processing.py", "orthohmm.py", "parser.py"}
    old_files = {p.relative_to(previous) for p in (previous / "orthohmm").rglob("*.py") if p.name not in excluded}
    new_files = {p.relative_to(source) for p in (source / "orthohmm").rglob("*.py") if p.name not in excluded}
    if old_files != new_files:
        raise ValueError("Downstream native source inventory changed")
    paths = sorted(old_files | {Path("benchmark_tools/replay_phylogeny.py")})
    checked = []
    for relative in paths:
        old, new = record(previous / relative), record(source / relative)
        if old["sha256"] != new["sha256"]:
            raise ValueError("Downstream native source changed: " + str(relative))
        checked.extend((old, new))
    for relative, allowed in (("orthohmm/orthohmm.py", None),
                              ("orthohmm/refinement.py", {"refine_cluster_indices", "_split_high_duplication_index_clusters"})):
        a, b = definitions(previous / relative), definitions(source / relative)
        names = {"_expand_phylogeny_candidates"} if allowed is None else set(a) - allowed
        if (allowed is not None and set(a) != set(b)) or any(a.get(name) != b.get(name) for name in names):
            raise ValueError("Candidate expansion or unchanged refinement helpers differ")
        checked.extend((record(previous / relative), record(source / relative)))
    if refinement_constants(previous / "orthohmm/refinement.py") != refinement_constants(source / "orthohmm/refinement.py"):
        raise ValueError("Refinement defaults or original deferral ceiling changed")
    return {"status": "downstream_native_bytes_and_candidate_helpers_unchanged", "checked": checked,
            "scope": "Does not assert equivalence of the intentionally changed copy-split branches; exact controls remain required."}


def native_command(original, source, folder, candidate, constraints, checkpoint):
    expected = {"--species-tree-mode": "infer", "--root-rule": "species_overlap",
        "--pair-rule": "positive_paralogy", "--species-tree-rooting": "min_variance",
        "--satellite-membership-rule": "high_confidence_pair", "--cpu": "16"}
    if any(argument(original, flag) != value for flag, value in expected.items()) or any(
            flag in original for flag in ("--species-tree", "--unconstrained-membership", "--official-benchmark")):
        raise ValueError("Require unchanged own-tree, duplication-aware, constrained parent inference")
    command = list(original)
    command[:2] = [sys.executable, str(source / "benchmark_tools/replay_phylogeny.py")]
    for flag, value in (("--candidate-clusters", candidate), ("--membership-constraints", constraints),
                        ("--output-directory", folder / "output"), ("--json", folder / "metrics.json"),
                        ("--checkpoint-source", checkpoint)):
        argument(command, flag)
        command[command.index(flag) + 1] = str(value)
    return command


def cases_for(parent, dataset):
    cases = parent["cases"] if dataset == "QfO" else next(c for c in parent["cases"] if c["dataset"] == dataset)["arms"]
    if tuple(c["label"] for c in cases) != PARENT_LABELS:
        raise ValueError("Require both ordered parent arms")
    return cases


def prepare(repo, source, output, dataset, parent_root, parent_job, fixed_root=None):
    if output.exists() or dataset not in EXPECTED_GENES or not str(parent_job).isdigit() or int(parent_job) <= 0:
        raise ValueError("Require fresh namespace, known dataset and explicit positive parent job")
    subprocess.run(["git", "-C", str(source), "diff", "--exit-code", "HEAD", "--", "orthohmm", "benchmark_tools"], check=True)
    parent_path = parent_root / "execution_plan.json"
    parent = json.loads(parent_path.read_text())
    if Path(parent["repo"]).resolve() != repo or parent["cpu"] != 16:
        raise ValueError("Wrong parent repository or resources")
    previous = Path(parent["source_root"])
    equivalence = replay_equivalence(previous, source)
    old_cases = cases_for(parent, dataset)
    inputs = [record(parent_path), *parent["inputs"], *parent["sources"], *parent["runtime"]["files"],
              *equivalence["checked"], record(source / "benchmark_tools" / PROTOCOL)]
    if dataset == "QfO":
        if fixed_root is None or parent["parent_plan"] != record(fixed_root / "execution_plan.json"):
            raise ValueError("Require the matching completed fixed-parent graph input")
        fixed = json.loads((fixed_root / "execution_plan.json").read_text())
        data = {k: fixed[k] for k in ("checkpoint", "checkpoint_manifest", "fasta_directory", "graph", "partition", "golden_refinement")}
        inputs.extend([record(fixed_root / "execution_plan.json"), *fixed["inputs"], *fixed["sources"]])
        checkpoint_source = Path(fixed["cases"][1]["output"]) / "output"
        parent_report = parent_root / "results.json"
    else:
        if fixed_root is not None:
            raise ValueError("Narrow-panel inputs do not take a fixed QfO parent")
        original = next(c for c in parent["cases"] if c["dataset"] == dataset)
        data = {k: original[k] for k in ("checkpoint", "checkpoint_manifest", "fasta_directory", "native_graph", "native_edges")}
        child = json.loads((parent_root / dataset / "results.json").read_text())
        if child.get("dataset") != dataset or child.get("status") != "finished_all_arms_verified":
            raise ValueError("Require a completed independently verified narrow-panel child")
        data["partition"] = child["graph_partition"]
        inputs.extend([record(parent_root / dataset / "results.json"), data["partition"]])
        checkpoint_source = Path(old_cases[1]["output"]) / "output"
        parent_report = parent_root / dataset / "results.json"
    candidates = []
    for case in old_cases:
        candidate = record(Path(argument(case["command"], "--candidate-clusters")))
        constraints = record(Path(argument(case["command"], "--membership-constraints")))
        inputs.extend((candidate, constraints))
        candidates.append({"candidate": candidate, "constraints": constraints})
    for directory in ("alignments", "checkpoints", "gene_trees"):
        paths = [p for p in sorted((checkpoint_source / "orthohmm_phylogeny" / directory).rglob("*")) if p.is_file()]
        if not paths:
            raise ValueError("Missing matching immutable raw-tree cache")
        inputs.extend(record(p) for p in paths)
    arms = []
    nonce = hashlib.sha256(str(output).encode()).hexdigest()[:16]
    for index, (label, profile) in enumerate(ARMS):
        folder = output / label
        candidate = folder / "candidate/orthohmm_working_res/orthohmm_edges_clustered.txt"
        constraints = candidate.with_name("phylogeny_candidate_merges.json")
        scoring_id = f"cb_{nonce}_{index}"
        if dataset == "QfO":
            require_darwin_work_path(repo / "qfo_benchmark/w" / ("hc1_" + scoring_id))
        arms.append({"label": label, "profile": profile, "output": str(folder), "scoring_id": scoring_id,
            "command": native_command(old_cases[0]["command"], source, folder, candidate, constraints, checkpoint_source)})
    current_runtime = runtime()
    if (parent["runtime"]["python"] != current_runtime["python"] or any(
            current_runtime["versions"].get(k) != v for k, v in parent["runtime"]["versions"].items())):
        raise ValueError("Scientific runtime versions changed")
    inputs = list({item["path"]: item for item in inputs}.values())
    for item in inputs:
        check(item)
    result = {"status": "prepared_waiting_for_terminal_parent", "repo": str(repo), "source_root": str(source),
        "source_revision": subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip(),
        "output": str(output), "dataset": dataset, "parent_root": str(parent_root), "parent_job": str(parent_job),
        "parent_plan": record(parent_path), "parent_results": str(parent_report), "parent_cases": old_cases,
        "data": data, "parent_candidates": candidates, "replay_equivalence": equivalence, "cpu": 16,
        "inputs": inputs, "sources": [*source_records(source), record(source / "benchmark_tools/run_copy_branch_ablation.slurm")],
        "runtime": current_runtime, "tool_path": parent["tool_path"], "arms": arms,
        "scope": "Development-exposed cached HMM graph branch ablation with own inferred trees; exact structural reuse for inactive narrow-panel branches. Not fresh full-pipeline runtime."}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", result)
    return result


def verify(plan, identities=()):
    for item in [*plan["inputs"], *plan["sources"], *plan["runtime"]["files"], *identities]:
        check(item)
    if runtime() != plan["runtime"] or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Runtime or CUDA visibility changed")


def resolve_parent(plan, output):
    accounting = subprocess.check_output(["sacct", "-j", plan["parent_job"], "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = terminal_job(accounting, plan["parent_job"])
    top_path = Path(plan["parent_root"]) / "results.json"
    top = json.loads(top_path.read_text())
    parent = json.loads(Path(plan["parent_plan"]["path"]).read_text())
    if (top.get("job_id") != plan["parent_job"] or top.get("execution_plan") != plan["parent_plan"]
            or top.get("source_revision") != parent["source_revision"]
            or not top.get("status", "").startswith(("finished", "failed"))):
        raise ValueError("Parent producer does not confirm the exact terminal attempt")
    child = json.loads(Path(plan["parent_results"]).read_text())
    if tuple(r.get("label") for r in child.get("arms", [])) != PARENT_LABELS:
        raise ValueError("Unexpected terminal parent arm inventory")
    if child["arms"][0].get("status") != "independent_content_readback_verified":
        raise ValueError("No successful terminal matched control")
    if plan["dataset"] != "QfO" and child.get("status") != "finished_all_arms_verified":
        raise ValueError("Narrow-panel structural reuse requires both verified parent arms")
    if plan["dataset"] == "QfO" and child.get("status") not in ("finished_all_arms_verified", "finished_with_explicit_gaps"):
        raise ValueError("QfO parent has no explicit terminal outcome")
    if plan["dataset"] == "QfO":
        if child["status"] == "finished_all_arms_verified" and scheduler["State"] != "COMPLETED":
            raise ValueError("QfO parent success disagrees with scheduler state")
        if child["status"] == "finished_with_explicit_gaps" and (scheduler["State"] != "FAILED"
                or child["arms"][1].get("status") != "failed" or not child["arms"][1].get("error")):
            raise ValueError("Require the specific terminal QfO variant failure, not an unexplained gap")
    snaps = output / "parent_readback_snapshots"
    snaps.mkdir()
    cells, captures = collect_deferred_cells(Path(plan["repo"]), "deferred_inferred", plan["dataset"],
                                            Path(plan["parent_root"]), snaps)
    if cells[0]["status"] != "independent_content_readback_verified" or (plan["dataset"] != "QfO" and
            cells[1]["status"] != "independent_content_readback_verified"):
        raise ValueError("Independent terminal parent readback failed")
    checked = [record(top_path), record(Path(plan["parent_results"])), plan["parent_plan"],
        *(c["snapshot"] for c in captures), *(r for cell in cells if cell["status"] == "independent_content_readback_verified"
                                                for r in cell["checked"])]
    for item in checked:
        check(item)
    resolution = {"status": "terminal_matched_parent_independently_admitted", "scheduler": scheduler,
        "cells": cells, "captures": captures, "checked": checked, "parent_rows": child["arms"]}
    path = output / "resolved_parent_inputs.json"
    write_json(path, resolution)
    return resolution, [record(path), *checked]


def load_graph(plan):
    from orthohmm.accuracy import load_accuracy_checkpoint, read_index_clusters
    data = plan["data"]
    checkpoint = Path(data["checkpoint"])
    audited = audit(checkpoint, data["checkpoint_manifest"]["sha256"])
    names, species, queries, targets, scores = load_accuracy_checkpoint(checkpoint, verify=False)
    owners = {}
    fasta = Path(data["fasta_directory"])
    for filename in sorted(fetch_fasta_files(str(fasta))):
        for sequence in SeqIO.parse(fasta / filename, "fasta"):
            if sequence.id in owners or not sequence.seq:
                raise ValueError("Duplicate or empty FASTA gene")
            owners[sequence.id] = filename
    ownership = validate_ownership(names, species, owners)
    if ownership["genes"] != EXPECTED_GENES[plan["dataset"]] or ownership["species"] != (78 if plan["dataset"] == "QfO" else 12):
        raise ValueError("Wrong full dataset gene/species universe")
    read_partition(Path(data["partition"]["path"]), set(names))
    clusters = read_index_clusters(data["partition"]["path"], names)
    if "graph" in data:
        payload = Path(data["graph"]["payload"])
        if (payload / "gene_names.txt").read_text().splitlines() != names or saved_fingerprint(payload) != data["graph"]["fingerprint"]:
            raise ValueError("Native graph fingerprint or checkpoint gene order changed")
        edges = tuple(np.load(payload / (name + ".npy"), mmap_mode="r", allow_pickle=False)
                      for name in ("sources", "targets", "weights"))
        hits = None
    else:
        edge = load_native_edges(Path(data["native_graph"]["path"]), names, data["native_edges"])
        edges = edge.sources, edge.targets, edge.weights
        hits = production_refinement_hits(queries, targets, scores, species)
    return names, species, (queries, targets, scores), clusters, edges, hits, audited, ownership


def construct(plan, resolution, output):
    from orthohmm.orthohmm import _expand_phylogeny_candidates
    from orthohmm.refinement import resolve_refinement_profile
    from benchmark_tools.replay_phylogeny import load_membership_constraints
    names, species, hits, clusters, edges, directed, audited, ownership = load_graph(plan)
    profiles = (("control", "default"), ("deferred500", "phylogeny_deferred"), *ARMS)
    rows, refined = {}, {}
    for label, profile in profiles:
        expected = ({} if label == "control" else {"defer_copy_split_max_size": 500}
                    if label == "deferred500" else {"defer_copy_split_max_size": 500,
                                                   "defer_copy_split_scope": "large" if label == "large_only500" else "medium"})
        if resolve_refinement_profile(profile) != expected:
            raise ValueError("Scoped profile differs from the prespecified branch/ceiling")
        started = time.monotonic()
        folder = output / label
        folder.mkdir(exist_ok=False)
        seed, trace = folder / "refined_seed.txt", []
        groups = regenerate_refinement(clusters, edges, species, names, profile, seed,
                                       hits=directed, copy_split_trace=trace)
        refined[label] = groups
        trace_path = folder / "copy_split_trace.json"
        write_json(trace_path, {"scope": "Executed post-merge copy-split eligibility and bounded deferral; not reference labels.", "events": trace})
        counts = Counter()
        for event in trace:
            if event["deferred"] and (event["size"] > 500 or label == "control"
                    or label == "large_only500" and not event["large_copy_eligible"]
                    or label == "medium_only500" and not event["medium_copy_eligible"]):
                raise ValueError("Executed deferral exceeded its branch/size safeguards")
            for branch in ("large", "medium"):
                if event[branch + "_copy_eligible"]:
                    counts[branch + "_eligible_families"] += 1
                    if event["deferred"]:
                        counts[branch + "_deferred_families"] += 1
                        counts[branch + "_deferred_genes"] += event["size"]
        candidate = prepare_partition(seed, folder / "candidate", names, species, hits, True,
                                      _expand_phylogeny_candidates, load_membership_constraints)
        rows[label] = {"label": label, "profile": profile, "status": "candidates_verified_unscored", "scores": {},
            "seed": record(seed), "candidate": candidate, "copy_split_trace": record(trace_path),
            "copy_split_counts": dict(counts), "candidate_wall_s": time.monotonic() - started,
            "copy_split_count_scope": "Executed eligibility/deferral predicates after merging; net partition effects are recorded separately."}
        if label in PARENT_LABELS:
            index = PARENT_LABELS.index(label)
            candidate_control_matches(Path(candidate["candidate_partition"]["path"]), Path(candidate["membership_constraints"]["path"]),
                                      plan["parent_candidates"][index]["candidate"], plan["parent_candidates"][index]["constraints"])
        if candidate["expansion"]["parameters"] != rows["control"]["candidate"]["expansion"]["parameters"]:
            raise ValueError("Satellite expansion parameters changed")
    if "golden_refinement" in plan["data"] and rows["control"]["seed"]["sha256"] != plan["data"]["golden_refinement"]["sha256"]:
        raise ValueError("Default refinement does not reproduce admitted bytes")
    if plan["dataset"] != "QfO":
        for label, reference in (("large_only500", "control"), ("medium_only500", "deferred500")):
            candidate = rows[label]["candidate"]
            candidate_control_matches(Path(candidate["candidate_partition"]["path"]), Path(candidate["membership_constraints"]["path"]),
                rows[reference]["candidate"]["candidate_partition"], rows[reference]["candidate"]["membership_constraints"])
    changes = {label: refinement_changes(clusters, refined["control"], refined[label], species, 500) for label, _ in ARMS}
    return rows, {"checkpoint_audit": audited, "ownership": ownership, "refinement_changes": changes}


def require_markers(summary):
    count = summary.get("species_tree_families")
    if type(count) is not int or count <= 0:
        raise ValueError("Own inferred species tree lacks valid native markers")


def run(path):
    plan = json.loads(path.read_text())
    source, output, repo = (Path(plan[k]) for k in ("source_root", "output", "repo"))
    if Path(__file__).resolve() != source / "benchmark_tools/run_copy_branch_ablation.py":
        raise ValueError("Execute frozen branch-ablation source")
    if (output / "results.json").exists() or tuple((a["label"], a["profile"]) for a in plan["arms"]) != ARMS:
        raise ValueError("Existing attempt or unexpected prespecified arms")
    identity = record(path)
    report = {"status": "running", "dataset": plan["dataset"], "job_id": os.environ.get("SLURM_JOB_ID"),
        "execution_plan": identity, "source_revision": plan["source_revision"], "scope": plan["scope"], "arms": [],
        "defaults_promoted": False, "goal_completed": False, "host_before": host_activity()}
    write_json(output / "results.json", report)
    checked = [identity]
    env = {**os.environ, "PATH": plan["tool_path"], "PYTHONPATH": str(source), "CUDA_VISIBLE_DEVICES": "",
           "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
           "NUMBA_NUM_THREADS": "16"}
    try:
        verify(plan, checked)
        resolution, pins = resolve_parent(plan, output)
        checked.extend(pins)
        report["parent_admission"] = pins[0]
        rows, diagnostics = construct(plan, resolution, output)
        report.update(diagnostics)
        report["structural_controls"] = {label: rows[label] for label in PARENT_LABELS}
        checked.extend(item for row in rows.values() for item in (row["seed"], row["copy_split_trace"],
            row["candidate"]["candidate_partition"], row["candidate"]["membership_constraints"]))
        for arm in plan["arms"]:
            row = rows[arm["label"]]
            report["arms"].append(row)
            verify(plan, checked)
            if plan["dataset"] != "QfO":
                index = 0 if arm["label"] == "large_only500" else 1
                old = resolution["parent_rows"][index]
                native_case = plan["parent_cases"][index]
                native = replay_native(old, native_case)
                score_folder = Path(native_case["output"])
                scoring_row = old
                row["execution_kind"] = "exact_candidate_structural_reuse_no_native_rerun"
                row["reused_parent_label"] = PARENT_LABELS[index]
            else:
                row["status"] = "native_running"
                write_json(output / "results.json", report)
                folder = Path(arm["output"])
                row["native_execution"] = timed_stage(arm["command"], folder, "native", source, env)
                if row["native_execution"]["returncode"]:
                    raise RuntimeError("Branch-specific native inferred phylogeny failed")
                phylo = folder / "output/orthohmm_phylogeny"
                prediction = folder / "output/orthohmm_working_res/orthohmm_edges_clustered.txt"
                row["coverage"] = validate_partition(Path(row["candidate"]["candidate_partition"]["path"]), prediction)
                row["outputs"] = [record(phylo / name) for name in PRODUCTS]
                native_case = arm
                native = replay_native(row, arm)
                require_markers(native["summary"])
                row["status"] = "native_verified_scoring_pending"
                write_json(output / "results.json", report)
                verify(plan, [*checked, *native["checked"]])
                row.update(score_qfo(repo, phylo, folder, arm["scoring_id"]))
                score_folder, scoring_row = folder, row
                row["execution_kind"] = "new_cached_graph_branch_inference_own_tree"
            scores, details, files = readback_scores(repo, plan["dataset"], scoring_row, native, score_folder,
                "ohmm_hillclimb_r1_" + arm["scoring_id"] if plan["dataset"] == "QfO" else "unused")
            if native["genes"] != EXPECTED_GENES[plan["dataset"]] or native["species"] != (78 if plan["dataset"] == "QfO" else 12):
                raise ValueError("Branch inference lost the full benchmark universe")
            require_markers(native["summary"])
            checked.extend([*native["checked"], *files])
            verify(plan, checked)
            row.update(scores=scores, details=details, native_case=native_case, native_products=native["products"],
                native_summary=native["summary"], native_readback=files, scoring_output=str(score_folder),
                native_coverage={k: native[k] for k in ("genes", "species", "groups", "pairs")},
                status="independent_content_readback_verified", inputs_reverified_unchanged=True)
            write_json(output / "results.json", report)
        report["status"] = "finished_all_branch_arms_verified"
    except Exception as error:
        report.update(status="failed_with_explicit_gaps", error_type=type(error).__name__, error=str(error))
        if report["arms"] and report["arms"][-1]["status"] != "independent_content_readback_verified":
            report["arms"][-1].update(status="failed", error_type=type(error).__name__, error=str(error))
    report["host_after"] = host_activity()
    write_json(output / "results.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    prep = modes.add_parser("prepare")
    for name in ("repo", "source-root", "output", "parent-root"):
        prep.add_argument("--" + name, type=Path, required=True)
    prep.add_argument("--dataset", choices=tuple(EXPECTED_GENES), required=True)
    prep.add_argument("--parent-job", required=True)
    prep.add_argument("--fixed-root", type=Path)
    worker = modes.add_parser("run")
    worker.add_argument("--execution-plan", type=Path, required=True)
    args = parser.parse_args()
    result = (prepare(args.repo.resolve(), args.source_root.resolve(), args.output.resolve(), args.dataset,
        args.parent_root.resolve(), args.parent_job, args.fixed_root.resolve() if args.fixed_root else None)
        if args.mode == "prepare" else run(args.execution_plan.resolve()))
    print(json.dumps({"status": result["status"]}))
    return 0 if result["status"] in ("prepared_waiting_for_terminal_parent", "finished_all_branch_arms_verified") else 1


if __name__ == "__main__":
    raise SystemExit(main())
