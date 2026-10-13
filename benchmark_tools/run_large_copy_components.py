"""Evaluate bounded HMM-supported components after completed copy-branch admission."""

import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from orthohmm.refinement import resolve_refinement_profile
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_orthobench_factorial import prepare_partition
from benchmark_tools.report_hillclimb import ENDPOINTS, EXPECTED_GENES, endpoint_vector, readback_scores, replay_native
from benchmark_tools.run_copy_branch_ablation import (
    ARMS, load_graph, native_command, replay_equivalence, require_markers, runtime, terminal_job,
)
from benchmark_tools.run_deferred_copy_split import regenerate_refinement
from benchmark_tools.run_held_tree_hillclimb import source_records
from benchmark_tools.run_membership_hillclimb import PRODUCTS, argument, score_qfo, validate_partition, write_json
from benchmark_tools.run_qfo_recovered_assessment import require_darwin_work_path
from benchmark_tools.run_search_recovery import host_activity
from benchmark_tools.run_three_kingdoms_membership_supplement import candidate_control_matches, timed_stage

PROTOCOL = "hillclimb_large_copy_components_20261013.json"
PROFILE = "phylogeny_components"
BRANCH_LABELS = ("control", "deferred500", "large_only500", "medium_only500")


def admitted_branch_cells(comparison, dataset):
    result = {}
    for label in BRANCH_LABELS:
        methods = [m for m in comparison["methods"] if m["key"] == "copy_branch:" + label]
        if len(methods) != 1 or set(methods[0]["scores"]) != set(ENDPOINTS):
            raise ValueError("Require complete independently admitted branch vectors before component inference")
        vector = {}
        for panel in EXPECTED_GENES:
            cells = [c for c in methods[0]["cells"] if c["dataset"] == panel]
            if len(cells) != 1 or cells[0]["status"] != "independent_content_readback_verified":
                raise ValueError("Require one independently admitted branch cell per dataset")
            endpoint_vector(panel, {k.split(": ", 1)[1]: v for k, v in cells[0]["scores"].items()})
            vector.update(cells[0]["scores"])
            if panel == dataset:
                result[label] = cells[0]
        if vector != methods[0]["scores"]:
            raise ValueError("Branch vector differs from its independently admitted cells")
    return result


def executed_component_counts(control, variant, dataset):
    if len(control) != len(variant):
        raise ValueError("Control and component trace inventories differ")
    counts = Counter()
    for baseline, event in zip(control, variant):
        if any(event.get(k) != v for k, v in baseline.items()) or event["deferred"]:
            raise ValueError("Component retention changed original eligibility or enabled whole-family deferral")
        if event["large_copy_eligible"]:
            if dataset != "QfO":
                raise ValueError("Large-copy predicate cannot activate on a narrow panel")
            sizes = event["component_sizes"]
            if (not sizes or any(type(s) is not int or s <= 0 for s in sizes)
                    or sum(sizes) != event["size"] or event["component_min_edge_weight"] != 1.5
                    or event["retain_large_copy_components_max_size"] != 500
                    or event["retained_component_sizes"] != [s for s in sizes if 1 < s <= 500]
                    or event["oversized_component_sizes"] != [s for s in sizes if s > 500]
                    or event["isolated_genes"] != sizes.count(1)):
                raise ValueError("Executed components violate the fixed score/size/fallback policy")
            counts.update(families=1, genes=event["size"], components=len(sizes),
                bounded_components=sum(1 < s <= 500 for s in sizes),
                bounded_component_genes=sum(s for s in sizes if 1 < s <= 500),
                oversized_components=sum(s > 500 for s in sizes),
                oversized_component_genes=sum(s for s in sizes if s > 500), isolated_genes=sizes.count(1))
        elif event != baseline:
            raise ValueError("Non-large refinement trace changed")
    return dict(counts)


def prepare(repo, source, output, dataset, branch_root, comparison_path):
    if output.exists() or dataset not in EXPECTED_GENES:
        raise ValueError("Require fresh output and a known dataset")
    if resolve_refinement_profile(PROFILE) != {"retain_large_copy_components_max_size": 500}:
        raise ValueError("Component profile differs from the prespecified bound")
    if Path(__file__).resolve() != source / "benchmark_tools/run_large_copy_components.py":
        raise ValueError("Prepare using the explicitly bound frozen source")
    subprocess.run(["git", "-C", str(source), "diff", "--exit-code", "HEAD", "--", "orthohmm", "benchmark_tools"], check=True)
    comparison = json.loads(comparison_path.read_text())
    cells = admitted_branch_cells(comparison, dataset)
    parent_path = branch_root / dataset / "execution_plan.json"
    producer_path = parent_path.with_name("results.json")
    parent = json.loads(parent_path.read_text())
    producer = json.loads(producer_path.read_text())
    accounting = subprocess.check_output(["sacct", "-j", producer["job_id"], "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = terminal_job(accounting, producer["job_id"])
    if (scheduler["State"] != "COMPLETED" or producer["status"] != "finished_all_branch_arms_verified"
            or producer["execution_plan"] != record(parent_path) or producer["source_revision"] != parent["source_revision"]
            or Path(parent["repo"]).resolve() != repo or parent["dataset"] != dataset
            or tuple((a["label"], a["profile"]) for a in parent["arms"]) != ARMS
            or tuple(a["label"] for a in producer["arms"]) != tuple(a[0] for a in ARMS)
            or any(row["status"] != "independent_content_readback_verified" for row in producer["arms"])):
        raise ValueError("Require the exact successful terminal branch producer")
    for cell in cells.values():
        check(cell["source_snapshot"])
        captured = json.loads(Path(cell["source_snapshot"]["path"]).read_text())
        if captured["execution_plan"] != producer["execution_plan"] or captured["job_id"] != producer["job_id"]:
            raise ValueError("Comparison cell belongs to another branch execution")
    equivalence = replay_equivalence(Path(parent["source_root"]), source)
    check(producer["parent_admission"])
    resolution = json.loads(Path(producer["parent_admission"]["path"]).read_text())
    if resolution["status"] != "terminal_matched_parent_independently_admitted" or resolution["parent_rows"][0]["label"] != "control":
        raise ValueError("Missing exact matched native control resolution")
    current_runtime = runtime()
    if parent["runtime"] != current_runtime:
        raise ValueError("Scientific runtime changed since matched branch inference")
    probe_path = repo / "benchmark_tools/results/large_copy_component_feasibility_v3_20261013.json"
    probe_readback = repo / "benchmark_tools/results/large_copy_component_independent_readback_20261013.json"
    probe = json.loads(probe_path.read_text())
    readback = json.loads(probe_readback.read_text())
    if (probe["status"] != "complete_all_large_copy_component_feasibility_verified"
            or readback["report"] != record(probe_path) or readback["counts"] != probe["counts"]):
        raise ValueError("Missing independently verified all-family component feasibility")
    inputs = [record(comparison_path), record(parent_path), record(producer_path),
        record(source / "benchmark_tools" / PROTOCOL), *parent["inputs"], *parent["sources"],
        *equivalence["checked"], producer["parent_admission"], record(probe_path), record(probe_readback),
        *[c["source_snapshot"] for c in cells.values()]]
    inputs = list({r["path"]: r for r in inputs}.values())
    for pin in inputs:
        check(pin)
    folder = output / "components500"
    candidate = folder / "candidate/orthohmm_working_res/orthohmm_edges_clustered.txt"
    constraints = candidate.with_name("phylogeny_candidate_merges.json")
    checkpoint = Path(argument(parent["arms"][0]["command"], "--checkpoint-source"))
    scoring_id = "lc_" + hashlib.sha256(str(output).encode()).hexdigest()[:16]
    if dataset == "QfO":
        require_darwin_work_path(repo / "qfo_benchmark/w" / ("hc1_" + scoring_id))
    plan = {"status": "prepared_after_complete_branch_admission", "repo": str(repo), "source_root": str(source),
        "source_revision": subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip(),
        "output": str(output), "dataset": dataset, "data": parent["data"], "cpu": 16,
        "branch_plan": record(parent_path), "branch_producer": record(producer_path), "branch_scheduler": scheduler,
        "comparison": record(comparison_path), "matched_control": cells["control"],
        "structural_control": producer["structural_controls"]["control"],
        "parent_native_case": parent["parent_cases"][0], "parent_native_row": resolution["parent_rows"][0],
        "checkpoint_source": str(checkpoint), "tool_path": parent["tool_path"],
        "expected_component_counts": {k: v for k, v in probe["counts"].items() if not k.startswith("over500_family_")},
        "replay_equivalence": equivalence, "inputs": inputs, "runtime": current_runtime,
        "sources": [*source_records(source), record(source / "benchmark_tools/run_large_copy_components.slurm"),
                    *[record(p) for p in sorted((source / "orthohmm/search/csrc").glob("*.c"))],
                    *[record(p) for p in sorted((source / "orthohmm/search/csrc").glob("*.so"))]],
        "case": {"label": "components500", "profile": PROFILE, "output": str(folder), "scoring_id": scoring_id,
                 "command": native_command(parent["parent_cases"][0]["command"], source, folder, candidate, constraints, checkpoint)},
        "scope": "Development-exposed cached HMM graph component ablation with own inferred trees. Exact inactive narrow-panel reuse; not fresh full-pipeline timing."}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", plan)
    return plan


def verify(plan, extra=()):
    source, folder = Path(plan["source_root"]), Path(plan["output"]) / "components500"
    candidate = folder / "candidate/orthohmm_working_res/orthohmm_edges_clustered.txt"
    expected = native_command(plan["parent_native_case"]["command"], source, folder, candidate,
        candidate.with_name("phylogeny_candidate_merges.json"), Path(plan["checkpoint_source"]))
    if (plan["cpu"] != 16 or plan["case"]["profile"] != PROFILE or plan["case"]["label"] != "components500"
            or plan["case"]["command"] != expected or Path(plan["case"]["output"]).resolve() != folder.resolve()
            or resolve_refinement_profile(PROFILE) != {"retain_large_copy_components_max_size": 500}):
        raise ValueError("Component execution changed the prespecified inference protocol")
    for pin in [*plan["inputs"], *plan["sources"], *plan["runtime"]["files"], *extra]:
        check(pin)
    if runtime() != plan["runtime"] or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Runtime or CUDA visibility changed")


def construct(plan, output):
    from orthohmm.orthohmm import _expand_phylogeny_candidates
    from benchmark_tools.replay_phylogeny import load_membership_constraints
    names, species, hits, clusters, edges, directed, audited, ownership = load_graph(plan)
    rows, traces = {}, {}
    for label, profile in (("control", "default"), ("components500", PROFILE)):
        started = time.monotonic()
        folder = output / label
        folder.mkdir(exist_ok=False)
        seed, trace = folder / "refined_seed.txt", []
        regenerate_refinement(clusters, edges, species, names, profile, seed, hits=directed, copy_split_trace=trace)
        traces[label] = trace
        trace_path = folder / "copy_split_trace.json"
        write_json(trace_path, {"scope": "Executed label-blind post-merge predicates and bounded components.", "events": trace})
        candidate = prepare_partition(seed, folder / "candidate", names, species, hits, True,
                                      _expand_phylogeny_candidates, load_membership_constraints)
        rows[label] = {"label": label, "profile": profile, "seed": record(seed), "copy_split_trace": record(trace_path),
                       "candidate": candidate, "candidate_wall_s": time.monotonic() - started}
    original = plan["structural_control"]
    if rows["control"]["seed"]["sha256"] != original["seed"]["sha256"]:
        raise ValueError("Default refinement failed exact seed reproduction")
    for label in (("control", "components500") if plan["dataset"] != "QfO" else ("control",)):
        candidate_control_matches(Path(rows[label]["candidate"]["candidate_partition"]["path"]),
            Path(rows[label]["candidate"]["membership_constraints"]["path"]),
            original["candidate"]["candidate_partition"], original["candidate"]["membership_constraints"])
    if rows["components500"]["candidate"]["expansion"]["parameters"] != rows["control"]["candidate"]["expansion"]["parameters"]:
        raise ValueError("Satellite expansion parameters changed")
    counts = executed_component_counts(traces["control"], traces["components500"], plan["dataset"])
    if plan["dataset"] == "QfO" and not counts.get("bounded_components"):
        raise ValueError("Component mechanism did not activate; refuse ineffective duplicate inference")
    if plan["dataset"] == "QfO" and counts != plan["expected_component_counts"]:
        raise ValueError("Executed all-family decomposition differs from admitted feasibility")
    return rows, {"checkpoint_audit": audited, "ownership": ownership, "component_counts": counts}


def native_inputs_unchanged(rows):
    return all(rows["components500"]["candidate"][key]["sha256"] == rows["control"]["candidate"][key]["sha256"]
               for key in ("candidate_partition", "membership_constraints"))


def run(path):
    plan = json.loads(path.read_text())
    source, repo, output = (Path(plan[k]) for k in ("source_root", "repo", "output"))
    if Path(__file__).resolve() != source / "benchmark_tools/run_large_copy_components.py" or (output / "results.json").exists():
        raise ValueError("Require frozen source and an untouched execution attempt")
    case = plan["case"]
    report = {"status": "running", "execution_plan": record(path), "dataset": plan["dataset"],
        "source_revision": plan["source_revision"], "job_id": os.environ.get("SLURM_JOB_ID"), "scope": plan["scope"],
        "scores": {}, "defaults_promoted": False, "goal_completed": False, "host_before": host_activity()}
    write_json(output / "results.json", report)
    checked = [record(path)]
    try:
        verify(plan, checked)
        rows, diagnostics = construct(plan, output)
        report.update(diagnostics, structural_controls=rows)
        row = rows["components500"]
        checked.extend(pin for r in rows.values() for pin in (r["seed"], r["copy_split_trace"],
            r["candidate"]["candidate_partition"], r["candidate"]["membership_constraints"]))
        if plan["dataset"] == "QfO" and not native_inputs_unchanged(rows):
            report["status"] = "native_running"
            write_json(output / "results.json", report)
            env = {**os.environ, "PATH": plan["tool_path"], "PYTHONPATH": str(source), "CUDA_VISIBLE_DEVICES": "",
                   "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "16"}
            report["native_execution"] = timed_stage(case["command"], Path(case["output"]), "native", source, env)
            if report["native_execution"]["returncode"]:
                raise RuntimeError("Component native inferred phylogeny failed")
            phylo = Path(case["output"]) / "output/orthohmm_phylogeny"
            report["coverage"] = validate_partition(Path(row["candidate"]["candidate_partition"]["path"]),
                phylo.parent / "orthohmm_working_res/orthohmm_edges_clustered.txt")
            report["outputs"] = [record(phylo / name) for name in PRODUCTS]
            native = replay_native(report, case)
            require_markers(native["summary"])
            report["status"] = "native_verified_scoring_pending"
            write_json(output / "results.json", report)
            verify(plan, [*checked, *native["checked"]])
            report.update(score_qfo(repo, phylo, Path(case["output"]), case["scoring_id"]))
            scores, details, files = readback_scores(repo, "QfO", report, native, Path(case["output"]),
                "ohmm_hillclimb_r1_" + case["scoring_id"])
            execution_kind = "new_cached_graph_component_inference_own_tree"
        else:
            native = replay_native(plan["parent_native_row"], plan["parent_native_case"])
            participant = ("ohmm_hillclimb_r1_" + plan["parent_native_case"]["scoring_id"]
                           if plan["dataset"] == "QfO" else "unused")
            scores, details, files = readback_scores(repo, plan["dataset"], plan["parent_native_row"], native,
                Path(plan["parent_native_case"]["output"]), participant)
            if endpoint_vector(plan["dataset"], scores) != plan["matched_control"]["scores"]:
                raise ValueError("Exact native-input reuse score readback differs from matched control")
            if (native["products"] != plan["matched_control"]["products"]
                    or native["summary"] != plan["matched_control"]["native_summary"]):
                raise ValueError("Exact native-input reuse products differ from admitted control")
            execution_kind = ("exact_inactive_candidate_reuse_no_native_rerun" if plan["dataset"] != "QfO"
                              else "exact_component_candidate_reuse_no_native_rerun")
        if native["genes"] != EXPECTED_GENES[plan["dataset"]] or native["species"] != (78 if plan["dataset"] == "QfO" else 12):
            raise ValueError("Component evaluation lost the full dataset universe")
        require_markers(native["summary"])
        checked.extend([*native["checked"], *files, *native["products"].values()])
        verify(plan, checked)
        report.update(status="independent_content_readback_verified", scores=endpoint_vector(plan["dataset"], scores),
            details=details, native_products=native["products"], native_summary=native["summary"],
            coverage={k: native[k] for k in ("genes", "species", "groups", "pairs")},
            execution_kind=execution_kind, inputs_reverified_unchanged=True)
    except Exception as error:
        report.update(status="failed_with_explicit_gap", error_type=type(error).__name__, error=str(error))
    report["host_after"] = host_activity()
    write_json(output / "results.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    prep = modes.add_parser("prepare")
    for name in ("repo", "source-root", "output", "branch-root", "comparison"):
        prep.add_argument("--" + name, required=True, type=Path)
    prep.add_argument("--dataset", required=True, choices=tuple(EXPECTED_GENES))
    worker = modes.add_parser("run")
    worker.add_argument("--execution-plan", required=True, type=Path)
    args = parser.parse_args()
    result = (prepare(args.repo.resolve(), args.source_root.resolve(), args.output.resolve(), args.dataset,
        args.branch_root.resolve(), args.comparison.resolve()) if args.mode == "prepare" else run(args.execution_plan.resolve()))
    print(json.dumps({"status": result["status"]}))
    return 0 if result["status"] in ("prepared_after_complete_branch_admission", "independent_content_readback_verified") else 1


if __name__ == "__main__":
    raise SystemExit(main())
