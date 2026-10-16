"""Require fresh protected-union HMM results to match the cached feasibility proof."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_search_recovery import read_raw
from benchmark_tools.diagnose_protected_candidate_union import context, protected_union
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.run_search_recovery import host_activity, runtime_records, source_records
from benchmark_tools.search_decision_trace import classify_search_result
from benchmark_tools.strict_native_cli import install_guard
from benchmark_tools.trace_ob_search_decisions import load_watched, subset_queries

PROTOCOL = "hillclimb_protected_union_fresh_pilot_20261016.json"
SLURM = "fresh_protected_union_search.slurm"
LIBRARIES = ("hmm_viterbi.so", "kmer_prefilter.so", "pair_align.so")
FIELDS = ("query_indices", "target_indices", "scores", "evalues")
SUCCESS = "fresh_control_and_protected_union_all_144_directions_verified_pending_independent_readback"


def require_same(actual, expected, queries, targets):
    for result in (actual, expected):
        classify_search_result(result, queries, targets, [], 0.0001)
    for field in FIELDS:
        a, b = getattr(actual, field), getattr(expected, field)
        floating = field in ("scores", "evalues")
        dtype = np.dtype("float64" if floating else "int32")
        if a.dtype != dtype or b.dtype != dtype or a.shape != b.shape:
            raise ValueError("Fresh/cached candidate layout differs: " + field)
        if not np.array_equal(a.view(np.uint64) if floating else a, b.view(np.uint64) if floating else b):
            raise ValueError("Fresh/cached candidate values differ: " + field)
    return True


def require_routes(directory, library, pairs, calls):
    paths = sorted(directory.glob("routes_*.json"))
    if len(paths) != 1:
        raise ValueError("Require one observed diagnostic CPU kernel process")
    report = json.loads(paths[0].read_text())
    if (report["library"] != library or report["routes"] != {
            "initial_search": {"calls": calls, "scored_pairs": pairs}}):
        raise ValueError("Actual kernel routes do not match all fresh candidates")
    check(library)
    return record(paths[0])


def prepare(repo, source, output):
    if source == repo or output.exists():
        raise ValueError("Require a separate frozen executor and fresh diagnostic namespace")
    _, execution, reports, original_pins, revision = context(repo, source)
    protocol = json.loads((source / "benchmark_tools" / PROTOCOL).read_text())
    prerequisite = {k: v for k, v in protocol["prerequisite"].items() if k != "status"}
    prerequisite["path"] = str(repo / prerequisite["path"])
    check(prerequisite)
    receipt = json.loads(Path(prerequisite["path"]).read_text())
    if receipt["status"] != protocol["prerequisite"]["status"] or receipt["advance_gate_passed"] is not True:
        raise ValueError("Require independently verified cached-union feasibility")
    check(receipt["producer"])
    cached = json.loads(Path(receipt["producer"]["path"]).read_text())
    sources = [*source_records(source), *[record(source / "benchmark_tools" / name) for name in
        ("fresh_protected_union_search.py", SLURM, PROTOCOL, "diagnose_protected_candidate_union.py", "strict_native_cli.py")]]
    for pin in sources:
        relative = Path(pin["path"]).relative_to(source)
        if relative.suffix != ".so":
            blob = subprocess.check_output(["git", "-C", str(source), "show", f"{revision}:{relative}"])
            if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
                raise ValueError("Fresh diagnostic differs from committed frozen source")
    libraries = []
    for name in LIBRARIES:
        originals = [p for p in execution["sources"] if Path(p["path"]).name == name]
        actual = record(source / "orthohmm/search/csrc" / name)
        if len(originals) != 1 or (actual["bytes"], actual["sha256"]) != (originals[0]["bytes"], originals[0]["sha256"]):
            raise ValueError("CPU library differs from original control: " + name)
        libraries.append({"original": originals[0], "copy": actual})
    inputs = [prerequisite, receipt["producer"], *original_pins, *cached["checked"],
              *[record(Path(execution["output"]) / r["label"] / "report.json") for r in reports]]
    for pin in (*inputs, *sources):
        check(pin)
    plan = {"status": "prepared_fresh_protected_union_search", "repo": str(repo), "output": str(output),
        "source_root": str(source), "source_revision": revision, "runtime": execution["runtime"],
        "original_execution": reports[0]["execution_plan"], "original_reports": inputs[-2:],
        "inputs": list({p["path"]: p for p in inputs}.values()), "sources": sources, "libraries": libraries,
        "threads": 8, "query_scope": "Reference-conditioned diagnostic queries against complete target proteomes, not accuracy or full-pipeline timing.",
        "defaults_promoted": False, "goal_completed": False}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", plan)
    return plan


def run(path):
    plan_pin = record(path)
    plan = json.loads(path.read_text())
    repo, source, output = (Path(plan[k]) for k in ("repo", "source_root", "output"))
    if (Path(__file__).resolve() != source / "benchmark_tools/fresh_protected_union_search.py"
            or path.resolve() != output / "execution_plan.json" or (output / "results.json").exists()):
        raise ValueError("Wrong frozen executor or previously attempted namespace")
    if (os.environ.get("SLURM_CPUS_PER_TASK"), os.environ.get("SLURM_MEM_PER_NODE"),
            os.environ.get("SLURMD_NODENAME"), os.environ.get("CUDA_VISIBLE_DEVICES")) != ("8", "65536", "bizon", ""):
        raise ValueError("Require CPU-only 8-CPU 64-GiB local allocation")
    report = {"status": "validating_fresh_protected_union_search", "execution_plan": plan_pin,
        "source_revision": plan["source_revision"], "job_id": os.environ.get("SLURM_JOB_ID"), "directions": [],
        "accuracy_evaluated": False, "defaults_promoted": False, "goal_completed": False,
        "query_scope": plan["query_scope"], "host_before": host_activity()}
    write_json(output / "results.json", report)
    checked = [plan_pin, *plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]
    started = time.monotonic()
    try:
        for pin in checked:
            check(pin)
        if runtime_records() != plan["runtime"]:
            raise ValueError("Original scientific runtime changed")
        revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
        if revision != plan["source_revision"] or plan["threads"] != 8:
            raise ValueError("Frozen executor revision or thread allocation changed")
        from orthohmm.search import engine, profile_expansion, viterbi
        from orthohmm.search.prefilter import is_c_prefilter_available
        if Path(engine.__file__).resolve() != source / "orthohmm/search/engine.py" or not is_c_prefilter_available():
            raise ValueError("Wrong scientific package or missing compiled prefilter")
        library = install_guard(engine, profile_expansion, viterbi, output / "kernel_routes")
        if library != next(p["copy"] for p in plan["libraries"] if Path(p["copy"]["path"]).name == "hmm_viterbi.so"):
            raise ValueError("Actual loaded kernel differs from frozen original control library")
        execution = json.loads(Path(plan["original_execution"]["path"]).read_text())
        retained = [json.loads(Path(p["path"]).read_text()) for p in plan["original_reports"]]
        species, owners = {}, {}
        for pin in execution["fastas"]:
            filename = Path(pin["path"]).name
            sp = engine.SpeciesSequences.from_fasta(pin["path"], filename)
            species[filename] = sp
            for gene in sp.ids:
                if gene in owners:
                    raise ValueError("Duplicate original FASTA gene")
                owners[gene] = filename
        if len(species) != 12 or len(owners) != 251378:
            raise ValueError("Wrong complete diagnostic input universe")
        watched = load_watched(execution["pair_trace"]["path"])
        directions = {}
        for pair in watched:
            directions.setdefault((owners[pair[0]], owners[pair[1]]), []).append(pair)
        expected = sorted(directions)
        if len(expected) != 144 or len(watched) != 81466 or any(
                [(d["query_species"], d["target_species"]) for d in r["directions"]] != expected for r in retained):
            raise ValueError("Wrong diagnostic direction inventory")
        candidates, kernel_calls = 0, 0
        for i, key in enumerate(expected):
            query = subset_queries(species[key[0]], {q for q, _ in directions[key]})
            target = species[key[1]]
            current = []
            row = {"query_species": key[0], "target_species": key[1], "query_count": len(query.ids), "target_count": len(target.ids)}
            # The retained candidate arrays are read only after both fresh searches.
            for label in ("control", "protected_union"):
                settings = {**retained[0]["settings"], **({"candidate_policy": "protected_union"} if label != "control" else {})}
                result = engine.search_species_pair_indexed(query, target, n_threads=plan["threads"], **settings)
                classify_search_result(result, query.ids, target.ids, [], 0.0001)
                raw = output / f"{label}_direction_{i:03d}.npz"
                np.savez_compressed(raw, query_ids=np.array(query.ids), target_ids=np.array(target.ids),
                    **{field: getattr(result, field) for field in FIELDS}, candidate_count=np.array(result.candidate_count))
                current.append(result)
                row[label] = {"raw": record(raw), "candidate_count": result.candidate_count}
                candidates += result.candidate_count
                kernel_calls += bool(result.candidate_count)
            originals = []
            for old in retained:
                old_result, queries, targets = read_raw(old["directions"][i]["raw"]["path"])
                if queries != query.ids or targets != target.ids:
                    raise ValueError("Fresh/cached query or complete target order differs")
                originals.append(old_result)
            expected_union, union_stats = protected_union(*originals, query.ids, target.ids)
            require_same(current[0], originals[0], query.ids, target.ids)
            require_same(current[1], expected_union, query.ids, target.ids)
            row.update(fresh_control_exact=True, fresh_union_exact=True, expected_union_counts=union_stats)
            report["directions"].append(row)
            report["status"] = "fresh_protected_union_search_in_progress"
            write_json(output / "results.json", report)
            print(f"Fresh protected union: {i + 1}/{len(expected)} directions verified", flush=True)
        routes = require_routes(output / "kernel_routes", library, candidates, kernel_calls)
        checked += [routes, library, *[row[label]["raw"] for row in report["directions"] for label in ("control", "protected_union")]]
        report.update(status=SUCCESS, actual_kernel_routes=routes, scored_candidates=candidates,
            actual_kernel_calls=kernel_calls, fresh_initial_search_executed=True,
            scope="Fresh reference-conditioned search only. All saved values match authenticated retained control/cached-union arrays; no full-proteome orthology accuracy claim.")
    except Exception as error:
        report.update(status="failed_fresh_protected_union_search", error_type=type(error).__name__, error=str(error))
    finally:
        try:
            for pin in checked:
                check(pin)
            if runtime_records() != plan["runtime"]:
                raise ValueError("Runtime changed during fresh diagnostic")
            report["postflight_pins_verified"] = True
        except Exception as error:
            report.update(status="failed_fresh_protected_union_search", postflight_error=str(error))
        report.update(total_diagnostic_wall_s=time.monotonic() - started, host_after=host_activity())
        write_json(output / "results.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="operation", required=True)
    preparation = sub.add_parser("prepare")
    for name in ("repo", "source", "output"):
        preparation.add_argument("--" + name, type=Path, required=True)
    execution = sub.add_parser("run")
    execution.add_argument("--execution-plan", type=Path, required=True)
    args = parser.parse_args()
    result = (prepare(args.repo.resolve(), args.source.resolve(), args.output.resolve())
              if args.operation == "prepare" else run(args.execution_plan.resolve()))
    print(result["status"])
    raise SystemExit(0 if result["status"] in ("prepared_fresh_protected_union_search", SUCCESS) else 1)
