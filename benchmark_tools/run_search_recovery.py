"""Fresh, reference-conditioned HMM candidate-recovery diagnostics."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import importlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import tempfile
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.search_decision_trace import classify_search_result
from benchmark_tools.trace_ob_search_decisions import load_watched, subset_queries

PLAN_NAME = "hillclimb_search_recovery_20261012.json"


def resolve_settings(plan, label):
    arms = {arm["label"]: arm for arm in plan["arms"]}
    if len(arms) != len(plan["arms"]) or label not in arms:
        raise ValueError("Unknown or duplicate search arm")
    overrides = arms[label]["overrides"]
    if set(overrides) - set(plan["base_settings"]):
        raise ValueError("Unknown search parameter")
    return {**plan["base_settings"], **overrides}


def link_changes(control, observed):
    def accepted(rows):
        pairs = [(row["query"], row["target"]) for row in rows]
        if len(pairs) != len(set(pairs)):
            raise ValueError("Duplicate diagnostic pair")
        if any(row["decision"] not in {
            "accepted", "scored_not_significant", "not_selected_by_prefilter"
        } for row in rows):
            raise ValueError("Unknown diagnostic decision")
        return set(pairs), {(row["query"], row["target"])
                            for row in rows if row["decision"] == "accepted"}
    left_universe, left = accepted(control)
    right_universe, right = accepted(observed)
    if left_universe != right_universe:
        raise ValueError("Diagnostic pair universes differ")
    return {"gained": len(right - left), "lost": len(left - right),
            "retained": len(left & right), "absent_both": len(left_universe - left - right)}


def source_records(source):
    paths = [source / "benchmark_tools" / name for name in (
        "run_search_recovery.py", "run_search_recovery.slurm", "audit_search_recovery.py",
        "search_decision_trace.py", "trace_ob_search_decisions.py",
        "prepare_ob_candidate_neighborhood.py", "run_membership_hillclimb.py")]
    paths += sorted((source / "orthohmm").rglob("*.py"))
    paths += sorted((source / "orthohmm/search/csrc").glob("*.c"))
    paths += sorted((source / "orthohmm/search/csrc").glob("*.so"))
    for name in ("kmer_prefilter.so", "hmm_viterbi.so"):
        if not (source / "orthohmm/search/csrc" / name).is_file():
            raise ValueError("Require compiled CPU extension: " + name)
    return [record(path) for path in paths]


def runtime_records():
    versions, files = {}, [record(sys.executable)]
    for name in ("numpy", "numba", "llvmlite", "Bio"):
        module = importlib.import_module(name)
        versions[name] = module.__version__
        files.append(record(module.__file__))
    return {"python": sys.version, "versions": versions, "files": files}


def prepare(repo, source, output, threads):
    if output.exists() or not 1 <= threads <= 8:
        raise ValueError("Require fresh output and 1-8 threads")
    revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    subprocess.run(["git", "-C", str(source), "diff", "--exit-code", "HEAD", "--",
                    "orthohmm", "benchmark_tools/run_search_recovery.py", "benchmark_tools/" + PLAN_NAME], check=True)
    plan_path = source / "benchmark_tools" / PLAN_NAME
    plan = json.loads(plan_path.read_text())
    trace_path = repo / plan["input_manifest"]
    if record(trace_path)["sha256"] != plan["input_manifest_sha256"]:
        raise ValueError("Changed reference manifest")
    trace = json.loads(trace_path.read_text())
    fastas = sorted([item for item in trace["inputs"] if Path(item["path"]).suffix == ".fa"],
                    key=lambda item: item["path"])
    if len(fastas) != 12:
        raise ValueError("Wrong FASTA inventory")
    inputs = [record(plan_path), record(trace_path), trace["pair_trace"], *fastas]
    historical_path = repo / "benchmark_tools/results/ob_search_decisions_audit_20260918.json"
    historical = None
    if historical_path.exists():
        audit = json.loads(historical_path.read_text())
        if audit["status"] != "search_decisions_independently_recounted":
            raise ValueError("Historical control has not been independently read back")
        historical = audit["source_report"]
        inputs.extend([record(historical_path), historical])
    for item in inputs:
        check(item)
    for arm in plan["arms"]:
        resolve_settings(plan, arm["label"])
    report = {"schema_version": 1, "status": "prepared_unrun", "source_revision": revision,
              "source_root": str(source), "repo": str(repo), "output": str(output),
              "threads": threads, "plan": plan, "inputs": inputs, "fastas": fastas,
              "pair_trace": trace["pair_trace"], "sources": source_records(source),
              "historical_control": historical,
              "runtime": runtime_records(), "accuracy_evaluated": False,
              "publication_baseline_replaced": False}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", report)
    return report


def assert_control_fixture(engine, sequences):
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "fixture.fa"
        path.write_text(">same\n" + "ACDEFGHIKLMNPQRSTVWY" * 4 + "\n>biased\n" + "A" * 80 + "\n")
        species = sequences.SpeciesSequences.from_fasta(str(path), "fixture.fa")
        implicit = engine.search_species_pair_indexed(species, species, "BLOSUM62",
                    max_candidates_per_query=100, n_threads=1)
        explicit = engine.search_species_pair_indexed(species, species, "BLOSUM62",
                    max_candidates_per_query=100, n_threads=1, prefilter_phase_a_topk=5000)
        for field in ("query_indices", "target_indices", "scores", "evalues"):
            if not np.array_equal(getattr(implicit, field), getattr(explicit, field)):
                raise ValueError("Explicit control differs from production fixture")
        if implicit.candidate_count != explicit.candidate_count:
            raise ValueError("Explicit control differs from production fixture")
    return True


def host_activity():
    return subprocess.check_output(["ps", "-eo", "pid,pcpu,comm"], text=True)


def run(execution_path, label):
    execution_identity = record(execution_path)
    execution = json.loads(execution_path.read_text())
    source = Path(execution["source_root"])
    if Path(__file__).resolve() != source / "benchmark_tools/run_search_recovery.py":
        raise ValueError("Wrong execution source")
    checked = [execution_identity, *execution["inputs"], *execution["sources"],
               *execution["runtime"]["files"]]
    for item in checked:
        check(item)
    if runtime_records() != execution["runtime"]:
        raise ValueError("Runtime differs from preparation")
    if os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Hide CUDA before starting the CPU diagnostic interpreter")
    from orthohmm.search import engine, sequences
    from orthohmm.search.prefilter import is_c_prefilter_available
    if Path(engine.__file__).resolve() != source / "orthohmm/search/engine.py":
        raise ValueError("Wrong scientific engine")
    if engine.is_cuda_available() or not engine.is_c_available() or not is_c_prefilter_available():
        raise ValueError("Require CPU HMM and C extension-ranking backends")
    settings = resolve_settings(execution["plan"], label)
    folder = Path(execution["output"]) / label
    if folder.exists():
        raise FileExistsError(folder)
    folder.mkdir()
    started, usage = time.monotonic(), resource.getrusage(resource.RUSAGE_SELF)
    report = {"status": "running", "label": label, "settings": settings,
              "execution_plan": execution_identity, "source_revision": execution["source_revision"],
              "job_id": os.environ.get("SLURM_JOB_ID"), "directions": [],
              "accuracy_evaluated": False, "defaults_promoted": False,
              "limitations": [execution["plan"]["query_scope"],
                  execution["plan"]["resources"]["timing_scope"],
                  "Unlabelled hits are not false positives; diagnostic link recovery is not orthology F1."]}
    try:
        report["host_activity_before"] = host_activity()
        report["explicit_control_fixture_exact"] = assert_control_fixture(engine, sequences)
        watched = load_watched(execution["pair_trace"]["path"])
        species, owners = {}, {}
        for item in execution["fastas"]:
            path = Path(item["path"])
            sp = sequences.SpeciesSequences.from_fasta(str(path), path.name)
            species[path.name] = sp
            for gene in sp.ids:
                if gene in owners:
                    raise ValueError("Duplicate input gene")
                owners[gene] = path.name
        if len(owners) != 251378 or any(g not in owners for pair in watched for g in pair):
            raise ValueError("Wrong full input universe")
        directions = {}
        for pair in watched:
            directions.setdefault((owners[pair[0]], owners[pair[1]]), []).append(pair)
        for index, ((qsp, tsp), pairs) in enumerate(sorted(directions.items())):
            pairs = sorted(pairs)
            query = subset_queries(species[qsp], {q for q, _ in pairs})
            target = species[tsp]
            result = engine.search_species_pair_indexed(query, target,
                       n_threads=execution["threads"], **settings)
            raw = folder / f"direction_{index:03d}.npz"
            np.savez_compressed(raw, query_ids=np.array(query.ids), target_ids=np.array(target.ids),
                query_indices=result.query_indices, target_indices=result.target_indices,
                scores=result.scores, evalues=result.evalues, candidate_count=np.array(result.candidate_count))
            thresholds = []
            for threshold in execution["plan"]["significance_thresholds"]:
                rows = classify_search_result(result, query.ids, target.ids, pairs, threshold)
                table = folder / f"direction_{index:03d}_e{threshold:g}.tsv"
                counts, recovered = Counter(), 0
                with table.open("x", newline="") as stream:
                    writer = csv.DictWriter(stream, delimiter="\t", fieldnames=[
                        "query", "target", "decision", "score", "evalue", "historical_present"])
                    writer.writeheader()
                    for row in rows:
                        row["historical_present"] = watched[row["query"], row["target"]]["historical_present"]
                        counts[row["decision"]] += 1
                        recovered += row["decision"] == "accepted" and not row["historical_present"]
                        writer.writerow(row)
                thresholds.append({"threshold": threshold, "decisions": dict(counts),
                    "historically_absent_accepted": recovered,
                    "all_significant_hits": int(np.count_nonzero(result.evalues < threshold)),
                    "table": record(table)})
            report["directions"].append({"query_species": qsp, "target_species": tsp,
                "query_count": len(query.ids), "target_count": len(target.ids),
                "candidate_count": int(result.candidate_count), "watched_pairs": len(pairs),
                "raw": record(raw), "thresholds": thresholds})
            write_json(folder / "progress.json", report)
            print(f"{label}: direction {index + 1}/{len(directions)}", flush=True)
        for item in checked:
            check(item)
        if runtime_records() != execution["runtime"]:
            raise ValueError("Runtime changed during search")
        report["host_activity_after"] = host_activity()
        report["status"] = "fresh_search_completed_pending_independent_audit"
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        after = resource.getrusage(resource.RUSAGE_SELF)
        report["resources"] = {"wall_s": time.monotonic() - started,
            "cpu_s": after.ru_utime + after.ru_stime - usage.ru_utime - usage.ru_stime,
            "peak_rss_GiB": after.ru_maxrss / 1024 ** 2}
        write_json(folder / "report.json", report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    preparation = commands.add_parser("prepare")
    for name in ("repo", "source", "output"):
        preparation.add_argument("--" + name, type=Path, required=True)
    preparation.add_argument("--threads", type=int, default=8)
    execution = commands.add_parser("run")
    execution.add_argument("--execution-plan", type=Path, required=True)
    execution.add_argument("--label", required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args.repo.resolve(), args.source.resolve(), args.output.resolve(), args.threads)
    else:
        run(args.execution_plan.resolve(), args.label)


if __name__ == "__main__":
    main()
