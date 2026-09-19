"""Observe frozen CPU search decisions for every retained reference-family pair."""

import argparse
from collections import Counter
import csv
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.build_publication_runtime import verify_runtime, COMMIT
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.search_decision_trace import classify_search_result
from benchmark_tools.trace_ob_initial_edges import TRACE_SHA

RUNTIME_SHA = "aebea83807356b02307473506fa30c2dbd2c511d7ba75a0655eb12180a474d74"
SETTINGS = dict(matrix_name="BLOSUM62", kmer_k=4, min_total_hits=4,
                min_diag_hits=1, diag_bin_width=10, max_candidates_per_query=100,
                band_width=64, use_reduced_alphabet=False)


def subset_queries(species, wanted):
    indices = np.array([i for i, gene in enumerate(species.ids) if gene in wanted], dtype=np.int64)
    lengths = species.lengths[indices].copy()
    offsets = np.zeros(len(indices), dtype=np.int64)
    offsets[1:] = np.cumsum(lengths[:-1])
    flat = (np.concatenate([species.get_sequence(i) for i in indices]) if len(indices)
            else np.empty(0, dtype=species.flat_sequences.dtype))
    return type(species)(species.species_file, [species.ids[i] for i in indices],
                         flat, offsets, lengths)


def load_watched(path):
    pairs, families = {}, set()
    with Path(path).open(newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            families.add(row["refog"])
            for q, t, field in ((row["left"], row["right"], "forward_normalized_hit"),
                                (row["right"], row["left"], "reverse_normalized_hit")):
                pair, present = (q, t), row[field] != "NA"
                if pair in pairs and pairs[pair]["historical_present"] != present:
                    raise ValueError("Conflicting historical hit presence")
                entry = pairs.setdefault(pair, {"historical_present": present, "families": set()})
                entry["families"].add(row["refog"])
    if len(families) != 70 or len({g for pair in pairs for g in pair}) != 1944:
        raise ValueError("Wrong reference-family universe")
    return pairs


def run(root, core, output, threads):
    if output.exists() or threads < 1:
        raise ValueError("Require fresh output and positive thread count")
    if subprocess.check_output(["git", "-C", str(core), "rev-parse", "HEAD"], text=True).strip() != COMMIT:
        raise ValueError("Wrong frozen core")
    subprocess.run(["git", "-C", str(core), "diff", "--exit-code", "HEAD", "--", "orthohmm"], check=True)
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    read_frozen(runtime_path, RUNTIME_SHA)
    runtime = verify_runtime(runtime_path, core)
    trace_path = root / "benchmark_tools/results/ob_family_trace_verified_20260916.json"
    trace = read_frozen(trace_path, TRACE_SHA)
    fastas = sorted((r for r in trace["inputs"] if Path(r["path"]).suffix == ".fa"), key=lambda r: r["path"])
    if len(fastas) != 12:
        raise ValueError("Wrong FASTA inventory")
    checked = [record(__file__), record(Path(__file__).with_name("search_decision_trace.py")),
               record(trace_path), record(runtime_path), trace["pair_trace"], *fastas,
               *[record(p) for p in sorted((core / "orthohmm").rglob("*.py"))]]
    for item in checked:
        check(item)
    watched = load_watched(trace["pair_trace"]["path"])
    if any(name == "orthohmm" or name.startswith("orthohmm.") for name in sys.modules):
        raise ValueError("Require fresh interpreter before frozen-core import")
    sys.path.insert(0, str(core))
    engine = importlib.import_module("orthohmm.search.engine")
    sequences = importlib.import_module("orthohmm.search.sequences")
    if Path(engine.__file__).resolve() != core / "orthohmm/search/engine.py" or engine.is_cuda_available():
        raise ValueError("Wrong imported engine or non-CPU runtime")
    species, owners = {}, {}
    for item in fastas:
        path = Path(item["path"])
        sp = sequences.SpeciesSequences.from_fasta(str(path), path.name)
        species[path.name] = sp
        for gene in sp.ids:
            if gene in owners:
                raise ValueError("Duplicate input gene")
            owners[gene] = path.name
    if len(owners) != 251378 or any(g not in owners for pair in watched for g in pair):
        raise ValueError("Wrong input sequence universe")
    directions = {}
    for pair in watched:
        directions.setdefault((owners[pair[0]], owners[pair[1]]), []).append(pair)
    output.mkdir(parents=True)
    report = {"status": "running", "source": record(__file__), "core_commit": COMMIT,
              "runtime": record(runtime_path), "settings": SETTINGS, "threads": threads,
              "job_id": os.environ.get("SLURM_JOB_ID"), "checked_records": checked,
              "directions": [], "accuracy_evaluated": False, "publication_ready": False,
              "limitations": ["Contemporary frozen CPU observation, not historical execution authentication.",
                  "Queries restricted to all reference-family members; full target proteomes and fixed cap retained.",
                  "Historical comparison is hit presence only, not numerical score equivalence.",
                  "No counterfactual scoring of candidates excluded by the prefilter.",
                  "Shared-host diagnostic is not a comparative timing run."]}
    try:
        for index, ((qsp, tsp), pairs) in enumerate(sorted(directions.items())):
            pairs = sorted(pairs)
            query = subset_queries(species[qsp], {q for q, _ in pairs})
            target = species[tsp]
            result = engine.search_species_pair_indexed(query, target, n_threads=threads, **SETTINGS)
            rows = classify_search_result(result, query.ids, target.ids, pairs, 1e-4)
            raw = output / f"direction_{index:03d}.npz"
            np.savez_compressed(raw, query_ids=np.array(query.ids), target_ids=np.array(target.ids),
                                query_indices=result.query_indices, target_indices=result.target_indices,
                                scores=result.scores, evalues=result.evalues,
                                candidate_count=np.array(result.candidate_count))
            table = output / f"direction_{index:03d}.tsv"
            counts, mismatches = Counter(), 0
            with table.open("x", newline="") as stream:
                writer = csv.DictWriter(stream, delimiter="\t", fieldnames=[
                    "query", "target", "decision", "score", "evalue", "historical_present", "families"])
                writer.writeheader()
                for row in rows:
                    prior = watched[row["query"], row["target"]]
                    row.update(historical_present=prior["historical_present"], families=",".join(sorted(prior["families"])))
                    counts[row["decision"]] += 1
                    mismatches += (row["decision"] == "accepted") != row["historical_present"]
                    writer.writerow(row)
            report["directions"].append({"query_species": qsp, "target_species": tsp,
                "query_count": len(query.ids), "target_count": len(target.ids),
                "candidate_count": result.candidate_count, "watched_pairs": len(rows),
                "decisions": dict(counts), "historical_presence_mismatches": mismatches,
                "raw": record(raw), "table": record(table)})
            print(f"Completed direction {index + 1}/{len(directions)}", flush=True)
        for item in checked:
            check(item)
        if verify_runtime(runtime_path, core) != runtime:
            raise ValueError("Runtime changed during diagnostic")
        report["status"] = "search_decisions_observed_pending_independent_audit"
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with (output / "report.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "core", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    run(args.root.resolve(), args.core.resolve(), args.output.resolve(), args.threads)
