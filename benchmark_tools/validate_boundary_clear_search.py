"""Validate private boundary-clear libraries through complete native pair searches."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
from benchmark_tools.benchmark_viterbi_center_recurrence import array_records, bind_library
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_search_stage_cost import fingerprint_result, selected_names, subset_species
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.run_search_recovery import host_activity, runtime_records
from orthohmm.search import engine, viterbi
from orthohmm.search.sequences import SpeciesSequences

ARMS = ("full_clear", "boundary_only")
ORDERS = (ARMS, ARMS[::-1], ARMS, ARMS[::-1])
PROTOCOL = "hillclimb_boundary_clear_search_validation_20261014.json"
SETTINGS = {"queries_per_pair":512, "k":4, "matrix":"BLOSUM62", "band_width":64, "threads":4,
    "caps":[100,400], "phase_a_topk":5000, "use_reduced_alphabet":False, "min_total_hits":4,
    "min_diag_hits":1, "diag_bin_width":10, "query_orders":["fasta","reversed"],
    "timed_orders":[list(order) for order in ORDERS]}


def reordered_species(species, order):
    indices = list(order)
    if sorted(indices) != list(range(species.num_sequences)):
        raise ValueError("Query order must be an exact permutation")
    lengths = np.array(species.lengths[indices], dtype=np.int32)
    offsets = np.zeros(len(indices), dtype=np.int64)
    offsets[1:] = np.cumsum(lengths[:-1])
    return SpeciesSequences(species.species_file, [species.ids[i] for i in indices],
        np.concatenate([species.get_sequence(i) for i in indices]), offsets, lengths)


def canonical_result(result, query, target):
    indices = (result.query_indices, result.target_indices)
    if (result.candidate_count != len(result.scores) or len(result.evalues) != len(result.scores)
            or result.scores.ndim != 1 or result.evalues.ndim != 1
            or any(len(values) != len(result.scores) for values in indices)
            or any(values.ndim != 1 or not np.issubdtype(values.dtype, np.integer) for values in indices)
            or np.any(indices[0] < 0) or np.any(indices[0] >= query.num_sequences)
            or np.any(indices[1] < 0) or np.any(indices[1] >= target.num_sequences)
            or not np.isfinite(result.scores).all() or not np.isfinite(result.evalues).all()
            or np.any(result.evalues < 0)):
        raise ValueError("Invalid native search outputs")
    records = sorted((query.ids[int(q)], target.ids[int(t)], float(score).hex(), float(evalue).hex())
        for q,t,score,evalue in zip(*indices, result.scores, result.evalues))
    if len({row[:2] for row in records}) != len(records):
        raise ValueError("Duplicate native search gene pair")
    return hashlib.sha256(json.dumps(records, separators=(",", ":")).encode()).hexdigest()


def search_with_library(query, target, cap, library, capture_inputs=False):
    previous = viterbi._c_hmm_lib, viterbi._c_hmm_available
    viterbi._c_hmm_lib, viterbi._c_hmm_available = library, True
    original = engine.batch_viterbi_multipair_c
    captures = []
    def strict_kernel(*args, **kwargs):
        if len(args) != 10 or args[-1] != 64 or kwargs != {"n_threads":4}:
            raise ValueError("Native search changed prespecified kernel route/settings")
        before = array_records(args[:9]) if capture_inputs else None
        scores = original(*args, **kwargs)
        if scores.dtype != np.int32 or scores.ndim != 1 or len(scores) != len(args[8]):
            raise ValueError("Kernel returned invalid integer score vector")
        if capture_inputs and array_records(args[:9]) != before:
            raise ValueError("Kernel modified search input arrays")
        captures.append((scores, before))
        return scores
    def forbidden(*args, **kwargs):
        raise AssertionError("Native search silently fell back from the required multipair kernel")
    wall, cpu = time.monotonic(), time.process_time()
    try:
        with patch.object(engine, "batch_viterbi_multipair_c", strict_kernel), \
                patch.object(engine, "batch_viterbi_c", forbidden), \
                patch.object(engine, "batch_viterbi_score", forbidden), \
                patch.object(engine, "batch_viterbi_cuda", forbidden):
            result = engine.search_species_pair_indexed(query, target, "BLOSUM62", band_width=64,
                kmer_k=4, use_reduced_alphabet=False, min_total_hits=4, min_diag_hits=1,
                diag_bin_width=10, max_candidates_per_query=cap, n_threads=4, prefilter_phase_a_topk=5000)
        wall, cpu = time.monotonic() - wall, time.process_time() - cpu
    finally:
        viterbi._c_hmm_lib, viterbi._c_hmm_available = previous
    if len(captures) != 1:
        raise ValueError("Expected one actual multipair-kernel call per sampled search")
    integer_scores, inputs = captures[0]
    return {"wall_s": wall, "process_cpu_s": cpu, "fingerprint": fingerprint_result(result),
        "integer_score_sha256": hashlib.sha256(integer_scores.tobytes()).hexdigest(),
        "canonical_gene_pair_sha256": canonical_result(result, query, target),
        "significant_hits": int(np.count_nonzero(result.evalues < 1e-4)),
        "kernel_calls": 1, "kernel_inputs": inputs}


def output_identity(row):
    return {key: row[key] for key in ("fingerprint", "integer_score_sha256", "canonical_gene_pair_sha256", "significant_hits")}


def run(repo, output, source_revision):
    if output.exists() or os.environ.get("CUDA_VISIBLE_DEVICES") != "" or engine.is_cuda_available():
        raise ValueError("Require new isolated CPU-only search validation namespace")
    source = Path(__file__).resolve().parent.parent
    revision = subprocess.check_output(["git", "rev-parse", "--verify", source_revision + "^{commit}"], cwd=repo, text=True).strip()
    source_names = ["benchmark_tools/validate_boundary_clear_search.py", "benchmark_tools/" + PROTOCOL,
        "orthohmm/search/csrc/hmm_viterbi.c", *[str(p.relative_to(source)) for p in sorted((source / "orthohmm/search").rglob("*.py"))]]
    sources = [record(source / name) for name in source_names]
    for name, pin in zip(source_names, sources):
        blob = subprocess.check_output(["git", "show", revision + ":" + name], cwd=repo)
        if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
            raise ValueError("Wrong frozen search validation source")
    live_path = repo / "benchmarks/results/full_search_hillclimb_r1_v3_20261012/execution_plan.json"
    live = json.loads(live_path.read_text())
    case = next(c for c in live["cases"] if c["dataset"] == "QfO" and c["label"] == "control100")
    pins = {Path(pin["path"]).name: pin for pin in case["inputs"]}
    if len(pins) != 78 or len(pins) != len(case["inputs"]):
        raise ValueError("Frozen QfO proteome inventory changed")
    names = selected_names(list(pins), 3)
    checked = [record(live_path), *[pins[name] for name in names]]
    for name in ("kmer_prefilter.so", "pair_align.so"):
        relative = Path("orthohmm/search/csrc") / name
        old, copied = record(Path(live["source_root"]) / relative), record(source / relative)
        if old["sha256"] != copied["sha256"]:
            raise ValueError("Frozen support extension differs from live search")
        checked.extend((old, copied))
    decision_path = repo / "benchmark_tools/results/viterbi_boundary_clear_terminal_decision_20261013.json"
    decision = json.loads(decision_path.read_text())
    if decision["status"] != "replicated_exact_kernel_gain_eligible_for_broader_output_validation_only":
        raise ValueError("Previous exact-kernel validation not admitted")
    previous_c = subprocess.check_output(["git", "show", decision["source_revision"] + ":orthohmm/search/csrc/hmm_viterbi.c"], cwd=repo)
    if hashlib.sha256(previous_c).hexdigest() != record(source / "orthohmm/search/csrc/hmm_viterbi.c")["sha256"]:
        raise ValueError("Kernel changed after the admitted microbenchmark")
    checked.append(record(decision_path))
    for pin in (*checked, *sources):
        check(pin)
    protocol = json.loads((source / "benchmark_tools" / PROTOCOL).read_text())
    if protocol["settings"] != SETTINGS:
        raise ValueError("Search validation differs from prespecified cohort/orders")
    output.mkdir(parents=True)
    runtime = runtime_records()
    report = {"status": "boundary_clear_search_validation_running", "source_revision": revision,
        "job_id": os.environ.get("SLURM_JOB_ID"), "sources": sources, "checked": checked,
        "runtime": runtime, "host_before": host_activity(), "cases": [], "settings": protocol["settings"],
        "truth_labels_loaded": False, "defaults_promoted": False, "goal_completed": False,
        "scope": "Warmed complete species-pair search on sampled queries against full targets, not full orthology inference or isolated timing. Input/output fingerprint checks are excluded from timed calls except required route tracking."}
    write_json(output / "report.json", report)
    started = time.monotonic()
    try:
        compiler = shutil.which("gcc")
        if compiler is None:
            raise FileNotFoundError("GCC is unavailable")
        libraries, library_pins, commands = {}, {}, {}
        for arm in ARMS:
            path = output / (arm + ".so")
            flag = ["-DHMM_VITERBI_BOUNDARY_CLEAR_ONLY=1"] if arm == "boundary_only" else []
            command = [compiler, "-O3", "-fopenmp", "-shared", "-fPIC", "-march=native", *flag,
                       str(source / "orthohmm/search/csrc/hmm_viterbi.c"), "-o", str(path)]
            subprocess.run(command, check=True, capture_output=True)
            libraries[arm], library_pins[arm], commands[arm] = bind_library(path), record(path), command
        report.update(kernels=library_pins, compiler=record(compiler), compiler_commands=commands)
        loaded = {name: SpeciesSequences.from_fasta(pins[name]["path"], name) for name in names}
        order_identities = {}
        for index, name in enumerate(names):
            subset, indices = subset_species(loaded[name], 512)
            target = loaded[names[(index + 1) % 3]]
            for query_order in ("fasta", "reversed"):
                query = subset if query_order == "fasta" else reordered_species(subset, range(511,-1,-1))
                sequence_arrays = (query.flat_sequences, query.offsets, query.lengths, target.flat_sequences, target.offsets, target.lengths)
                sequence_before = array_records(sequence_arrays)
                for cap in (100,400):
                    row = {"query": name, "target": target.species_file, "query_order": query_order, "cap": cap,
                        "query_ids": query.ids, "query_indices": indices if query_order == "fasta" else indices[::-1],
                        "target_genes": target.num_sequences, "input_sequences": sequence_before, "warm": {}, "repeats": []}
                    report["cases"].append(row)
                    for arm in ARMS:
                        row["warm"][arm] = search_with_library(query, target, cap, libraries[arm], capture_inputs=True)
                    reference = output_identity(row["warm"][ARMS[0]])
                    if (output_identity(row["warm"][ARMS[1]]) != reference
                            or row["warm"][ARMS[0]]["kernel_inputs"] != row["warm"][ARMS[1]]["kernel_inputs"]):
                        raise ValueError("Boundary clearing changed complete warm search outputs/inputs")
                    order_key = (name,cap)
                    if query_order == "fasta":
                        order_identities[order_key] = reference["canonical_gene_pair_sha256"]
                    elif order_identities[order_key] != reference["canonical_gene_pair_sha256"]:
                        raise ValueError("Canonical gene-pair outputs changed across query orders")
                    for repeat, order in enumerate(ORDERS):
                        entry = {"repeat":repeat, "order":list(order), "arms":{}}
                        for arm in order:
                            actual = search_with_library(query, target, cap, libraries[arm])
                            if output_identity(actual) != reference:
                                raise ValueError("Boundary clearing changed complete timed search outputs")
                            entry["arms"][arm] = actual
                        row["repeats"].append(entry)
                        write_json(output / "report.json", report)
                        if time.monotonic() - started > 35*60:
                            raise TimeoutError("Prespecified native-search wall limit exceeded")
                    row["ratios"] = {field:float(np.median([r["arms"][ARMS[1]][field]/r["arms"][ARMS[0]][field] for r in row["repeats"]]))
                                     for field in ("wall_s", "process_cpu_s")}
                    row["outputs_and_order_identity_preserved"] = True
                    if array_records(sequence_arrays) != sequence_before:
                        raise ValueError("Search modified query/target sequence arrays")
        for pin in (*checked, *sources, *runtime["files"], *library_pins.values(), report["compiler"]):
            check(pin)
        if runtime_records() != runtime:
            raise ValueError("Scientific runtime changed during search validation")
        usage = resource.getrusage(resource.RUSAGE_SELF)
        if usage.ru_maxrss > 16*1024**2:
            raise MemoryError("Prespecified native-search memory limit exceeded")
        report.update(status="complete_boundary_clear_native_search_output_equivalence",
            resources={"wall_s":time.monotonic()-started, "cpu_s":usage.ru_utime+usage.ru_stime, "peak_rss_kib":usage.ru_maxrss})
    except Exception as error:
        report.update(status="failed_boundary_clear_native_search_validation", error_type=type(error).__name__, error=str(error))
    report["host_after"] = host_activity()
    write_json(output / "report.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-revision", required=True)
    args = parser.parse_args()
    result = run(args.repo.resolve(), args.output.resolve(), args.source_revision)
    print(result["status"])
    raise SystemExit(0 if result["status"] == "complete_boundary_clear_native_search_output_equivalence" else 1)
