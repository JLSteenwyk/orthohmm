"""Measure label-free native search stages without changing inference behavior."""

import argparse
from contextlib import ExitStack
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from orthohmm.search import engine
from orthohmm.search.sequences import SpeciesSequences
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.run_search_recovery import host_activity, runtime_records

PROTOCOL = "hillclimb_search_stage_probe_20261013.json"
STAGES = ("prepare_kmer_index", "prefilter_candidates", "build_profiles_batch",
          "batch_viterbi_multipair_c", "batch_viterbi_c", "batch_viterbi_score", "batch_estimate_evalues")
SETTINGS = {"matrix": "BLOSUM62", "band_width": 64, "kmer_k": 3, "prefilter_phase_a_topk": 5000,
            "min_total_hits": 4, "min_diag_hits": 1, "diag_bin_width": 10, "candidate_caps": [100, 400],
            "evalue_threshold": 0.0001, "threads": 4, "cuda_visible_devices": ""}


def selected_names(names, count):
    if type(count) is not int or count <= 0 or len(names) < count or len(names) != len(set(names)):
        raise ValueError("Require enough unique identifiers and a positive integer sample size")
    return sorted(names, key=lambda name: (hashlib.sha256(name.encode()).hexdigest(), name))[:count]


def subset_species(species, count):
    selected = set(selected_names(species.ids, count))
    indices = [i for i, name in enumerate(species.ids) if name in selected]
    lengths = np.array(species.lengths[indices], dtype=np.int32)
    if np.any(lengths <= 0):
        raise ValueError("Empty sampled protein")
    offsets = np.zeros(count, dtype=np.int64)
    offsets[1:] = np.cumsum(lengths[:-1])
    subset = SpeciesSequences(species.species_file, [species.ids[i] for i in indices],
        np.concatenate([species.get_sequence(i) for i in indices]), offsets, lengths)
    return subset, indices


def fingerprint_result(result):
    arrays = {}
    for key in ("query_indices", "target_indices", "scores", "evalues"):
        array = np.ascontiguousarray(getattr(result, key))
        arrays[key] = {"shape": list(array.shape), "dtype": str(array.dtype),
                       "sha256": hashlib.sha256(array.tobytes()).hexdigest()}
    return {"query_species": result.query_species, "target_species": result.target_species,
            "candidate_count": int(result.candidate_count), "arrays": arrays}


def timed_wrapper(function, measurements, name):
    def wrapped(*args, **kwargs):
        wall, cpu = time.monotonic(), time.process_time()
        try:
            return function(*args, **kwargs)
        finally:
            row = measurements.setdefault(name, {"calls": 0, "wall_s": 0., "process_cpu_s": 0.})
            row["calls"] += 1
            row["wall_s"] += time.monotonic() - wall
            row["process_cpu_s"] += time.process_time() - cpu
    return wrapped


def search_pair(query, target, cap, threads, instrumented):
    measurements = {}
    wall, cpu = time.monotonic(), time.process_time()
    with ExitStack() as stack:
        if instrumented:
            for name in STAGES:
                original = getattr(engine, name)
                stack.enter_context(patch.object(engine, name, timed_wrapper(original, measurements, name)))
        result = engine.search_species_pair_indexed(query, target, "BLOSUM62", band_width=64,
            kmer_k=3, max_candidates_per_query=cap, n_threads=threads, prefilter_phase_a_topk=5000)
    return {"instrumented": instrumented, "wall_s": time.monotonic() - wall,
            "process_cpu_s": time.process_time() - cpu, "inclusive_stages": measurements,
            "fingerprint": fingerprint_result(result), "significant_hits": int(np.count_nonzero(result.evalues < 1e-4)),
            "raw_returned_hits": len(result.scores)}


def search_source_records(source, live_source):
    expected = {p.relative_to(live_source) for p in (live_source / "orthohmm/search").rglob("*.py")}
    actual = {p.relative_to(source) for p in (source / "orthohmm/search").rglob("*.py")}
    if actual != expected:
        raise ValueError("Search source inventory differs from the live frozen search")
    records = []
    for relative in sorted(expected):
        old, new = record(live_source / relative), record(source / relative)
        if old["sha256"] != new["sha256"]:
            raise ValueError("Search code differs from live frozen search: " + str(relative))
        records.extend((old, new))
    for name in ("hmm_viterbi.so", "kmer_prefilter.so", "pair_align.so"):
        relative = Path("orthohmm/search/csrc") / name
        old, new = record(live_source / relative), record(source / relative)
        if old["sha256"] != new["sha256"]:
            raise ValueError("CPU extension differs from live frozen search: " + name)
        records.extend((old, new))
    return records


def run(repo, source, output):
    if output.exists() or Path(__file__).resolve() != source / "benchmark_tools/probe_search_stage_cost.py":
        raise ValueError("Require untouched output and explicitly bound source")
    if os.environ.get("CUDA_VISIBLE_DEVICES") != "" or engine.is_cuda_available():
        raise ValueError("Require the prespecified CPU-only diagnostic")
    protocol_path = source / "benchmark_tools" / PROTOCOL
    protocol = json.loads(protocol_path.read_text())
    if protocol["settings"] != SETTINGS or protocol["cohort"]["queries_per_pair"] != 32:
        raise ValueError("Diagnostic settings differ from the prespecified sample/caps")
    live_path = repo / protocol["live_frozen_search_plan"]
    plan = json.loads(live_path.read_text())
    case = next(c for c in plan["cases"] if c["dataset"] == "QfO" and c["label"] == "control100")
    source_records = search_source_records(source, Path(plan["source_root"]))
    checked = [record(protocol_path), record(live_path), record(Path(__file__)), *source_records]
    pins = {Path(r["path"]).name: r for r in case["inputs"]}
    if len(pins) != len(case["inputs"]) or len(pins) != 78:
        raise ValueError("Require the complete frozen 78-species input inventory")
    filenames = selected_names(list(pins), 3)
    for name in filenames:
        check(pins[name])
        checked.append(pins[name])
    for pin in checked:
        check(pin)
    runtime = runtime_records()
    output.mkdir(parents=True)
    report = {"status": "running", "protocol": record(protocol_path), "live_search_plan": record(live_path),
        "source_revision": subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip(),
        "selected_species": filenames, "settings": protocol["settings"], "checked": checked,
        "runtime": runtime, "cases": [], "host_before": host_activity(), "truth_labels_loaded": False,
        "limitations": protocol["limitations"], "defaults_promoted": False, "goal_completed": False}
    write_json(output / "report.json", report)
    started = time.monotonic()
    try:
        loaded, load_costs = {}, {}
        for name in filenames:
            wall, cpu = time.monotonic(), time.process_time()
            loaded[name] = SpeciesSequences.from_fasta(pins[name]["path"], name)
            load_costs[name] = {"wall_s": time.monotonic() - wall, "process_cpu_s": time.process_time() - cpu,
                                "genes": loaded[name].num_sequences, "residues": int(loaded[name].lengths.sum())}
        report["fasta_load_costs"] = load_costs
        for i, name in enumerate(filenames):
            target_name = filenames[(i + 1) % len(filenames)]
            query, indices = subset_species(loaded[name], 32)
            target = loaded[target_name]
            for cap in (100, 400):
                row = {"query_species": name, "target_species": target_name, "candidate_cap": cap,
                    "selected_query_ids": query.ids, "original_query_indices": indices,
                    "query_lengths": query.lengths.tolist(), "target_genes": target.num_sequences, "runs": []}
                report["cases"].append(row)
                for instrumented in (False, True):
                    row["runs"].append(search_pair(query, target, cap, 4, instrumented))
                    write_json(output / "report.json", report)
                if row["runs"][0]["fingerprint"] != row["runs"][1]["fingerprint"]:
                    raise ValueError("Timing instrumentation changed native search outputs")
                row["instrumentation_byte_equivalent"] = True
                if time.monotonic() - started > 20 * 60:
                    raise TimeoutError("Prespecified diagnostic wall limit exceeded")
        for pin in [*checked, *runtime["files"]]:
            check(pin)
        if runtime_records() != runtime:
            raise ValueError("Scientific runtime changed during diagnostic")
        fields = dict(line.split(":", 1) for line in Path("/proc/self/status").read_text().splitlines() if ":" in line)
        peak = int(fields["VmHWM"].split()[0])
        if peak > 16 * 1024**2:
            raise MemoryError("Prespecified post-exec memory budget exceeded")
        report.update(status="complete_label_free_search_stage_probe_byte_equivalent",
            resources={"wall_s": time.monotonic() - started, "linux_postexec_peak_rss_kib": peak},
            inputs_reverified_unchanged=True)
    except Exception as error:
        report.update(status="failed_with_explicit_partial_diagnostic", error_type=type(error).__name__, error=str(error))
    report["host_after"] = host_activity()
    write_json(output / "report.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "source-root", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    result = run(args.repo.resolve(), args.source_root.resolve(), args.output.resolve())
    print(json.dumps({"status": result["status"]}))
    raise SystemExit(0 if result["status"].startswith("complete_") else 1)
