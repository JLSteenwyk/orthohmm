"""Reconstruct kernel batches and independently summarize terminal timings."""

import argparse
import csv
import ctypes
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json


def bind_library(path):
    lib = ctypes.CDLL(str(path))
    lib.hmm_set_num_threads.argtypes = [ctypes.c_int32]
    lib.hmm_set_num_threads.restype = None
    for name in ("batch_hmm_viterbi_multipair_avx2_c", "batch_hmm_viterbi_multipair_recurrence_c"):
        if hasattr(lib, name):
            function = getattr(lib, name)
            function.argtypes = [ctypes.c_void_p] * 9 + [ctypes.c_int32, ctypes.c_int32, ctypes.c_void_p]
            function.restype = None
    return lib


def timing_summary(case):
    rows = case["repeats"]
    arms = ("original_division", "division", "recurrence")
    if len(rows) != 6 or [tuple(row["order"]) for row in rows] != list(itertools.permutations(arms)):
        raise ValueError("Incomplete or reordered three-arm timing schedule")
    for index, row in enumerate(rows):
        if row["repeat"] != index or set(row["modes"]) != set(arms):
            raise ValueError("Invalid timing arm inventory")
        for value in row["modes"].values():
            for field in ("wall_s", "process_cpu_s"):
                if not math.isfinite(value[field]) or value[field] <= 0:
                    raise ValueError("Invalid measured time")
    result = {"query": case["query"], "target": case["target"], "cap": case["cap"],
              "candidate_pairs": case["candidate_pairs"]}
    for baseline in ("division", "original_division"):
        for field in ("wall_s", "process_cpu_s"):
            ratios = [r["modes"]["recurrence"][field] / r["modes"][baseline][field] for r in rows]
            result[f"median_recurrence_over_{baseline}_{field}"] = statistics.median(ratios)
    if (result["median_recurrence_over_division_wall_s"] != case["median_recurrence_over_division_wall"]
            or result["median_recurrence_over_original_division_wall_s"] != case["median_recurrence_over_original_wall"]):
        raise ValueError("Producer median does not match raw repeat measurements")
    return result


def run(path, job_id, output):
    if output.exists() or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Require fresh readback namespace and CPU-only runtime")
    scheduler_text = subprocess.check_output(["sacct", "-j", str(job_id), "-X",
        "--format=JobID,State,ExitCode", "--parsable2"], text=True)
    scheduler = list(csv.DictReader(scheduler_text.splitlines(), delimiter="|"))
    if len(scheduler) != 1 or scheduler[0]["JobID"] != str(job_id) or scheduler[0]["State"] != "COMPLETED" or scheduler[0]["ExitCode"] != "0:0":
        raise ValueError("Kernel comparison job is not successfully terminal")
    raw_pin = record(path)
    report = json.loads(path.read_text())
    if report["status"] != "complete_exact_kernel_comparison_no_default_promotion" or len(report["cases"]) != 6:
        raise ValueError("Incomplete kernel comparison report")
    if any(report[field] for field in ("truth_labels_loaded", "defaults_promoted", "goal_completed")):
        raise ValueError("Kernel report exceeds its scientific scope")
    for pin in [*report["checked"], *report["runtime"]["files"], report["kernel"], report["compiler"]]:
        check(pin)
    source = Path(report["compiler_command"][6]).parent.parent.parent.parent
    revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    if revision != report["source_revision"]:
        raise ValueError("Frozen source revision changed")
    sys.path.insert(0, str(source))
    from orthohmm.search import engine, viterbi
    from orthohmm.search.sequences import SpeciesSequences
    if Path(engine.__file__).resolve() != source / "orthohmm/search/engine.py":
        raise ValueError("Readback imported the wrong inference source")
    live_path = next(Path(p["path"]) for p in report["checked"] if Path(p["path"]).name == "execution_plan.json")
    live = json.loads(live_path.read_text())
    inputs = next(c["inputs"] for c in live["cases"] if c["dataset"] == "QfO" and c["label"] == "control100")
    pins = {Path(p["path"]).name: p for p in inputs}
    names = sorted(pins, key=lambda n: (hashlib.sha256(n.encode()).hexdigest(), n))[:3]
    expected = [(n, names[(i + 1) % 3], cap) for i, n in enumerate(names) for cap in (100, 400)]
    if [(c["query"], c["target"], c["cap"]) for c in report["cases"]] != expected:
        raise ValueError("Kernel sample differs from prespecified hash selection")
    settings = {"k": 3, "use_reduced_alphabet": False, "min_total_hits": 4, "min_diag_hits": 1,
                "diag_bin_width": 10, "phase_a_topk": 5000, "band_width": 64,
                "threads": 4, "caps": [100, 400], "matrix": "BLOSUM62"}
    if report["settings"] != settings:
        raise ValueError("Kernel settings differ from prespecified pipeline settings")
    current = bind_library(report["kernel"]["path"])
    original = bind_library(Path(live["source_root"]) / "orthohmm/search/csrc/hmm_viterbi.so")
    loaded = {n: SpeciesSequences.from_fasta(pins[n]["path"], n) for n in names}
    matrix = engine.get_matrix("BLOSUM62")
    summaries = []
    for case in report["cases"]:
        full = loaded[case["query"]]
        chosen = set(sorted(full.ids, key=lambda n: (hashlib.sha256(n.encode()).hexdigest(), n))[:32])
        indices = [i for i, name in enumerate(full.ids) if name in chosen]
        ids = [full.ids[i] for i in indices]
        if ids != case["query_ids"] or indices != case["query_indices"]:
            raise ValueError("Query identities or original ordering changed")
        lengths = np.array(full.lengths[indices], dtype=np.int32)
        offsets = np.concatenate(([0], np.cumsum(lengths[:-1], dtype=np.int64)))
        query = SpeciesSequences(full.species_file, ids, np.concatenate([full.get_sequence(i) for i in indices]), offsets, lengths)
        target = loaded[case["target"]]
        profiles = engine.build_profiles_batch(query, matrix, engine.get_background_freqs("BLOSUM62"))
        candidates, boundaries = engine.prefilter_candidates(query, target, k=3,
            use_reduced_alphabet=False, min_total_hits=4, min_diag_hits=1, diag_bin_width=10,
            max_candidates_per_query=case["cap"], n_threads=4, sub_matrix=matrix, phase_a_topk=5000)
        pairs = np.column_stack((np.repeat(np.arange(32, dtype=np.int32), np.diff(boundaries)), candidates))
        data = (*profiles, target.flat_sequences, target.offsets, target.lengths, pairs)
        fingerprints = [{"shape": list(a.shape), "dtype": str(a.dtype), "sha256": hashlib.sha256(a.tobytes()).hexdigest()} for a in data]
        if len(pairs) != case["candidate_pairs"] or fingerprints != case["input_arrays"]:
            raise ValueError("Reconstructed native batch changed")
        for lib, mode in ((original, "division"), (current, "division"), (current, "recurrence")):
            with patch.object(viterbi, "_c_hmm_lib", lib), patch.object(viterbi, "_c_hmm_available", True):
                actual = viterbi.batch_viterbi_multipair_c(*data, 64, n_threads=4, band_center_mode=mode)
            if actual.dtype != np.int32 or actual.shape != (len(pairs),) or hashlib.sha256(actual.tobytes()).hexdigest() != case["score_sha256"]:
                raise ValueError("Reconstructed integer scores changed")
        after = [{"shape": list(a.shape), "dtype": str(a.dtype), "sha256": hashlib.sha256(a.tobytes()).hexdigest()} for a in data]
        if after != fingerprints:
            raise ValueError("Native replay modified reconstructed input arrays")
        summaries.append(timing_summary(case))
    for pin in [*report["checked"], *report["runtime"]["files"], report["kernel"], report["compiler"]]:
        check(pin)
    check(raw_pin)
    result = {"status": "terminal_kernel_arrays_scores_and_timings_independently_verified", "raw_report": raw_pin,
        "scheduler": scheduler[0], "source_revision": revision, "cases": summaries,
        "reader": record(Path(__file__)), "python_executable": record(sys.executable),
        "scope": "Exact kernel-output replay and raw timing arithmetic; no new accuracy or full-pipeline speed claim.",
        "defaults_promoted": False, "goal_completed": False}
    output.mkdir(parents=True)
    write_json(output / "report.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--job-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(run(args.report.resolve(), args.job_id, args.output.resolve())["status"])
