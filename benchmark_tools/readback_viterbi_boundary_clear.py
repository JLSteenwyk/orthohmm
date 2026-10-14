"""Rebuild matched native batches and independently admit boundary-clear timings."""

import argparse
import csv
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json

ARMS = ("original_full", "compiled_full", "boundary_only")


def timing_summary(case):
    repeats = case["repeats"]
    if len(repeats) != 6 or [tuple(row["order"]) for row in repeats] != list(itertools.permutations(ARMS)):
        raise ValueError("Incomplete or reordered boundary-clear timing schedule")
    for index, row in enumerate(repeats):
        if row["repeat"] != index or set(row["arms"]) != set(ARMS):
            raise ValueError("Invalid boundary-clear arm inventory")
        for arm in row["arms"].values():
            if any(not math.isfinite(arm[key]) or arm[key] <= 0 for key in ("wall_s", "process_cpu_s")):
                raise ValueError("Invalid measured boundary-clear cost")
    ratios = {}
    for baseline in ARMS[:2]:
        for key in ("wall_s", "process_cpu_s"):
            ratios[f"boundary_over_{baseline}_{key}"] = statistics.median(
                row["arms"]["boundary_only"][key] / row["arms"][baseline][key] for row in repeats)
    if ratios != case["ratios"]:
        raise ValueError("Producer ratios do not reproduce from raw repeats")
    return {"query": case["query"], "target": case["target"], "cap": case["cap"],
            "candidate_pairs": case["candidate_pairs"], "ratios": ratios}


def run(path, job_id, output):
    if output.exists() or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Require new CPU-only independent readback namespace")
    text = subprocess.check_output(["sacct", "-j", str(job_id), "-X",
        "--format=JobID,State,ExitCode,AllocCPUS", "--parsable2"], text=True)
    jobs = list(csv.DictReader(text.splitlines(), delimiter="|"))
    if len(jobs) != 1 or jobs[0] != dict(JobID=str(job_id), State="COMPLETED", ExitCode="0:0", AllocCPUS="4"):
        raise ValueError("Boundary-clear benchmark is not successfully terminal")
    raw_pin = record(path)
    report = json.loads(path.read_text())
    if (report["status"] != "complete_boundary_clear_exact_kernel_comparison" or len(report["cases"]) != 6
            or any(report[key] for key in ("truth_labels_loaded", "defaults_promoted", "goal_completed"))):
        raise ValueError("Incomplete boundary-clear benchmark or exceeded scientific scope")
    commands = report["compiler_commands"]
    full, boundary = commands["compiled_full"], commands["boundary_only"]
    if (boundary[:6] != full[:6] or boundary[6] != "-DHMM_VITERBI_BOUNDARY_CLEAR_ONLY=1"
            or boundary[7] != full[6] or len(full) != 9 or len(boundary) != 10):
        raise ValueError("Boundary-clear builds differ beyond the prespecified opt-in flag")
    source = Path(full[6]).parent.parent.parent.parent
    revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    if revision != report["source_revision"]:
        raise ValueError("Frozen benchmark source revision changed")
    checked = [*report["checked"], *report["runtime"]["files"], *report["kernels"].values(), report["compiler"]]
    for pin in checked:
        check(pin)
    sys.path.insert(0, str(source))
    import numpy as np
    from orthohmm.search import engine, viterbi
    from orthohmm.search.sequences import SpeciesSequences
    from benchmark_tools.readback_viterbi_center_recurrence import bind_library
    if Path(engine.__file__).resolve() != source / "orthohmm/search/engine.py":
        raise ValueError("Independent readback loaded the wrong frozen engine")
    from unittest.mock import patch
    libraries = {arm: bind_library(report["kernels"][arm]["path"]) for arm in ARMS}
    live_path = next(Path(pin["path"]) for pin in report["checked"] if Path(pin["path"]).name == "execution_plan.json")
    live = json.loads(live_path.read_text())
    inputs = next(c["inputs"] for c in live["cases"] if c["dataset"] == "QfO" and c["label"] == "control100")
    pins = {Path(pin["path"]).name: pin for pin in inputs}
    names = sorted(pins, key=lambda n: (hashlib.sha256(n.encode()).hexdigest(), n))[:3]
    expected = [(n, names[(i + 1) % 3], cap) for i, n in enumerate(names) for cap in (100, 400)]
    if [(c["query"], c["target"], c["cap"]) for c in report["cases"]] != expected:
        raise ValueError("Boundary-clear native batch sample changed")
    expected_settings = {"k": 4, "band_width": 64, "matrix": "BLOSUM62", "threads": 4,
        "use_reduced_alphabet": False, "min_total_hits": 4, "min_diag_hits": 1, "diag_bin_width": 10,
        "phase_a_topk": 5000, "caps": [100, 400], "center_mode": "division"}
    if report["settings"] != expected_settings:
        raise ValueError("Boundary-clear experiment changed matched k4 settings")
    loaded = {n: SpeciesSequences.from_fasta(pins[n]["path"], n) for n in names}
    matrix = engine.get_matrix("BLOSUM62")
    summaries = []
    for case in report["cases"]:
        full = loaded[case["query"]]
        chosen = set(sorted(full.ids, key=lambda n: (hashlib.sha256(n.encode()).hexdigest(), n))[:32])
        indices = [i for i, name in enumerate(full.ids) if name in chosen]
        ids = [full.ids[i] for i in indices]
        if ids != case["query_ids"] or indices != case["query_indices"]:
            raise ValueError("Boundary-clear query identities or original ordering changed")
        lengths = np.array(full.lengths[indices], dtype=np.int32)
        offsets = np.concatenate(([0], np.cumsum(lengths[:-1], dtype=np.int64)))
        query = SpeciesSequences(full.species_file, ids, np.concatenate([full.get_sequence(i) for i in indices]), offsets, lengths)
        target = loaded[case["target"]]
        profiles = engine.build_profiles_batch(query, matrix, engine.get_background_freqs("BLOSUM62"))
        candidates, boundaries = engine.prefilter_candidates(query, target, k=4, use_reduced_alphabet=False,
            min_total_hits=4, min_diag_hits=1, diag_bin_width=10, max_candidates_per_query=case["cap"],
            n_threads=4, sub_matrix=matrix, phase_a_topk=5000)
        pairs = np.column_stack((np.repeat(np.arange(32, dtype=np.int32), np.diff(boundaries)), candidates))
        data = (*profiles, target.flat_sequences, target.offsets, target.lengths, pairs)
        def fingerprints():
            return [{"shape": list(a.shape), "dtype": str(a.dtype), "sha256": hashlib.sha256(a.tobytes()).hexdigest()} for a in data]
        if fingerprints() != case["input_arrays"] or len(pairs) != case["candidate_pairs"]:
            raise ValueError("Reconstructed boundary-clear input arrays changed")
        for lib in libraries.values():
            with patch.object(viterbi, "_c_hmm_lib", lib), patch.object(viterbi, "_c_hmm_available", True):
                actual = viterbi.batch_viterbi_multipair_c(*data, 64, n_threads=4)
            if actual.dtype != np.int32 or actual.shape != (len(pairs),) or hashlib.sha256(actual.tobytes()).hexdigest() != case["score_sha256"]:
                raise ValueError("Reconstructed boundary-clear integer scores changed")
        if fingerprints() != case["input_arrays"]:
            raise ValueError("Boundary-clear replay modified input arrays")
        summaries.append(timing_summary(case))
    for pin in [raw_pin, *checked]:
        check(pin)
    result = {"status": "terminal_boundary_clear_batches_scores_and_timings_independently_verified",
        "raw_report": raw_pin, "scheduler": jobs[0], "source_revision": revision, "reader": record(Path(__file__)),
        "cases": summaries, "scope": "Exact matched-k4 kernel replay and timing arithmetic; no eight-endpoint accuracy or isolated full-pipeline speed claim.",
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
