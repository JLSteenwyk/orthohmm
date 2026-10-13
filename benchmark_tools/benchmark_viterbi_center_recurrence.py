"""Benchmark opt-in exact band centers on fixed label-free native batches."""

import argparse
import ctypes
import hashlib
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from orthohmm.search import engine, viterbi
from orthohmm.search.sequences import SpeciesSequences
from benchmark_tools.probe_search_stage_cost import selected_names, subset_species
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.run_search_recovery import host_activity, runtime_records

PROTOCOL = "hillclimb_viterbi_center_recurrence_20261013.json"


def array_records(data):
    return [{"shape": list(a.shape), "dtype": str(a.dtype),
             "sha256": hashlib.sha256(a.tobytes()).hexdigest()} for a in data]


def bind_library(path):
    lib = ctypes.CDLL(str(path))
    lib.hmm_have_avx2.restype = ctypes.c_int32
    if not lib.hmm_have_avx2():
        raise RuntimeError("Kernel experiment requires native AVX2")
    lib.hmm_set_num_threads.argtypes = [ctypes.c_int32]
    lib.hmm_set_num_threads.restype = None
    for name in ("batch_hmm_viterbi_multipair_avx2_c", "batch_hmm_viterbi_multipair_recurrence_c"):
        if hasattr(lib, name):
            fn = getattr(lib, name)
            fn.argtypes = [ctypes.c_void_p] * 9 + [ctypes.c_int32, ctypes.c_int32, ctypes.c_void_p]
            fn.restype = None
    return lib


def score(data, lib, mode):
    previous = viterbi._c_hmm_lib, viterbi._c_hmm_available
    viterbi._c_hmm_lib, viterbi._c_hmm_available = lib, True
    try:
        return viterbi.batch_viterbi_multipair_c(*data, 64, n_threads=4, band_center_mode=mode)
    finally:
        viterbi._c_hmm_lib, viterbi._c_hmm_available = previous


def paired_times(data, current, original):
    before = [a.tobytes() for a in data]
    reference = score(data, original, "division")
    for mode in ("division", "recurrence"):
        if not np.array_equal(score(data, current, mode), reference):
            raise ValueError("Warm kernel scores differ from the original frozen kernel")
    measurements = []
    arms = ("original_division", "division", "recurrence")
    for repeat, order in enumerate(itertools.permutations(arms)):
        row = {"repeat": repeat, "order": order, "modes": {}}
        for arm in order:
            lib, mode = (original, "division") if arm == "original_division" else (current, arm)
            wall, cpu = time.monotonic(), time.process_time()
            actual = score(data, lib, mode)
            row["modes"][arm] = {"wall_s": time.monotonic() - wall,
                                  "process_cpu_s": time.process_time() - cpu}
            if not np.array_equal(actual, reference):
                raise ValueError("Band-center optimization changed integer scores")
        measurements.append(row)
    if [a.tobytes() for a in data] != before:
        raise ValueError("Kernel modified input arrays")
    ratios = [r["modes"]["recurrence"]["wall_s"] / r["modes"]["division"]["wall_s"] for r in measurements]
    original_ratios = [r["modes"]["recurrence"]["wall_s"] / r["modes"]["original_division"]["wall_s"] for r in measurements]
    return {"repeats": measurements, "median_recurrence_over_division_wall": float(np.median(ratios)),
            "median_recurrence_over_original_wall": float(np.median(original_ratios)),
            "exact_scores_match_live_kernel": True, "inputs_unchanged": True,
            "score_sha256": hashlib.sha256(reference.tobytes()).hexdigest()}


def run(repo, source, output):
    if output.exists() or Path(__file__).resolve() != source / "benchmark_tools/benchmark_viterbi_center_recurrence.py":
        raise ValueError("Require a fresh namespace and bound frozen source")
    if os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Require CPU-only kernel measurement")
    live_path = repo / "benchmarks/results/full_search_hillclimb_r1_v3_20261012/execution_plan.json"
    live = json.loads(live_path.read_text())
    case = next(c for c in live["cases"] if c["dataset"] == "QfO" and c["label"] == "control100")
    pins = {Path(p["path"]).name: p for p in case["inputs"]}
    if len(pins) != 78 or len(pins) != len(case["inputs"]):
        raise ValueError("Require the frozen complete 78-species input inventory")
    names = selected_names(list(pins), 3)
    old_library = Path(live["source_root"]) / "orthohmm/search/csrc/hmm_viterbi.so"
    checked = [record(live_path), record(old_library), record(source / "benchmark_tools" / PROTOCOL),
               record(source / "benchmark_tools/hillclimb_viterbi_center_benchmark_addendum_20261013.json"),
               record(source / "orthohmm/search/csrc/hmm_viterbi.c"), record(source / "orthohmm/search/viterbi.py"),
               record(Path(__file__)), *[pins[n] for n in names]]
    checked.extend(record(p) for p in sorted((source / "orthohmm/search").rglob("*.py")))
    checked.extend(record(source / "orthohmm/search/csrc" / name)
                   for name in ("kmer_prefilter.so", "pair_align.so"))
    checked.extend(record(source / "benchmark_tools" / name) for name in (
        "probe_search_stage_cost.py", "prepare_ob_candidate_neighborhood.py",
        "run_membership_hillclimb.py", "run_search_recovery.py"))
    for pin in checked:
        check(pin)
    subprocess.run(["git", "-C", str(source), "diff", "--exit-code", "HEAD", "--",
                    "orthohmm", "benchmark_tools"], check=True, capture_output=True)
    output.mkdir(parents=True)
    report = {"status": "running", "checked": checked, "source_revision": subprocess.check_output(
            ["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip(), "cases": [],
        "runtime": runtime_records(), "host_before": host_activity(), "truth_labels_loaded": False,
        "settings": {"k": 3, "use_reduced_alphabet": False, "min_total_hits": 4,
            "min_diag_hits": 1, "diag_bin_width": 10, "phase_a_topk": 5000,
            "band_width": 64, "threads": 4, "caps": [100, 400], "matrix": "BLOSUM62"},
        "defaults_promoted": False, "goal_completed": False}
    write_json(output / "report.json", report)
    started = time.monotonic()
    try:
        compiler = shutil.which("gcc")
        if compiler is None:
            raise FileNotFoundError("GCC unavailable")
        library = output / "hmm_viterbi.so"
        command = [compiler, "-O3", "-fopenmp", "-shared", "-fPIC", "-march=native",
                   str(source / "orthohmm/search/csrc/hmm_viterbi.c"), "-o", str(library)]
        subprocess.run(command, check=True, capture_output=True)
        current, original = bind_library(library), bind_library(old_library)
        report.update(compiler_command=command, compiler=record(compiler), kernel=record(library))
        loaded = {n: SpeciesSequences.from_fasta(pins[n]["path"], n) for n in names}
        matrix = engine.get_matrix("BLOSUM62")
        for i, name in enumerate(names):
            query, indices = subset_species(loaded[name], 32)
            target = loaded[names[(i + 1) % 3]]
            profiles = engine.build_profiles_batch(query, matrix, engine.get_background_freqs("BLOSUM62"))
            for cap in (100, 400):
                candidates, offsets = engine.prefilter_candidates(query, target, k=3, max_candidates_per_query=cap,
                    use_reduced_alphabet=False, min_total_hits=4, min_diag_hits=1, diag_bin_width=10,
                    n_threads=4, sub_matrix=matrix, phase_a_topk=5000)
                pairs = np.column_stack((np.repeat(np.arange(32, dtype=np.int32), np.diff(offsets)), candidates))
                data = (*profiles, target.flat_sequences, target.offsets, target.lengths, pairs)
                input_arrays = array_records(data)
                result = paired_times(data, current, original)
                report["cases"].append({"query":name,"target":target.species_file,"cap":cap,
                    "query_ids":query.ids,"query_indices":indices,"candidate_pairs":len(pairs),
                    "input_arrays": input_arrays, **result})
                write_json(output / "report.json", report)
                if time.monotonic() - started > 20 * 60:
                    raise TimeoutError("Prespecified kernel benchmark wall limit exceeded")
        for pin in [*checked, *report["runtime"]["files"]]:
            check(pin)
        if runtime_records() != report["runtime"]:
            raise ValueError("Scientific runtime changed during kernel benchmark")
        subprocess.run(["git", "-C", str(source), "diff", "--exit-code", "HEAD", "--",
                        "orthohmm", "benchmark_tools"], check=True, capture_output=True)
        fields = dict(line.split(":", 1) for line in Path("/proc/self/status").read_text().splitlines() if ":" in line)
        peak = int(fields["VmHWM"].split()[0])
        if peak > 16 * 1024**2:
            raise MemoryError("Prespecified post-exec memory budget exceeded")
        report["resources"] = {"wall_s": time.monotonic() - started, "linux_postexec_peak_rss_kib": peak}
        report["status"] = "complete_exact_kernel_comparison_no_default_promotion"
    except Exception as error:
        report.update(status="failed_with_partial_kernel_comparison", error_type=type(error).__name__, error=str(error))
    report["host_after"] = host_activity()
    write_json(output / "report.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "source-root", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = run(args.repo.resolve(), args.source_root.resolve(), args.output.resolve())
    print(json.dumps({"status": result["status"]}))
    raise SystemExit(0 if result["status"].startswith("complete_") else 1)
