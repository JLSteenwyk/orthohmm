"""Compare opt-in boundary clearing against compiled and frozen full-clear controls."""

import argparse
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
from benchmark_tools.benchmark_viterbi_center_recurrence import array_records, bind_library, score
from benchmark_tools.probe_search_stage_cost import selected_names, subset_species
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.run_search_recovery import host_activity, runtime_records
from orthohmm.search import engine
from orthohmm.search.sequences import SpeciesSequences

PROTOCOL = "hillclimb_viterbi_boundary_clear_20261013.json"
ARMS = ("original_full", "compiled_full", "boundary_only")


def paired_times(data, libraries):
    if tuple(libraries) != ARMS:
        raise ValueError("Require the three prespecified boundary-clear arms")
    before = array_records(data)
    reference = score(data, libraries["original_full"], "division")
    for arm in ARMS[1:]:
        if not np.array_equal(score(data, libraries[arm], "division"), reference):
            raise ValueError("Warm boundary-clear scores differ from frozen control")
    repeats = []
    for index, order in enumerate(itertools.permutations(ARMS)):
        row = {"repeat": index, "order": order, "arms": {}}
        for arm in order:
            wall, cpu = time.monotonic(), time.process_time()
            actual = score(data, libraries[arm], "division")
            row["arms"][arm] = {"wall_s": time.monotonic() - wall,
                                 "process_cpu_s": time.process_time() - cpu}
            if not np.array_equal(actual, reference):
                raise ValueError("Boundary-clear optimization changed integer scores")
        repeats.append(row)
    if array_records(data) != before:
        raise ValueError("Kernel modified benchmark input arrays")
    ratios = {}
    for baseline in ARMS[:2]:
        for field in ("wall_s", "process_cpu_s"):
            values = [r["arms"]["boundary_only"][field] / r["arms"][baseline][field] for r in repeats]
            ratios[f"boundary_over_{baseline}_{field}"] = float(np.median(values))
    return {"repeats": repeats, "ratios": ratios, "inputs_unchanged": True,
            "exact_scores_match_original": True, "score_sha256": hashlib.sha256(reference.tobytes()).hexdigest()}


def run(repo, source, output):
    if (output.exists() or Path(__file__).resolve() != source / "benchmark_tools/benchmark_viterbi_boundary_clear.py"
            or os.environ.get("CUDA_VISIBLE_DEVICES") != ""):
        raise ValueError("Require new CPU-only namespace and explicitly bound frozen source")
    live_path = repo / "benchmarks/results/full_search_hillclimb_r1_v3_20261012/execution_plan.json"
    live = json.loads(live_path.read_text())
    case = next(c for c in live["cases"] if c["dataset"] == "QfO" and c["label"] == "control100")
    pins = {Path(pin["path"]).name: pin for pin in case["inputs"]}
    if len(pins) != 78 or len(pins) != len(case["inputs"]):
        raise ValueError("Incomplete frozen QfO input inventory")
    names = selected_names(list(pins), 3)
    original_library = Path(live["source_root"]) / "orthohmm/search/csrc/hmm_viterbi.so"
    checked = [record(live_path), record(original_library), record(source / "benchmark_tools" / PROTOCOL),
               record(Path(__file__)), record(source / "orthohmm/search/csrc/hmm_viterbi.c"), *[pins[n] for n in names]]
    checked.extend(record(p) for p in sorted((source / "orthohmm/search").rglob("*.py")))
    checked.extend(record(source / "orthohmm/search/csrc" / name) for name in ("kmer_prefilter.so", "pair_align.so"))
    checked.extend(record(source / "benchmark_tools" / name) for name in (
        "benchmark_viterbi_center_recurrence.py", "probe_search_stage_cost.py",
        "prepare_ob_candidate_neighborhood.py", "run_membership_hillclimb.py", "run_search_recovery.py"))
    for pin in checked:
        check(pin)
    subprocess.run(["git", "-C", str(source), "diff", "--exit-code", "HEAD", "--",
                    "orthohmm", "benchmark_tools"], check=True, capture_output=True)
    output.mkdir(parents=True)
    report = {"status": "boundary_clear_benchmark_running", "checked": checked, "cases": [],
        "source_revision": subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip(),
        "runtime": runtime_records(), "host_before": host_activity(), "truth_labels_loaded": False,
        "settings": {"k": 4, "band_width": 64, "matrix": "BLOSUM62", "threads": 4,
            "use_reduced_alphabet": False, "min_total_hits": 4, "min_diag_hits": 1,
            "diag_bin_width": 10, "phase_a_topk": 5000, "caps": [100, 400], "center_mode": "division"},
        "defaults_promoted": False, "goal_completed": False}
    write_json(output / "report.json", report)
    started = time.monotonic()
    try:
        compiler = shutil.which("gcc")
        if compiler is None:
            raise FileNotFoundError("GCC unavailable")
        commands, kernel_pins = {}, {"original_full": record(original_library)}
        libraries = {"original_full": bind_library(original_library)}
        for arm in ARMS[1:]:
            path = output / (arm + ".so")
            flag = ["-DHMM_VITERBI_BOUNDARY_CLEAR_ONLY=1"] if arm == "boundary_only" else []
            command = [compiler, "-O3", "-fopenmp", "-shared", "-fPIC", "-march=native", *flag,
                       str(source / "orthohmm/search/csrc/hmm_viterbi.c"), "-o", str(path)]
            subprocess.run(command, check=True, capture_output=True)
            commands[arm], kernel_pins[arm] = command, record(path)
            libraries[arm] = bind_library(path)
        report.update(compiler=record(compiler), compiler_commands=commands, kernels=kernel_pins)
        loaded = {n: SpeciesSequences.from_fasta(pins[n]["path"], n) for n in names}
        matrix = engine.get_matrix("BLOSUM62")
        for i, name in enumerate(names):
            query, indices = subset_species(loaded[name], 32)
            target = loaded[names[(i + 1) % 3]]
            profiles = engine.build_profiles_batch(query, matrix, engine.get_background_freqs("BLOSUM62"))
            for cap in (100, 400):
                candidates, offsets = engine.prefilter_candidates(query, target, k=4, use_reduced_alphabet=False,
                    min_total_hits=4, min_diag_hits=1, diag_bin_width=10, max_candidates_per_query=cap,
                    n_threads=4, sub_matrix=matrix, phase_a_topk=5000)
                pairs = np.column_stack((np.repeat(np.arange(32, dtype=np.int32), np.diff(offsets)), candidates))
                data = (*profiles, target.flat_sequences, target.offsets, target.lengths, pairs)
                fingerprints = array_records(data)
                result = paired_times(data, libraries)
                report["cases"].append({"query": name, "target": target.species_file, "cap": cap,
                    "query_ids": query.ids, "query_indices": indices, "candidate_pairs": len(pairs),
                    "input_arrays": fingerprints, **result})
                write_json(output / "report.json", report)
                if time.monotonic() - started > 20 * 60:
                    raise TimeoutError("Prespecified boundary benchmark wall budget exceeded")
        for pin in [*checked, *report["runtime"]["files"], *kernel_pins.values(), report["compiler"]]:
            check(pin)
        if runtime_records() != report["runtime"]:
            raise ValueError("Scientific runtime changed during boundary benchmark")
        fields = dict(line.split(":", 1) for line in Path("/proc/self/status").read_text().splitlines() if ":" in line)
        peak = int(fields["VmHWM"].split()[0])
        if peak > 16 * 1024**2:
            raise MemoryError("Prespecified post-exec memory budget exceeded")
        report.update(status="complete_boundary_clear_exact_kernel_comparison", inputs_reverified_unchanged=True,
            resources={"wall_s": time.monotonic() - started, "linux_postexec_peak_rss_kib": peak})
    except Exception as error:
        report.update(status="failed_boundary_clear_kernel_comparison", error_type=type(error).__name__, error=str(error))
    report["host_after"] = host_activity()
    write_json(output / "report.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "source-root", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = run(args.repo.resolve(), args.source_root.resolve(), args.output.resolve())
    print(result["status"])
    raise SystemExit(0 if result["status"] == "complete_boundary_clear_exact_kernel_comparison" else 1)
