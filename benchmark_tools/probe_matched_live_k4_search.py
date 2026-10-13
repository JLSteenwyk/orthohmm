"""Confirm label-free stage costs using the untouched live k4 search source."""

import argparse
from contextlib import ExitStack
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from unittest.mock import patch


def search_pair(engine, helpers, query, target, cap, instrumented):
    costs = {}
    wall, cpu = time.monotonic(), time.process_time()
    with ExitStack() as stack:
        if instrumented:
            for name in helpers.STAGES:
                stack.enter_context(patch.object(engine, name,
                    helpers.timed_wrapper(getattr(engine, name), costs, name)))
        result = engine.search_species_pair_indexed(query, target, "BLOSUM62", band_width=64,
            kmer_k=4, max_candidates_per_query=cap, n_threads=4, prefilter_phase_a_topk=5000)
    import numpy as np
    return {"wall_s": time.monotonic() - wall, "process_cpu_s": time.process_time() - cpu,
        "inclusive_stages": costs, "fingerprint": helpers.fingerprint_result(result),
        "significant_hits": int(np.count_nonzero(result.evalues < 1e-4)), "instrumented": instrumented}


def run(repo, output):
    if output.exists() or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Require new CPU-only diagnostic namespace")
    live_path = repo / "benchmarks/results/full_search_hillclimb_r1_v3_20261012/execution_plan.json"
    live = json.loads(live_path.read_text())
    source = Path(live["source_root"])
    sys.path.insert(0, str(source))
    from orthohmm.search import engine
    from orthohmm.search.sequences import SpeciesSequences
    from orthohmm.accuracy import resolve_accuracy_profile
    if (Path(engine.__file__).resolve() != source / "orthohmm/search/engine.py"
            or resolve_accuracy_profile("high_sensitivity").kmer_k != 4 or engine.is_cuda_available()):
        raise ValueError("Wrong live source, actual high-sensitivity k or backend")
    helper_root = repo / "benchmarks/work/search_stage_probe_r1_v1_20261013"
    helper_path = helper_root / "benchmark_tools/probe_search_stage_cost.py"
    spec = importlib.util.spec_from_file_location("original_label_free_probe_helpers", helper_path)
    helpers = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helpers)
    case = next(c for c in live["cases"] if c["dataset"] == "QfO" and c["label"] == "control100")
    pins = {Path(pin["path"]).name: pin for pin in case["inputs"]}
    if len(pins) != 78 or len(pins) != len(case["inputs"]):
        raise ValueError("Incomplete frozen QfO input inventory")
    names = helpers.selected_names(list(pins), 3)
    protocol_path = Path(__file__).with_name("hillclimb_matched_k4_stage_probe_20261013.json")
    checked = [helpers.record(live_path), helpers.record(protocol_path), helpers.record(Path(__file__)),
               helpers.record(helper_path), *helpers.search_source_records(helper_root, source),
               helpers.record(source / "orthohmm/accuracy.py"),
               *[pins[name] for name in names]]
    for pin in checked:
        helpers.check(pin)
    runtime = helpers.runtime_records()
    output.mkdir(parents=True)
    report = {"status": "matched_k4_stage_probe_running", "checked": checked, "runtime": runtime,
        "job_id": os.environ.get("SLURM_JOB_ID"), "source_revision": live["source_revision"],
        "engine_source": helpers.record(engine.__file__), "actual_high_sensitivity_k": 4,
        "cases": [], "host_before": helpers.host_activity(), "truth_labels_loaded": False,
        "defaults_promoted": False, "goal_completed": False}
    helpers.write_json(output / "report.json", report)
    started = time.monotonic()
    try:
        loaded = {name: SpeciesSequences.from_fasta(pins[name]["path"], name) for name in names}
        for i, name in enumerate(names):
            query, indices = helpers.subset_species(loaded[name], 32)
            target = loaded[names[(i + 1) % 3]]
            for cap in (100, 400):
                row = {"query_species": name, "target_species": target.species_file, "candidate_cap": cap,
                    "selected_query_ids": query.ids, "original_query_indices": indices, "runs": []}
                report["cases"].append(row)
                for instrumented in (False, True):
                    row["runs"].append(search_pair(engine, helpers, query, target, cap, instrumented))
                    helpers.write_json(output / "report.json", report)
                if row["runs"][0]["fingerprint"] != row["runs"][1]["fingerprint"]:
                    raise ValueError("Matched-k4 instrumentation changed native search outputs")
                row["instrumentation_byte_equivalent"] = True
                if time.monotonic() - started > 20 * 60:
                    raise TimeoutError("Prespecified diagnostic wall budget exceeded")
        for pin in [*checked, *runtime["files"]]:
            helpers.check(pin)
        if helpers.runtime_records() != runtime:
            raise ValueError("Scientific runtime changed during matched diagnostic")
        fields = dict(line.split(":", 1) for line in Path("/proc/self/status").read_text().splitlines() if ":" in line)
        peak = int(fields["VmHWM"].split()[0])
        if peak > 16 * 1024**2:
            raise MemoryError("Prespecified post-exec memory budget exceeded")
        report.update(status="complete_matched_live_k4_stages_byte_equivalent", inputs_reverified_unchanged=True,
                      resources={"wall_s": time.monotonic() - started, "linux_postexec_peak_rss_kib": peak})
    except Exception as error:
        report.update(status="failed_matched_k4_stage_probe", error_type=type(error).__name__, error=str(error))
    report["host_after"] = helpers.host_activity()
    helpers.write_json(output / "report.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.repo.resolve(), args.output.resolve())
    print(result["status"])
    raise SystemExit(0 if result["status"] == "complete_matched_live_k4_stages_byte_equivalent" else 1)
