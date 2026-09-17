"""Compose fresh preparation, runtime/input verification and native measurement.

This library does not authorize or submit runs. The caller must validate a frozen
execution authorization and supply the pinned collector and native enumerator.
"""

import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_verified_slurm_measurement import check_manifests, run_checked
from benchmark_tools.prepare_native_scaling_run import check_copies, check_original_inputs, prepare, run_directory


def measure_run(run, order, runtime_specs, enumerate_native, collector, job_id,
                cpus=20, memory_gib=96, timeout_s=86400, runtime_checker=check_manifests):
    root = run_directory(run)
    if Path.cwd().resolve() != Path(run["cwd"]).resolve():
        raise ValueError("Caller working directory differs from frozen native cwd")
    prepared = None

    def checker(specifications):
        runtime = runtime_checker(specifications)
        originals = check_original_inputs(run, order)
        observed_order = enumerate_native(run["dataset"]["input_directory"])
        if observed_order != order["native_order"]:
            raise ValueError("Actual frozen OrthoHMM input enumeration changed")
        result = {"runtime": runtime, "original_inputs": originals, "native_order": observed_order}
        if prepared is not None and run["native_method"] == "orthofinder_full":
            result["copied_inputs"] = check_copies(run)
        return result

    def measurement(directory):
        nonlocal prepared
        if str(directory) != run["measurement_directory"]:
            raise ValueError("Collector output differs from frozen run")
        start = time.monotonic()
        prepared = prepare(run, order)
        # Input copying can take time; recheck native enumeration at the handoff.
        if enumerate_native(run["dataset"]["input_directory"]) != order["native_order"]:
            raise ValueError("Native input enumeration changed during preparation")
        prepared["preparation_wall_s"] = time.monotonic() - start
        with (root / "preparation.json").open("x") as handle:
            json.dump(prepared, handle, indent=2, sort_keys=True)
            handle.write("\n")
        return collector(prepared["measured_argv"], directory, job_id, cpus,
                         memory_gib * 1024 ** 3, timeout_s, 1., monitor_host=True, host_interval_s=30.)

    return run_checked(runtime_specs, root, measurement, checker)
