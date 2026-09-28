"""Compose local preparation and measurement; never authorize or submit a run.

Caller must bind the run, baseline, runtime manifests and collector to a frozen
recipe and verify allocation and environmental eligibility before invocation.
There is deliberately no execution CLI.
"""

import os
from pathlib import Path
import shutil

from benchmark_tools.measure_native_scaling_run import run_checked, check_manifests
from benchmark_tools.prepare_threadripper_run import paths, prepare, check_prepared
from benchmark_tools.run_simulation_methods import verify_environment
from benchmark_tools.prepare_ob_candidate_neighborhood import check
from benchmark_tools.isolated_numba_cache import fresh_cache
from benchmark_tools.isolated_native_tmp import fresh_tmp
from benchmark_tools.snapshot_orthohmm_input_order import record
from benchmark_tools.probe_dgx_step_separation import save


def assert_environment(run, baseline):
    if Path.cwd().resolve() != Path(run["cwd"]).resolve():
        raise ValueError("Caller cwd differs from native run")
    for key, expected in baseline["environment_overrides"].items():
        if os.environ.get(key) != expected:
            raise ValueError("Caller environment differs: " + key)
    for name in ("mafft", "FastTree", "diamond"):
        actual = shutil.which(name)
        expected = baseline["tool_entrypoints"][name]["absolute_path"]
        if actual is None or Path(actual).resolve() != Path(expected).resolve():
            raise ValueError("Actual native PATH resolution differs: " + name)
    for key in ("LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"):
        if os.environ.get(key):
            raise ValueError("Unexpected loader override: " + key)
    if Path("/etc/ld.so.preload").exists():
        raise ValueError("System preload configuration needs separate review")


def measure_run(run, baseline, runtime_specs, collector, job_id, runtime_checker=check_manifests):
    if type(job_id) is not int or job_id <= 0 or not runtime_specs or not callable(collector):
        raise ValueError("Require explicit job identity, runtime pins and collector")
    root, target = paths(run)
    assert_environment(run, baseline)
    checks = 0

    def checker(specifications):
        nonlocal checks
        checks += 1
        assert_environment(run, baseline)
        runtime = runtime_checker(specifications)
        verification_cache = root.parent / (root.name + "_verification_caches") / str(checks)
        with fresh_cache(verification_cache):
            verify_environment(baseline)
        for row in run["dataset"]["inputs"]:
            check(row)
        result = dict(runtime=runtime, original_inputs=run["dataset"]["inputs"],
                      verification_numba_cache=str(verification_cache))
        if target.exists():
            result["prepared_inputs"] = check_prepared(run, baseline)
        return result

    def measurement(directory):
        if str(directory) != run["measurement_directory"]:
            raise ValueError("Unexpected collector output directory")
        prepare(run, baseline)
        # Copying/hashing is outside native timing; recheck immediately at handoff.
        check_prepared(run, baseline)
        assert_environment(run, baseline)
        cache = root / "native_numba_cache"
        with fresh_tmp(root / "native_tmp"), fresh_cache(cache):
            try:
                return collector(run["native_argv"], directory, job_id, 32, 128*1024**3,
                                 85800, 1., monitor_host=True, host_interval_s=30.)
            finally:
                save(root / "numba_cache.json", dict(directory=str(cache), initially_empty=True,
                    files=[record(p) for p in sorted(cache.rglob("*")) if p.is_file()],
                    native_tmp=str(root / "native_tmp"),
                    policy="fresh_native_cache_compilation_inside_native_timer",
                    limitations=["Cache file inventory is after inference; not a trace of every cache lookup.",
                                 "Verification probes use different fresh caches; historical caches are untouched."]))

    return run_checked(runtime_specs, root, measurement, checker)
