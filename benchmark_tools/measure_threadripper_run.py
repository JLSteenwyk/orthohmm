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

    def checker(specifications):
        assert_environment(run, baseline)
        runtime = runtime_checker(specifications)
        verify_environment(baseline)
        for row in run["dataset"]["inputs"]:
            check(row)
        result = dict(runtime=runtime, original_inputs=run["dataset"]["inputs"])
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
        return collector(run["native_argv"], directory, job_id, 32, 128*1024**3,
                         85800, 1., monitor_host=True, host_interval_s=30.)

    return run_checked(runtime_specs, root, measurement, checker)
