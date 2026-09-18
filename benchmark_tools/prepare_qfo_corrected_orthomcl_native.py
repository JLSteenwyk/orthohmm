"""Prepare isolated corrected OrthoMCL native sources without launching inference."""

import argparse
import json
from pathlib import Path
import re
import shutil
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.configure_orthomcl_1_4 import configure_module
from benchmark_tools.parallelize_orthomcl_pairs import parallelize_source
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_orthomcl_bpo_parity import TOOL
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.snapshot_runtime_trees import verify

RUNTIME_SHA = "9cc29e84777f47064c23096ce36e32d5465616979b05febe510219d58b6b3239"
HELPERS_SHA = "79ef5d6297914cf43e318858de4de42e58f31d601159a9e27bb8515734efaa2c"


def write_sources(destination, data, threads, parallel=True):
    if type(threads) is not int or threads < 1 or type(parallel) is not bool:
        raise ValueError("Invalid native resource/configuration settings")
    if any(not path.is_absolute() or not re.fullmatch(r"[A-Za-z0-9_./-]+", str(path))
           for path in (destination, data)):
        raise ValueError("Require absolute shell-safe native paths")
    if destination.exists():
        raise FileExistsError(destination)
    originals = [record(TOOL / name) for name in ("orthomcl.pl", "orthomcl_module.pm")]
    destination.mkdir(parents=True, exist_ok=False)
    for name in ("orthomcl.pl", "orthomcl_module.pm"):
        shutil.copyfile(TOOL / name, destination / name)
    configure_module(destination / "orthomcl_module.pm", destination, data, threads)
    if parallel:
        script = destination / "orthomcl.pl"
        script.write_text(parallelize_source(script.read_text()))
    for item in originals:
        check(item)
    return {"originals": originals, "configured_sources": [record(destination / name)
            for name in ("orthomcl.pl", "orthomcl_module.pm")], "threads": threads,
            "pair_parallel_patch": parallel, "data_directory": str(data), "tool_directory": str(destination)}


def prepare(root, output):
    if output.exists():
        raise FileExistsError(output)
    results = root / "benchmark_tools/results"
    paths = [results / "qfo_corrected_orthomcl_perl_runtime_20260918.json",
             results / "qfo_corrected_orthomcl_system_helpers_20260918.json"]
    manifests = [read_frozen(path, sha) for path, sha in zip(paths, (RUNTIME_SHA, HELPERS_SHA))]
    for manifest in manifests:
        verify(manifest)
    checked = [*[record(path) for path in paths], record(__file__),
               *[record(Path(__file__).with_name(name)) for name in (
                   "configure_orthomcl_1_4.py", "parallelize_orthomcl_pairs.py", "run_orthomcl_perl_script.pl")]]
    base = root / "benchmarks/results/qfo_corrected_orthomcl_v1"
    sources = write_sources(base / "native_tool", base / "work", 180)
    for manifest in manifests:
        verify(manifest)
    for item in checked:
        check(item)
    report = {"status": "corrected_orthomcl_native_sources_prepared_unrun", **sources,
              "checked_records": checked, "pair_workers_planned": 64,
              "execution_authorized": False, "accuracy_admitted": False, "publication_ready": False,
              "limitations": [
                  "Configured paths/thread count and existing opt-in pair-parallel patch; no scientific default changes.",
                  "Source preparation is not syntax, native parity or production-output validation.",
                  "Guarded runtime, search admission, conversion/index evidence and scheduler binding remain required.",
                  "No existing pair caches or inference outputs may be reused implicitly."]}
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = prepare(args.root.resolve(), args.output.resolve())
    print(json.dumps({"status": report["status"]}))
