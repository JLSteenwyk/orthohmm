"""Freeze biological application commands without executing inference."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.launch_dgx_native_run import read_pinned
from benchmark_tools.snapshot_orthohmm_input_order import record

CORE_COMMIT = "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806"


def commands(core, input_dir, output, python, orthofinder, sonicparanoid):
    common = [str(python), str(core / "benchmark_tools/benchmark_production.py")]
    high = common + [str(input_dir), str(output / "high_sensitivity/output"),
                     str(output / "high_sensitivity/metrics.json"), "--cpu", "32",
                     "--threads-per-worker", "8", "--accuracy-profile", "high_sensitivity"]
    satellite = common + [str(input_dir), str(output / "satellite_v2/output"),
                          str(output / "satellite_v2/metrics.json"), *high[5:],
                          "--phylogeny", "reconcile", "--species-tree-mode", "infer",
                          "--phylogeny-candidates", "satellite_v2", "--phylogeny-root-rule",
                          "species_overlap", "--phylogeny-pair-rule", "positive_paralogy",
                          "--species-tree-rooting", "min_variance"]
    return [
        {"method": "orthohmm_high_sensitivity", "argv": high, "output_semantics": "native_orthogroups"},
        {"method": "orthohmm_satellite_v2", "argv": satellite, "output_semantics": "root_hogs"},
        {"method": "orthofinder_full", "argv": [str(orthofinder), "-f", str(output / "orthofinder/input"),
                                                 "-t", "32", "-a", "8", "-S", "diamond"],
         "copy_inputs_to": str(output / "orthofinder/input"), "output_semantics": "N0_root_hogs",
         "diagnostic_output": "MCL checkpoint, not a separate sequence-only run"},
        {"method": "sonicparanoid", "argv": [str(sonicparanoid), "-i", str(output / "sonicparanoid/input"),
                                               "-o", str(output / "sonicparanoid/output"), "-t", "32"],
         "copy_inputs_to": str(output / "sonicparanoid/input"), "output_semantics": "native_ortholog_groups"},
    ]


def prepare(repo):
    results = repo / "benchmark_tools/results"
    inputs = read_pinned(results / "biological_wgd_inputs_20260917.json",
                         "fe32a36d3e310c1589e48c1860410bcba277929422fcf6a39e0856d79b902042")
    contrast = results / "BIOLOGICAL_WGD_AUDIT_AND_CONTRASTS_20260917.md"
    if record(contrast)["sha256"] != "85ee0066d659285b1d89fb2ccd93cfbd726af74c902be20f7e039aaf444e58cc":
        raise ValueError("Changed contrast freeze")
    for item in inputs["inputs"] + [inputs["reference"], inputs["protocol"], inputs["cohort"]]:
        if record(item["path"]) != item:
            raise ValueError("Changed prepared application input")
    directory = repo / "benchmarks/work/biological_wgd_application_v1/input"
    if set(directory.iterdir()) != {Path(r["path"]) for r in inputs["inputs"]}:
        raise ValueError("Changed input membership")
    core = repo / "benchmarks/work/publication_method_native_v2"
    commit = subprocess.check_output(["git", "-C", str(core), "rev-parse", "HEAD"], text=True).strip()
    if commit != CORE_COMMIT:
        raise ValueError("Changed frozen inference commit")
    subprocess.run(["git", "-C", str(core), "diff", "--exit-code", "HEAD", "--",
                    "orthohmm", "benchmark_tools/benchmark_production.py"], check=True)
    output = repo / "benchmarks/results/biological_wgd_application_v1"
    if output.exists():
        raise FileExistsError(output)
    python = Path("/home/bizon/anaconda3/bin/python")
    sonic = python.with_name("sonicparanoid")
    of = repo.parents[1] / "SOFTWARE/orthofinder_3.1.5/bin/orthofinder"
    rows = commands(core, directory, output, python, of, sonic)
    eligible = [r for r in inputs["cohort_pairs"] if r["reference_eligible"]]
    if len(eligible) != 231 or any(not sum(n for s, n in r["available_members_by_species"].items()
                                         if s != "Scerevisiae") for r in eligible):
        raise ValueError("Changed fixed comparison population")
    return {"status": "commands_frozen_not_execution_authorized", "execution_authorized": False,
            "source": record(__file__), "core_commit": commit, "core_root": str(core),
            "inputs": record(results / "biological_wgd_inputs_20260917.json"),
            "contrast_protocol": record(contrast), "output_root": str(output), "runs": rows,
            "command_sources": [record(repo / "benchmark_tools/run_ygob_validation.slurm"),
                                record(repo / "qfo_benchmark/run_tool.slurm")],
            "entrypoints": [record(p) for p in (python, sonic, of)],
            "resources": {"cpus": 32, "memory_gib": 128, "host": "bizon", "sequential": True},
            "requirements": ["Fresh output paths and fresh comparator input copies preserving basenames/bytes.",
                             "Verify full runtime/dependencies and effective child PATH before and after each run.",
                             "Record actual OrthoHMM native input enumeration without changing frozen ordering.",
                             "Pin launcher and execution environment before authorization.",
                             "Native completion/IDs/group semantics must pass separate output admission.",
                             "Keep command failures; do not erase/restart or silently change configurations.",
                             "No controlled timing claim; run on original host, not dedicated DGX."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = prepare(args.repo.resolve())
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
