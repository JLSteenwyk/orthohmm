"""Freeze the corrected-input Proteinortho command and container identity."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_qfo_corrected_primary import verify

PRIMARY_SHA = "dbebd3a6915fddeb2b89e0e591a2ac5c798ee6bccec9925fef485260f41baa5a"
REGISTRY_SHA = "094f842ad8211450519d4238c0b3a5020a7969049d4b7306f5465d5acbf914a3"
IMAGE = Path("/mnt/ca1e2e99-718e-417c-9ba6-62421455971a/SOFTWARE/proteinortho_6.3.6--h2b77389_0.sif")
RUNTIME = Path("/usr/local/bin/singularity")


def native_command(runtime, image, names):
    names = list(names)
    if len(names) != 78 or len(set(names)) != 78:
        raise ValueError("Require exactly 78 distinct input filenames")
    if any(Path(n).name != n or not n.endswith(".fasta") or n.startswith("-") for n in names):
        raise ValueError("Require canonical FASTA basenames")
    return [str(runtime), "exec", "--bind", "/mnt", str(image), "proteinortho",
            "-project=qfo", "-cpus=32", *sorted(names)]


def probe(argv, cwd):
    done = subprocess.run(argv, cwd=cwd, capture_output=True, text=True, timeout=60, check=True)
    return {"argv": argv, "exit_code": done.returncode, "stdout": done.stdout, "stderr": done.stderr}


def prepare(root, output, destination):
    if output.exists() or destination.exists():
        raise FileExistsError("Require fresh output root and manifest")
    results = root / "benchmark_tools/results"
    primary_path = results / "qfo_corrected_primary_commands_20260918.json"
    primary = read_frozen(primary_path, PRIMARY_SHA)
    _, inputs = verify(primary)
    registry = results / "publication_comparison_orthomcl_complete_20260916.json"
    read_frozen(registry, REGISTRY_SHA)
    sources = [record(root / path) for path in (
        "qfo_benchmark/run_tool.slurm", "benchmark_tools/proteinortho_to_pairwise.py",
        "bacterial_scaling/parse_proteinortho.py", "qfo_benchmark/filter_qfo_pairs.py",
        "benchmark_tools/run_qfo_corrected_primary.py")]
    fixed = [record(IMAGE), record(RUNTIME), record(primary_path), record(registry), *sources, *inputs]
    probes = {
        "runtime": probe([str(RUNTIME), "--version"], root),
        "proteinortho": probe([str(RUNTIME), "exec", "--bind", "/mnt", str(IMAGE), "proteinortho", "-version"], root),
        "diamond": probe([str(RUNTIME), "exec", "--bind", "/mnt", str(IMAGE), "diamond", "version"], root),
    }
    if probes["proteinortho"]["stdout"].strip() != "6.3.6":
        raise ValueError("Unexpected Proteinortho version")
    for item in fixed:
        check(item)
    report = {
        "status": "corrected_proteinortho_command_frozen_unrun", "execution_authorized": False,
        "accuracy_admitted": False, "source": record(__file__), "checked_records": fixed,
        "primary_manifest": record(primary_path), "image": record(IMAGE), "runtime": record(RUNTIME),
        "probes": probes, "input_directory": primary["input_directory"], "input_fastas": inputs,
        "output_root": str(output), "cwd": str(output / "input"),
        "native_argv": native_command(RUNTIME, IMAGE, [Path(r["path"]).name for r in inputs]),
        "observed_environment": {k: v for k, v in os.environ.items()
                                 if k.startswith(("SINGULARITY", "APPTAINER")) or k in
                                 ("PATH", "LD_LIBRARY_PATH", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
        "resources": {"cpus": 32, "memory_gib": 192, "node": "bizon", "time_limit_hours": 72},
        "native_outputs": {"groups": "input/qfo.proteinortho.tsv", "pairs": "input/qfo.proteinortho-graph"},
        "pair_semantics": "Native post-clustering graph; not blast-graph edges or group cliques.",
        "search_reuse": False,
        "remaining_gates": [
            "Pinned runner, runtime configuration and inherited/container environment verification before launch.",
            "Fresh byte-identical input copies; no original-release intermediate reuse.",
            "Pre/post-run provenance checks, native completeness and graph validation before conversion/scoring.",
            "Freeze exact conversion and corrected participant/scoring commands before evaluation.",
        ],
        "limitations": ["Current version probes are not historical child-process execution proof.",
                        "Shared-host resources are descriptive, not dedicated matched timings.",
                        "Container image identity does not by itself pin host runtime configuration or bind-mounted dependencies."],
    }
    with destination.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "output-root", "manifest"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.output_root.resolve(), args.manifest.resolve())
