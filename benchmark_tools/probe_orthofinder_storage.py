"""Bounded preparation-only check of OrthoFinder input/output storage separation."""

import argparse
from pathlib import Path
import subprocess

from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_simulation_methods import read_frozen, execution_environment, verify_environment
from benchmark_tools.run_simulation_tree_mode_control import METHOD_SHA
from benchmark_tools.snapshot_orthohmm_input_order import record


def validate(inputs, output, before):
    after = [record(p) for p in sorted(inputs.iterdir())]
    if before != after:
        raise ValueError("Input inventory changed or acquired generated files")
    working = list(output.glob("Results_*/WorkingDirectory"))
    if len(working) != 1 or not working[0].is_dir():
        raise ValueError("Require one explicit output WorkingDirectory")
    for path in output.rglob("*"):
        if path.is_symlink() or not path.resolve().is_relative_to(output.resolve()):
            raise ValueError("Indirect output path")
    mapping = working[0] / "SpeciesIDs.txt"
    names = []
    for line in mapping.read_text().splitlines():
        index, name = line.split(": ", 1)
        if int(index) != len(names):
            raise ValueError("Noncontiguous species mapping")
        names.append(name)
    if names != sorted(Path(r["path"]).name for r in before):
        raise ValueError("Native species enumeration differs")
    if any(not (working[0] / f"Species{i}.fa").is_file() for i in range(len(names))):
        raise ValueError("Missing prepared native FASTA")
    return dict(native_order=names, working_directory=str(working[0]),
                unchanged_inputs=after, outputs=[record(p) for p in sorted(output.rglob("*")) if p.is_file()])


def run(inputs, directory):
    inputs, directory = inputs.absolute(), directory.absolute()
    if inputs.exists() or directory.exists():
        raise FileExistsError("Require fresh input and diagnostic directories")
    if not inputs.is_relative_to(Path("/dev/shm")) or directory.is_relative_to(Path("/dev/shm")):
        raise ValueError("Require tmpfs input candidate and separate disk output")
    baseline_path = Path(__file__).parent / "results/publication_variable_native_methods_20260916.json"
    baseline = read_frozen(baseline_path, METHOD_SHA)
    verify_environment(baseline)
    env, _ = execution_environment(baseline)
    inputs.mkdir()
    directory.mkdir()
    # Synthetic placement fixture, not biological evidence or a timing dataset.
    for name in ("species_b.fa", "species_a.fa"):
        (inputs / name).write_text(">protein_1\nMALWMRLLPLLALLALWGPGPGACDEFGHIKLMNPQRSTVWY\n")
    before = [record(p) for p in sorted(inputs.iterdir())]
    output = directory / "results"
    command = [baseline["tool_entrypoints"]["orthofinder"]["absolute_path"],
               "-f", str(inputs), "-o", str(output), "-t", "1", "-a", "1", "-S", "diamond", "-op"]
    save(directory / "command.json", dict(argv=command, baseline=record(baseline_path), source=record(__file__)))
    with (directory / "native.log").open("x") as log:
        result = subprocess.run(command, cwd=baseline["core_root"], env=env,
                                stdout=log, stderr=subprocess.STDOUT, timeout=120)
    save(directory / "exit.json", dict(exit_code=result.returncode))
    if result.returncode:
        raise ValueError("Preparation-only control failed; retained without retry")
    evidence = validate(inputs, output, before)
    save(directory / "result.json", dict(status="explicit_disk_output_preparation_verified", **evidence,
        scientific_timings_admitted=False, full_pipeline_validated=False,
        limitations=["Preparation-only synthetic fixture; no searches or phylogeny.",
                     "Does not measure tmpfs charge ownership, overhead or full-run temporary-file placement."]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", required=True, type=Path)
    parser.add_argument("--directory", required=True, type=Path)
    args = parser.parse_args()
    run(args.inputs, args.directory)
