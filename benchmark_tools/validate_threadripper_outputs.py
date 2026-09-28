"""Local scaling output checks, including exact native species/sequence numbering."""

from pathlib import Path

from Bio import SeqIO

from benchmark_tools.prepare_threadripper_run import check_prepared
from benchmark_tools.snapshot_orthohmm_input_order import record
from benchmark_tools.simulation_method_outputs import unique_path
from benchmark_tools.validate_scaling_outputs import validate as validate_native


def orthofinder_mapping(run):
    output = Path(run["configuration"]["output"])
    species_path = unique_path(output, "**/SpeciesIDs.txt")
    sequence_path = unique_path(output, "**/SequenceIDs.txt")
    if species_path.is_symlink() or sequence_path.is_symlink():
        raise ValueError("Indirect native mapping file")
    names = run["expected_native_order"]
    source_names = sorted(Path(r["path"]).name for r in run["dataset"]["inputs"])
    if names != source_names or len(set(names)) != len(names):
        raise ValueError("Expected sorted species inventory differs")
    rows = species_path.read_text().splitlines()
    if rows != [f"{i}: {name}" for i, name in enumerate(names)]:
        raise ValueError("Native species IDs differ from frozen enumeration")
    expected = {}
    directory = Path(run["prepared_input_directory"])
    for i, name in enumerate(names):
        for j, sequence in enumerate(SeqIO.parse(directory / name, "fasta")):
            expected[f"{i}_{j}"] = sequence.id
    observed = {}
    with sequence_path.open() as handle:
        for line in handle:
            key, value = line.rstrip("\n").split(": ", 1)
            fields = value.split()
            if key in observed or not fields:
                raise ValueError("Duplicate or empty native sequence mapping")
            observed[key] = fields[0]
    if not expected or observed != expected:
        raise ValueError("Native sequence IDs differ from FASTA order and identifiers")
    return dict(species=len(names), sequences=len(expected), native_order=names,
                checked_files=[record(species_path), record(sequence_path)])


def validate(run, measurement, baseline):
    prepared = check_prepared(run, baseline)
    result = validate_native(run, measurement)
    mapping = orthofinder_mapping(run) if run["native_method"] == "orthofinder_full" else None
    return dict(status="threadripper_native_outputs_checked", native=result,
                prepared_inputs=prepared, orthofinder_native_mapping=mapping,
                scientific_timings_admitted=False, accuracy_evaluated=False,
                limitations=["Caller must validate before/after runtime, command plan and preparation provenance.",
                             "Valid outputs and native enumeration do not establish host isolation or timing admission."])
