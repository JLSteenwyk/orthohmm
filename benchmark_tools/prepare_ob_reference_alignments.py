"""Rebuild outcome-independent RefOG alignments with explicit normalization."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
import numpy as np
from benchmark_tools.assemble_orthobench_factorial import load_reference_snapshot
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.prepare_ob_sequence_features import CANONICAL
from benchmark_tools.recover_factorial_postflight import PREPARED_HASH, ENVIRONMENT_HASH
from benchmark_tools.run_simulation_methods import read_frozen


def normalize(sequence):
    sequence = sequence.upper()
    if not sequence or not set(sequence) <= set(CANONICAL + "X*"):
        raise ValueError("Unexpected reference protein alphabet or empty sequence")
    normalized = sequence.replace("*", "")
    if not normalized:
        raise ValueError("No residues after explicit stop removal")
    return normalized


def validate_alignment(path, expected):
    rows = list(SeqIO.parse(path, "fasta"))
    if len(rows) != len(expected) or {row.id for row in rows} != set(expected):
        raise ValueError("Alignment gene inventory differs or contains duplicates")
    if len({len(row.seq) for row in rows}) != 1 or not rows or not len(rows[0].seq):
        raise ValueError("Empty or unequal alignment lengths")
    aligned = {}
    for row in rows:
        sequence = str(row.seq).upper()
        if sequence.replace("-", "") != expected[row.id]:
            raise ValueError("Alignment changed normalized input residues")
        aligned[row.id] = sequence
    return aligned


def mean_canonical_identity(aligned):
    sequences = [aligned[name] for name in sorted(aligned)]
    if len(sequences) < 2 or len({len(s) for s in sequences}) != 1 or not sequences[0]:
        raise ValueError("Need at least two nonempty equally aligned sequences")
    encoded = np.frombuffer("".join(sequences).encode("ascii"), dtype=np.uint8).reshape(len(sequences), -1)
    valid = np.isin(encoded, np.frombuffer(CANONICAL.encode("ascii"), dtype=np.uint8))
    identities, missing = [], 0
    for i in range(len(sequences) - 1):
        comparable = valid[i] & valid[i + 1:]
        denominator = comparable.sum(axis=1)
        numerator = ((encoded[i] == encoded[i + 1:]) & comparable).sum(axis=1)
        missing += int(np.count_nonzero(denominator == 0))
        identities.extend((numerator[denominator > 0] / denominator[denominator > 0]).tolist())
    return {"mean_pairwise_identity": None if missing else float(np.mean(identities)),
            "pairs": len(sequences) * (len(sequences) - 1) // 2,
            "pairs_without_canonical_overlap": missing}


def tool_inventory(entry):
    prefix = entry.parent.parent / "libexec/mafft"
    if not prefix.is_dir():
        raise ValueError("Missing explicit MAFFT binary directory")
    paths = sorted({entry.resolve(), *[path.resolve() for path in prefix.rglob("*") if path.is_file()]})
    if len(paths) < 2:
        raise ValueError("Incomplete MAFFT executable inventory")
    return prefix, [file_provenance(path) for path in paths]


def align_family(name, genes, sequences, output, mafft, env):
    directory = output / Path(name).stem
    directory.mkdir()
    raw = {gene: sequences[gene] for gene in sorted(genes)}
    expected = {gene: normalize(sequence) for gene, sequence in raw.items()}
    input_path, alignment = directory / "input.faa", directory / "aligned.faa"
    with input_path.open("x") as handle:
        for gene, sequence in expected.items():
            handle.write(f">{gene}\n{sequence}\n")
    command = [str(mafft), "--amino", "--auto", "--thread", "1", str(input_path)]
    result = {"refog": name, "status": "running", "accuracy_evaluated": False, "command": command,
              "input": file_provenance(input_path), "removed_stop_symbols": {gene: sequence.count("*")
                  for gene, sequence in raw.items() if "*" in sequence}}
    started = time.monotonic()
    try:
        with alignment.open("x") as stdout, (directory / "mafft.stderr").open("x") as stderr:
            run = subprocess.run(command, env=env, cwd=directory, stdout=stdout, stderr=stderr)
        result["exit_code"] = run.returncode
        if run.returncode != 0:
            raise RuntimeError("MAFFT failed; no automatic retry")
        aligned = validate_alignment(alignment, expected)
        result.update(status="alignment_validated", genes=len(aligned), columns=len(next(iter(aligned.values()))),
                      alignment=file_provenance(alignment), identity=mean_canonical_identity(aligned))
        verify_file(input_path, result["input"])
    except Exception as error:
        result.update(status="failed", error_type=type(error).__name__, error=str(error))
    result["wall_s"] = time.monotonic() - started
    (directory / "status.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def prepare(root, output):
    root, output = root.resolve(), output.resolve()
    if output.exists():
        raise FileExistsError(output)
    results = root / "benchmark_tools/results"
    prepared_path = results / "orthobench_factorial_prepared_20260916.json"
    env_path = results / "publication_variable_native_methods_20260916.json"
    prepared = read_frozen(prepared_path, PREPARED_HASH)
    environment = read_frozen(env_path, ENVIRONMENT_HASH)
    references, _, _, reference_records = load_reference_snapshot(results / "orthobench_paired_uncertainty_20260916.json")
    wanted = set().union(*references.values())
    sequences, seen = {}, set()
    for item in prepared["fasta_inputs"]:
        verify_file(Path(item["path"]), item)
        for protein in SeqIO.parse(item["path"], "fasta"):
            if protein.id in seen:
                raise ValueError("Duplicate input gene identifier")
            seen.add(protein.id)
            if protein.id in wanted:
                sequences[protein.id] = str(protein.seq).upper()
                normalize(sequences[protein.id])
    if set(sequences) != wanted or len(seen) != 251378:
        raise ValueError("Incomplete input/reference gene coverage")
    mafft_record = environment["tool_entrypoints"]["mafft"]
    mafft = Path(mafft_record["absolute_path"])
    verify_file(mafft, mafft_record)
    prefix, tools_before = tool_inventory(mafft)
    overrides = {"MAFFT_BINARIES": str(prefix), "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                 "MKL_NUM_THREADS": "1", "LC_ALL": "C"}
    env = {**os.environ, **overrides}
    version = subprocess.run([str(mafft), "--version"], env=env, capture_output=True, text=True, check=True)
    version_text = (version.stdout + version.stderr).strip()
    if version_text != "v7.525 (2024/Mar/13)":
        raise ValueError("Unexpected MAFFT version")
    output.mkdir(parents=True)
    preflight = {"status": "ready_unscored", "accuracy_evaluated": False, "job_id": os.environ.get("SLURM_JOB_ID"),
                 "source": file_provenance(Path(__file__)), "input_manifest": file_provenance(prepared_path),
                 "environment_manifest": file_provenance(env_path), "tools": tools_before,
                 "version": version_text, "environment_overrides": overrides,
                 "reference_records": reference_records, "fasta_inputs": prepared["fasta_inputs"],
                 "reference_genes": len(wanted), "families": sorted(references), "workers": 8,
                 "normalization": "Uppercase; remove stop symbols explicitly, preserve X; do not modify inference FASTAs."}
    (output / "preflight.json").write_text(json.dumps(preflight, indent=2, sort_keys=True) + "\n")
    with ThreadPoolExecutor(max_workers=8) as executor:
        runs = list(executor.map(lambda name: align_family(name, references[name], sequences, output, mafft, env), sorted(references)))
    failed = [run["refog"] for run in runs if run["status"] != "alignment_validated"]
    if tool_inventory(mafft)[1] != tools_before:
        raise ValueError("MAFFT inventory changed during execution")
    for item in [*prepared["fasta_inputs"], *reference_records]:
        verify_file(Path(item["path"]), item)
    report = {"status": "failed_families_preserved" if failed else "reference_alignments_prepared_unscored",
              "accuracy_evaluated": False, "preflight": preflight, "runs": runs, "failed_families": failed,
              "limitations": ["Reference membership is used; no prediction/error outcomes define these alignments.",
                  "Canonical identity is a sequence descriptor, not evolutionary distance or duplication history.",
                  "Any pair without comparable canonical positions makes its family identity missing.",
                  "Shared-node alignment times are preparation costs, not orthology-inference performance."]}
    (output / "manifest.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if failed:
        raise RuntimeError("Reference alignment failures retained; no family silently dropped")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    prepare(args.root, args.output)


if __name__ == "__main__":
    main()
