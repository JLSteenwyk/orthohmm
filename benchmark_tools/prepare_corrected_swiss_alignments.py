"""Prepare outcome-independent corrected SwissTrees alignments and identity bins."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
import os
from pathlib import Path
import statistics
import subprocess

from Bio import SeqIO
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_ob_reference_alignments import align_family, normalize, tool_inventory
from benchmark_tools.run_simulation_methods import read_frozen

INVENTORY_SHA = "912363276fa5456f11aeff9f398fce71aec8d377c8c8eda7b144105945e554f1"
ENVIRONMENT_SHA = "bf677728cf9312edd9755d9ac1ceefbce3abb8c3848d43413631dc00ca00eb2f"
HELPER_SHA = "9a787bffec2a9481c8987a7e7c5066982c1ec3fe4c3f571cf7a02e69162e885d"


def identity_bins(runs):
    names = [r["refog"] for r in runs]
    if not names or len(set(names)) != len(names) or any(r["status"] != "alignment_validated" for r in runs):
        raise ValueError("Require complete unique successful family panel")
    values = {r["refog"]: r["identity"]["mean_pairwise_identity"] for r in runs}
    known = [value for value in values.values() if value is not None]
    if any(isinstance(v, bool) or not math.isfinite(v) or not 0 <= v <= 1 for v in known):
        raise ValueError("Invalid identity")
    cutoff = statistics.median(known) if known else None
    bins = {key: [] for key in ("lower_identity", "higher_identity", "missing_identity")}
    for name, value in sorted(values.items()):
        bins["missing_identity" if value is None else "lower_identity" if value <= cutoff else "higher_identity"].append(name)
    return dict(median_family_identity=cutoff, strata=bins)


def extract_sequences(inventory):
    families = inventory["family_memberships"]
    wanted = {gene for genes in families.values() for gene in genes}
    sequences = {}
    for item in inventory["fasta_inputs"]:
        check(item)
        for row in SeqIO.parse(item["path"], "fasta"):
            fields = row.id.split("|")
            if len(fields) != 3:
                raise ValueError("Unexpected accession header")
            gene = fields[1]
            if gene in wanted:
                if gene in sequences:
                    raise ValueError("Duplicate reference accession")
                sequence = str(row.seq).upper()
                if len(sequence) != inventory["genes"][gene]["length"]:
                    raise ValueError("Reference length changed")
                normalize(sequence)
                sequences[gene] = sequence
        check(item)
    if set(sequences) != wanted:
        raise ValueError("Incomplete reference sequence coverage")
    return sequences


def prepare(root, output, protocol_sha256):
    root, output = root.resolve(), output.absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    results = root / "benchmark_tools/results"
    inventory_path = results / "corrected_swiss_sequence_strata_20260918.json"
    environment_path = results / "publication_variable_native_methods_20260916.json"
    inventory = read_frozen(inventory_path, INVENTORY_SHA)
    environment = read_frozen(environment_path, ENVIRONMENT_SHA)
    protocol = record(results / "CORRECTED_SWISS_IDENTITY_PROTOCOL_20260923.md")
    helper = record(Path(__file__).with_name("prepare_ob_reference_alignments.py"))
    if protocol["sha256"] != protocol_sha256 or helper["sha256"] != HELPER_SHA:
        raise ValueError("Frozen protocol or helper changed")
    families = inventory["family_memberships"]
    if len(families) != 18 or len(inventory["fasta_inputs"]) != 78:
        raise ValueError("Wrong corrected reference panel")
    if any(not genes or len(set(genes)) != len(genes) or Path(name).name != name for name, genes in families.items()):
        raise ValueError("Invalid family membership")
    sequences = extract_sequences(inventory)
    if len(sequences) != 563 or sum(map(len, families.values())) != 563:
        raise ValueError("Wrong corrected reference universe")
    mafft_record = environment["tool_entrypoints"]["mafft"]
    mafft = Path(mafft_record["absolute_path"])
    if any(record(mafft)[k] != mafft_record[k] for k in ("bytes", "sha256")):
        raise ValueError("MAFFT entrypoint changed")
    prefix, tools = tool_inventory(mafft)
    overrides = dict(MAFFT_BINARIES=str(prefix), OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                     MKL_NUM_THREADS="1", LC_ALL="C")
    env = {**os.environ, **overrides}
    version = subprocess.run([str(mafft), "--version"], capture_output=True, text=True, env=env, check=True)
    if (version.stdout + version.stderr).strip() != "v7.525 (2024/Mar/13)":
        raise ValueError("Wrong MAFFT version")
    inputs = [record(inventory_path), record(environment_path), protocol, helper, record(__file__)]
    preflight = dict(status="ready_unscored", inputs=inputs, tools=tools, families=families,
        fasta_inputs=inventory["fasta_inputs"], workers=4, environment_overrides=overrides,
        job_id=os.environ.get("SLURM_JOB_ID"), prediction_statistics_evaluated=False)
    output.mkdir(parents=True)
    (output / "preflight.json").write_text(json.dumps(preflight, indent=2, sort_keys=True) + "\n")
    with ThreadPoolExecutor(max_workers=4) as executor:
        runs = list(executor.map(lambda name: align_family(name, families[name], sequences, output, mafft, env), sorted(families)))
    failed = [r["refog"] for r in runs if r["status"] != "alignment_validated"]
    for item in [*inputs, *inventory["fasta_inputs"]]:
        check(item)
    if tool_inventory(mafft)[1] != tools:
        raise ValueError("MAFFT inventory changed")
    report = dict(status="failed_families_preserved" if failed else "corrected_swiss_alignments_prepared_unscored",
        preflight=preflight, runs=runs, failed_families=failed, prediction_statistics_evaluated=False,
        publication_ready=False, **({} if failed else identity_bins(runs)))
    (output / "manifest.json").write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    if failed:
        raise RuntimeError("Failed families retained without retry")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--protocol-sha256", required=True)
    args = parser.parse_args()
    prepare(args.root, args.output, args.protocol_sha256)
