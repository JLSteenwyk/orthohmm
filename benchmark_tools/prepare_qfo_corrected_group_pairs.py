"""Convert admitted corrected R-off candidates to cross-species QfO pairs."""

import argparse
from collections import Counter
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.bootstrap_qfo_factorial import CELLS
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_qfo_corrected_factorial_cell import verify_admission
from benchmark_tools.prepare_qfo_corrected_comparator_pairs import ENV_SHA
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.qfo_filter_pairs import load_mapping, filter_pairs
from qfo_benchmark.og_to_pairwise import index_genes, _strip_to_uniprot


def select(manifest, index):
    if type(index) is not int or index not in (0, 2, 4, 6):
        raise ValueError("Only the four R-off cells use group-derived pairs")
    if [c["label"] for c in manifest["cells"]] != list(CELLS):
        raise ValueError("Wrong corrected factorial cell inventory")
    cell = manifest["cells"][index]
    if cell["reconciliation"] is not False:
        raise ValueError("Cannot replace native phylogenetic pairs with group cliques")
    arm = manifest["candidate_arms"][cell["label"].rsplit("_", 1)[0]]
    if cell["candidate_partition"] != arm["candidate_partition"]["path"]:
        raise ValueError("Cell partition differs from admitted arm")
    return cell, arm


def expected_pairs(partition, owners):
    """Count cross-species pairs independently using per-family species sizes."""
    observed, total = set(), 0
    with partition.open() as stream:
        for line in stream:
            counts = Counter()
            for gene in line.split():
                accession = _strip_to_uniprot(gene)
                if accession in observed or gene not in owners:
                    raise ValueError("Duplicate or unknown partition identifier")
                observed.add(accession)
                counts[owners[gene]] += 1
            n = sum(counts.values())
            total += (n * n - sum(v * v for v in counts.values())) // 2
    universe = {_strip_to_uniprot(g) for g in owners}
    if observed != universe:
        raise ValueError("Incomplete candidate gene coverage")
    return total


def prepare(root, index, admission_path, admission_sha, admission_job):
    if type(index) is not int or index not in (0, 2, 4, 6):
        raise ValueError("Require R-off index 0, 2, 4 or 6")
    if not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "2":
        raise ValueError("Require scheduled 2-CPU conversion")
    admission, manifest, _, _, _, scheduler = verify_admission(
        root, admission_path, admission_sha, admission_job, index // 2)
    cell, arm = select(manifest, index)
    environment_path = root / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment = read_frozen(environment_path, ENV_SHA)
    mappings = [r for r in environment["reference_files"] if Path(r["path"]).name == "mapping.json.gz"]
    if len(mappings) != 1:
        raise ValueError("Require one frozen QfO reference mapping")
    mapping = mappings[0]
    fastas = manifest["input_fastas"]
    parents = {Path(r["path"]).parent for r in fastas}
    if len(parents) != 1:
        raise ValueError("Ambiguous input FASTA directory")
    fasta = next(iter(parents))
    if {str(p) for p in fasta.glob("*.fasta")} != {r["path"] for r in fastas}:
        raise ValueError("Uninventoried FASTA files")
    converter = Path(__file__).resolve().parent.parent / "qfo_benchmark/og_to_pairwise.py"
    sources = [record(__file__), record(converter), record(Path(__file__).with_name("qfo_filter_pairs.py")),
               record(Path(__file__).with_name("run_qfo_corrected_factorial_cell.py"))]
    checked = [record(admission_path), admission["prepared_manifest"], record(environment_path), mapping,
               arm["candidate_partition"], *fastas, *sources]
    for item in checked:
        check(item)
    partition = Path(arm["candidate_partition"]["path"])
    expected = expected_pairs(partition, index_genes(fasta))
    output = root / "benchmarks/results/qfo_corrected_factorial_pairs_v1" / cell["label"]
    output.mkdir(parents=True, exist_ok=False)
    command = [sys.executable, str(converter), str(partition), str(fasta)]
    report = {"status": "preparing", "index": index, "cell": cell["label"],
        "participant": "ohmm_qfo_corrected_factorial_" + cell["label"], "job_id": os.environ["SLURM_JOB_ID"],
        "candidate_admission": record(admission_path), "admission_scheduler": scheduler,
        "prepared": admission["prepared_manifest"], "candidate_partition": arm["candidate_partition"],
        "input_fastas": fastas, "mapping": mapping, "command": command, "checked_records": checked,
        "semantics": "cross-species group-derived clique pairs", "expected_pairs": expected,
        "accuracy_evaluated": False, "publication_ready": False}
    try:
        pairs, filtered = output / "pairs.partial.tsv", output / "pairs.qfo.partial.tsv"
        with pairs.open("x") as stream, (output / "conversion.log").open("x") as log:
            subprocess.run(command, stdout=stream, stderr=log, check=True)
        total, retained = filter_pairs(pairs, filtered, load_mapping(Path(mapping["path"])))
        if not 0 < total == retained == expected:
            raise ValueError("Pair count mismatch or unexpected corrected-input mapping loss")
        for item in checked:
            check(item)
        verify_admission(root, admission_path, admission_sha, admission_job, index // 2)
        final, mapped = output / "pairs.tsv", output / "pairs.qfo.tsv"
        pairs.rename(final)
        filtered.rename(mapped)
        report.update(status="corrected_factorial_group_pairs_prepared_unscored", pairs=record(final),
            filtered_pairs=record(mapped), total_pairs=total, retained_pairs=retained, removed_mapping_pairs=0,
            conversion_log=record(output / "conversion.log"))
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with (output / "results.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "admission"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("admission-sha256", "admission-job"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--index", type=int, choices=(0, 2, 4, 6), required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.index, args.admission.resolve(), args.admission_sha256, args.admission_job)
