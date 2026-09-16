"""Bind QfO numeric hits to historical source, FASTAs and target partition."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from orthohmm.accuracy import load_accuracy_checkpoint
from benchmark_tools.audit_accuracy_checkpoint import audit
from benchmark_tools.orthobench_stage_diagnostics import file_provenance


def verify_species_partition(names, species, fasta_records):
    index = {name: i for i, name in enumerate(names)}
    seen, file_codes, code_files = set(), {}, {}
    for filename, identifiers in fasta_records:
        codes = set()
        for gene in identifiers:
            if gene not in index or gene in seen:
                raise ValueError("FASTA gene absent from checkpoint or duplicated")
            seen.add(gene)
            codes.add(int(species[index[gene]]))
        if len(codes) != 1:
            raise ValueError("Proteome does not map to exactly one species code")
        code = next(iter(codes))
        if code in code_files or filename in file_codes:
            raise ValueError("Species code or filename reused across proteomes")
        file_codes[filename], code_files[code] = code, filename
    if seen != set(names):
        raise ValueError("Checkpoint genes absent from FASTAs")
    return file_codes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    if args.output.exists():
        raise FileExistsError(args.output)
    run = root / "qfo_benchmark/results/orthohmm_high_sensitivity_isolated"
    metrics_path = run / "metrics.json"
    metrics_record = file_provenance(metrics_path)
    if metrics_record["sha256"] != "fb6b8d7e6824e1801c7cdb794e319d30ee9c16d9eada9213ef345b8de54d8893":
        raise ValueError("Historical metrics hash changed")
    metrics = json.loads(metrics_path.read_text())
    harness = metrics["harness"]
    if metrics["status"] != "complete" or harness["exit_code"] != 0:
        raise ValueError("Historical inference not complete")
    revision = "694a77fe56167754ca949751bca88aa7d11353dc"
    if harness["git_commit"] != revision:
        raise ValueError("Historical source revision mismatch")
    source_changes = []
    frozen = root / "benchmarks/work/publication_method_native_v2"
    for record in harness["source_manifest"]:
        data = subprocess.check_output(["git", "-C", str(root), "show", f"{revision}:{record['path']}"])
        if len(data) != record["bytes"] or hashlib.sha256(data).hexdigest() != record["sha256"]:
            raise ValueError("Historical source manifest differs from recorded commit")
        current = file_provenance(frozen / record["path"])
        if current["sha256"] != record["sha256"]:
            source_changes.append({"path": record["path"], "historical_sha256": record["sha256"], "frozen_sha256": current["sha256"]})
    fasta = root / "qfo_benchmark/input"
    inputs = harness["input_manifest"]
    if {p.name for p in fasta.glob("*.fasta")} != {r["path"] for r in inputs}:
        raise ValueError("FASTA file set mismatch")
    records = []
    for record in inputs:
        observed = file_provenance(fasta / record["path"])
        if (observed["bytes"], observed["sha256"]) != (record["bytes"], record["sha256"]):
            raise ValueError("Historical FASTA hash mismatch")
        records.append(observed)
    checkpoint = run / "output/orthohmm_working_res/high_sensitivity_checkpoint"
    numeric = audit(checkpoint, "b90c787f050a087adeeb81b1cace18c9cbb30e1724926a40a1df6d7e49bd549c")
    names, species, *_ = load_accuracy_checkpoint(checkpoint, verify=False)
    mapping = verify_species_partition(names, species, ((r["path"], (s.id for s in SeqIO.parse(fasta / r["path"], "fasta"))) for r in inputs))
    target = run / "output/orthohmm_working_res/orthohmm_edges_clustered.txt"
    target_record = file_provenance(target)
    if target_record["sha256"] != "63ade2f317c6343bd0d1e98af0fe0e299dd61530fb6094dec62d97dab30a4df3":
        raise ValueError("Historical target partition hash mismatch")
    report = {"schema_version": 1, "status": "historical_inputs_verified", "historical_revision": revision,
              "metrics": metrics_record, "input_fastas": records, "species_code_mapping": mapping,
              "numeric_checkpoint": numeric, "target_partition": target_record,
              "changed_recorded_sources_since_historical_run": source_changes,
              "auditor": file_provenance(Path(__file__)), "accuracy_evaluated": False,
              "replay_equivalence": "unproven; source differences require explicit frozen replay validation"}
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"Verified {len(names)} genes in {len(mapping)} FASTAs; {len(source_changes)} recorded source files differ from frozen core")


if __name__ == "__main__":
    main()
