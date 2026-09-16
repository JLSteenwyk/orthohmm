"""Validate sequence-search hits and create deterministic numeric replay checkpoints."""

import argparse
import json
import math
from pathlib import Path
import sqlite3
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_accuracy_checkpoint import audit
from benchmark_tools.audit_historical_profile_ablation import verify_file
from benchmark_tools.orthobench_stage_diagnostics import file_provenance
from benchmark_tools.run_sequence_search_control import MANIFEST_SHA, verify_plan
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_ygob_validation import require_completed_job
from orthohmm.accuracy import write_accuracy_checkpoint


def parse_hits(path, metadata, indices, target_species, species_index):
    with path.open() as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 7:
                raise ValueError(f"Expected seven fields at {path}:{number}")
            query, target = fields[:2]
            if query not in metadata or target not in metadata:
                raise ValueError("Unknown query or target ID")
            qlen, tlen = map(int, fields[2:4])
            if qlen != metadata[query]["length"] or tlen != metadata[target]["length"] or min(qlen, tlen) <= 0:
                raise ValueError("Sequence length differs from input")
            if metadata[target]["species"] != target_species:
                raise ValueError("Target belongs to a different species database")
            raw, bits, evalue = map(float, fields[4:])
            if not all(math.isfinite(x) for x in (raw, bits, evalue)) or raw <= 0 or bits <= 0 or not 0 <= evalue <= 1e-4:
                raise ValueError("Invalid score or significance threshold")
            yield indices[query], indices[target], species_index, raw, raw / math.sqrt(qlen * tlen)


def initialize(database):
    database.execute("CREATE TABLE hits(q INTEGER, t INTEGER, species INTEGER, raw REAL, score REAL, PRIMARY KEY(q,t)) WITHOUT ROWID")


def ingest(database, rows, batch_size=10000):
    batch = []
    for row in rows:
        batch.append(row)
        if len(batch) >= batch_size:
            database.executemany("INSERT INTO hits VALUES(?,?,?,?,?)", batch)
            batch.clear()
    database.executemany("INSERT INTO hits VALUES(?,?,?,?,?)", batch)
    database.commit()


def selected_sql(cap):
    if cap is None:
        return "SELECT q,t,score FROM hits ORDER BY q,t"
    if cap != 100:
        raise ValueError("Only the frozen top100 diagnostic is supported")
    return ("SELECT q,t,score FROM (SELECT q,t,score, ROW_NUMBER() OVER "
            "(PARTITION BY q,species ORDER BY raw DESC,t) AS rank FROM hits) WHERE rank<=100 ORDER BY q,t")


def write_variant(database, directory, names, species, cap):
    sql = selected_sql(cap)
    count = database.execute("SELECT COUNT(*) FROM (" + sql + ")").fetchone()[0]
    directory.mkdir()
    arrays = [np.lib.format.open_memmap(directory / filename, mode="w+", dtype=dtype, shape=(count,))
              for filename, dtype in (("queries.npy", np.int32), ("targets.npy", np.int32), ("scores.npy", np.float64))]
    cursor, offset = database.execute(sql), 0
    while batch := cursor.fetchmany(100000):
        values = np.asarray(batch)
        for column, array in enumerate(arrays):
            array[offset:offset + len(batch)] = values[:, column]
        offset += len(batch)
    if offset != count:
        raise ValueError("Hit count changed during conversion")
    for array in arrays:
        array.flush()
    checkpoint = write_accuracy_checkpoint(str(directory), names, species, *arrays)
    evidence = file_provenance(checkpoint / "manifest.json")
    return {"checkpoint": str(checkpoint), "manifest": evidence, "audit": audit(checkpoint, evidence["sha256"]), "cap": cap}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = read_frozen(args.manifest, MANIFEST_SHA)
    root = verify_plan(report)
    status_path = root / "execution.json"
    status_record = file_provenance(status_path)
    status = json.loads(status_path.read_text())
    accounting = subprocess.check_output(["sacct", "-j", "21291", "--parsable2",
                                         "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 21291)
    if status["job_id"] != "21291" or status["status"] != "complete_pending_numeric_validation" or len(status["targets"]) != 12:
        raise ValueError("Full frozen search panel has not completed")
    verify_file(args.manifest, status["manifest"])
    verify_file(Path(status["executor"]["path"]), status["executor"])
    for cell, result in zip(report["searches"], status["targets"], strict=True):
        if result["index"] != cell["index"] or result["status"] != "complete_pending_numeric_validation":
            raise ValueError("Unexpected target outcome")
        for phase in ("makedb", "search"):
            evidence = result["phases"][phase]
            if evidence["argv"] != cell[phase] or evidence["exit_code"] != 0:
                raise ValueError("Search phase failed or command changed")
        if result["hits"]["path"] != cell["output"]:
            raise ValueError("Hit path differs")
        verify_file(Path(cell["output"]), result["hits"])
    metadata = json.loads(Path(report["gene_metadata"]["path"]).read_text())
    names = sorted(metadata)
    indices = {gene: i for i, gene in enumerate(names)}
    species_names = sorted({r["species"] for r in metadata.values()})
    species_ids = {name: i for i, name in enumerate(species_names)}
    species = np.array([species_ids[metadata[g]["species"]] for g in names], dtype=np.int32)
    args.output.mkdir(parents=True)
    with sqlite3.connect(args.output / "hits.sqlite") as database:
        initialize(database)
        for cell in report["searches"]:
            target = Path(cell["target_fasta"]["path"]).name
            ingest(database, parse_hits(Path(cell["output"]), metadata, indices, target, species_ids[target]))
        database.execute("CREATE INDEX hit_rank ON hits(q,species,raw DESC,t)")
        variants = {label: write_variant(database, args.output / label, names, species, cap)
                    for label, cap in (("all_hits", None), ("top100", 100))}
    for result in status["targets"]:
        verify_file(Path(result["hits"]["path"]), result["hits"])
    verify_plan(report)
    verify_file(status_path, status_record)
    output = {"schema_version": 1, "status": "numeric_checkpoints_verified", "accuracy_evaluated": False,
              "scheduler": scheduler, "execution": status_record, "variants": variants,
              "source": file_provenance(Path(__file__)), "sqlite_version": sqlite3.sqlite_version,
              "normalization": report["normalization"], "tie_break": "descending raw score, then lexical target gene ID"}
    (args.output / "manifest.json").write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
