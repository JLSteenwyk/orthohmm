"""Reconstruct source-hit tuples independently of the numeric checkpoint writer."""

import csv
import math
from pathlib import Path
import sqlite3

import numpy as np

from orthohmm.accuracy import load_accuracy_checkpoint


def reconstruct(database, sources, metadata, names):
    if names != sorted(metadata) or len(names) != len(set(names)):
        raise ValueError("Require exact sorted source gene identities")
    indices = {name: i for i, name in enumerate(names)}
    database.execute("CREATE TABLE expected(q INTEGER,t INTEGER,species TEXT,raw REAL,score REAL,PRIMARY KEY(q,t)) WITHOUT ROWID")
    count, batch = 0, []
    for path, species in sources:
        with Path(path).open(newline="") as stream:
            for fields in csv.reader(stream, delimiter="\t"):
                if len(fields) != 7:
                    raise ValueError("Source hit requires seven fields")
                query, target = fields[:2]
                if query not in metadata or target not in metadata:
                    raise ValueError("Unknown source hit ID")
                qlength, tlength = int(fields[2]), int(fields[3])
                if (qlength != metadata[query]["length"] or tlength != metadata[target]["length"]
                        or min(qlength, tlength) <= 0 or metadata[target]["species"] != species):
                    raise ValueError("Source length or target ownership differs")
                raw, bits, evalue = map(float, fields[4:])
                if (not all(math.isfinite(v) for v in (raw, bits, evalue))
                        or min(raw, bits) <= 0 or not 0 <= evalue <= 1e-4):
                    raise ValueError("Invalid source score")
                batch.append((indices[query], indices[target], species, raw, raw / math.sqrt(qlength * tlength)))
                count += 1
                if len(batch) == 10000:
                    database.executemany("INSERT INTO expected VALUES(?,?,?,?,?)", batch)
                    batch.clear()
        database.executemany("INSERT INTO expected VALUES(?,?,?,?,?)", batch)
        batch.clear()
    database.commit()
    if count <= 0:
        raise ValueError("Empty source search panel")
    database.execute("CREATE INDEX expected_rank ON expected(q,species,raw DESC,t)")
    return count


def expected_top100(database):
    # Buffer at most 100 hits per target species for one query, then restore q,t order.
    current_query, current_species, rank, selected = None, None, 0, []
    for query, target, species, score in database.execute("SELECT q,t,species,score FROM expected ORDER BY q,species,raw DESC,t"):
        if query != current_query:
            yield from sorted(selected)
            current_query, current_species, selected = query, None, []
        if species != current_species:
            current_species, rank = species, 0
        rank += 1
        if rank <= 100:
            selected.append((query, target, score))
    yield from sorted(selected)


def compare_rows(rows, queries, targets, scores, chunk_size=100000):
    if type(chunk_size) is not int or chunk_size < 1:
        raise ValueError("Positive chunk size required")
    if queries.shape != targets.shape or queries.shape != scores.shape or queries.ndim != 1:
        raise ValueError("Checkpoint hit array shapes differ")
    if queries.dtype != np.int32 or targets.dtype != np.int32 or scores.dtype != np.float64:
        raise ValueError("Checkpoint hit dtypes differ")
    offset, batch = 0, []

    def compare(batch, offset):
        end = offset + len(batch)
        if end > len(scores):
            raise ValueError("Checkpoint omits source hits")
        expected = np.asarray(batch)
        if (not np.array_equal(queries[offset:end], expected[:, 0])
                or not np.array_equal(targets[offset:end], expected[:, 1])
                or not np.array_equal(scores[offset:end], expected[:, 2])):
            raise ValueError("Checkpoint tuples differ from source reconstruction")
        return end

    for row in rows:
        batch.append(row)
        if len(batch) == chunk_size:
            offset = compare(batch, offset)
            batch.clear()
    if batch:
        offset = compare(batch, offset)
    if offset != len(scores):
        raise ValueError("Checkpoint contains additional hits")
    return offset


def verify(database, checkpoint, metadata, cap):
    if cap not in (None, 100):
        raise ValueError("Only all-hit and top100 variants are specified")
    names, species, queries, targets, scores = load_accuracy_checkpoint(checkpoint, verify=True)
    if names != sorted(metadata):
        raise ValueError("Checkpoint gene identities differ")
    owners = {name: i for i, name in enumerate(sorted({r["species"] for r in metadata.values()}))}
    expected_species = np.asarray([owners[metadata[name]["species"]] for name in names], dtype=np.int32)
    if species.dtype != np.int32 or not np.array_equal(species, expected_species):
        raise ValueError("Checkpoint species indices differ")
    rows = database.execute("SELECT q,t,score FROM expected ORDER BY q,t") if cap is None else expected_top100(database)
    count = compare_rows(rows, queries, targets, scores)
    return {"status": "checkpoint_matches_reconstructed_source_hits", "cap": cap,
            "genes": len(names), "hits": count, "accuracy_evaluated": False}
