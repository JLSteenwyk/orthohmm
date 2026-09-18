"""Create a sorted diagnostic copy of an admitted numeric search checkpoint."""

import argparse
import json
from pathlib import Path
import sqlite3
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_accuracy_checkpoint import audit
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.summarize_sorted_search_hits import blocks
from orthohmm.accuracy import load_accuracy_checkpoint, write_accuracy_checkpoint


def ingest(database, queries, targets, scores, n, chunk_size=100000):
    if type(chunk_size) is not int or chunk_size < 1:
        raise ValueError("Positive integer chunk size required")
    if (type(n) is not int or n < 1 or any(a.ndim != 1 for a in (queries, targets, scores))
            or not len(queries) == len(targets) == len(scores)
            or queries.dtype != np.int32 or targets.dtype != np.int32 or scores.dtype != np.float64):
        raise ValueError("Invalid native hit arrays")
    database.execute("CREATE TABLE hits(q INTEGER,t INTEGER,score REAL,PRIMARY KEY(q,t)) WITHOUT ROWID")
    for start in range(0, len(scores), chunk_size):
        q, t, s = (a[start:start + chunk_size] for a in (queries, targets, scores))
        if (np.any(q < 0) or np.any(t < 0) or np.any(q >= n) or np.any(t >= n)
                or not np.isfinite(s).all() or np.any(s <= 0)):
            raise ValueError("Invalid checkpoint tuple")
        database.executemany("INSERT INTO hits VALUES(?,?,?)",
                             ((int(a), int(b), float(c)) for a, b, c in zip(q, t, s)))
    database.commit()
    count = database.execute("SELECT COUNT(*) FROM hits").fetchone()[0]
    if count != len(scores) or database.execute("PRAGMA quick_check").fetchall() != [("ok",)]:
        raise ValueError("Canonical database integrity/count differs")
    return count


def write_sorted(database, names, species, output, count, chunk_size=100000):
    paths = [output / name for name in ("sorted_q.npy", "sorted_t.npy", "sorted_s.npy")]
    arrays = [np.lib.format.open_memmap(path, mode="w+", dtype=dtype, shape=(count,))
              for path, dtype in zip(paths, (np.int32, np.int32, np.float64))]
    cursor = database.execute("SELECT q,t,score FROM hits ORDER BY q,t")
    offset = 0
    while rows := cursor.fetchmany(chunk_size):
        for column, array in enumerate(arrays):
            array[offset:offset + len(rows)] = [r[column] for r in rows]
        offset += len(rows)
    if offset != count:
        raise ValueError("Sorted database count changed")
    for array in arrays:
        array.flush()
    for _ in blocks(*arrays, len(names), chunk_size):
        pass
    return write_accuracy_checkpoint(str(output), names, species, *arrays)


def canonicalize(source, expected_sha, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    admission = audit(source, expected_sha)
    names, species, q, t, s = load_accuracy_checkpoint(source, verify=False)
    if names != sorted(names):
        raise ValueError("Canonical diagnostic requires lexical gene IDs")
    checked = [record(path) for path in sorted(source.iterdir())]
    output.mkdir(parents=True, exist_ok=False)
    result = {"status": "canonicalizing", "source_checkpoint_audit": admission,
        "source": record(__file__), "checked_records": checked, "accuracy_evaluated": False,
        "publication_ready": False, "helpers": [record(Path(__file__).with_name(name)) for name in
            ("audit_accuracy_checkpoint.py", "summarize_sorted_search_hits.py")],
        "checkpoint_writer": record(Path(__file__).resolve().parent.parent / "orthohmm/accuracy.py")}
    try:
        with sqlite3.connect(output / "canonical.sqlite") as database:
            database.execute("PRAGMA temp_store=FILE")
            database.execute("PRAGMA cache_size=-131072")
            count = ingest(database, q, t, s, len(names))
            checkpoint = write_sorted(database, names, species, output, count)
        manifest = record(checkpoint / "manifest.json")
        result.update(status="canonical_search_checkpoint_written", hits=count, genes=len(names),
                      checkpoint=str(checkpoint), manifest=manifest, audit=audit(checkpoint, manifest["sha256"]),
                      limitations=["Diagnostic copy only; original native checkpoint remains unchanged.",
                          "Requires separate upstream scientific provenance admission.",
                          "Duplicate ordered pairs fail; no deduplication, score adjustment or gene relabeling.",
                          "Scratch SQLite and sorted arrays are retained; disk usage exceeds output checkpoint size."])
        for item in [*checked, result["source"], *result["helpers"], result["checkpoint_writer"]]:
            check(item)
    except BaseException as error:
        result.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with (output / "manifest.json").open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    canonicalize(args.checkpoint.resolve(), args.sha256, args.output.absolute())
