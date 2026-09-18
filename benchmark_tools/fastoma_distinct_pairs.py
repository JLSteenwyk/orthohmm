"""Disk-backed distinct conversion of FastOMA's native orthologs.tsv.gz.

The caller supplies accession-to-species ownership from validated inputs and
owns output publication and provenance admission. This is not a HOG converter.
"""

from pathlib import Path
import sqlite3
import tempfile
from typing import Mapping, TextIO

from benchmark_tools.fastoma_to_pairwise import iter_pairs


def write_pairs(source: Path, owners: Mapping[str, str], output: TextIO,
                scratch: Path) -> dict[str, int]:
    """Validate every native row, then emit sorted distinct cross-species pairs.

    A temporary disk index bounds pair-storage memory. No output is emitted
    until parsing, ownership checks, and gzip integrity validation complete.
    Identifiers must already be exact input accessions: no heuristic remapping.
    """
    if not owners or any(not gene or any(c.isspace() for c in gene)
                         or not species for gene, species in owners.items()):
        raise ValueError("Require nonempty accession-to-species ownership")
    rows = 0
    with tempfile.TemporaryDirectory(prefix="fastoma-pairs-", dir=scratch) as tmp:
        connection = sqlite3.connect(str(Path(tmp) / "pairs.sqlite"))
        try:
            connection.execute("PRAGMA cache_size = -8192")
            connection.execute("CREATE TABLE pairs (a TEXT NOT NULL, b TEXT NOT NULL, "
                               "PRIMARY KEY (a, b)) WITHOUT ROWID")
            for a, b in iter_pairs(source, owners):
                connection.execute("INSERT OR IGNORE INTO pairs VALUES (?, ?)", (a, b))
                rows += 1
                if rows % 100000 == 0:
                    connection.commit()
            if not rows:
                raise ValueError("Empty FastOMA native pair file")
            connection.commit()
            distinct = connection.execute("SELECT COUNT(*) FROM pairs").fetchone()[0]
            for a, b in connection.execute("SELECT a, b FROM pairs ORDER BY a, b"):
                output.write(f"{a}\t{b}\n")
            return {"native_rows": rows, "distinct_pairs": distinct,
                    "duplicate_relations": rows - distinct}
        finally:
            connection.close()
