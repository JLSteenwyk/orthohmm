#!/usr/bin/env python3
"""Convert compact BLAST output to OrthoMCL 1.4's BPO format."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import BinaryIO


def read_fasta_lengths(path: Path) -> dict[bytes, int]:
    """Return sequence lengths keyed by the first FASTA-header token."""
    lengths: dict[bytes, int] = {}
    sequence_id: bytes | None = None
    with path.open("rb") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(b">"):
                fields = line[1:].split(None, 1)
                if not fields:
                    raise ValueError(f"Empty FASTA header at {path}:{line_number}")
                sequence_id = fields[0]
                if sequence_id in lengths:
                    raise ValueError(f"Duplicate FASTA ID: {sequence_id.decode()}")
                lengths[sequence_id] = 0
            elif sequence_id is None:
                if line.strip():
                    raise ValueError(f"Sequence before first header at {path}:{line_number}")
            else:
                lengths[sequence_id] += len(line.strip())
    return lengths


def _normalize_evalue(value: bytes) -> bytes:
    if value.startswith(b"e-"):
        return b"1" + value
    return value


def _numeric_evalue(value: bytes) -> float:
    return float(_normalize_evalue(value))


def _write_hit(
    output: BinaryIO,
    similarity_id: int,
    query: bytes,
    query_length: int,
    subject: bytes,
    subject_length: int,
    evalue: bytes,
    weighted_identity: float,
    aligned_subject_length: int,
    spans: list[bytes],
) -> None:
    percent_identity = int(weighted_identity / aligned_subject_length)
    output.write(
        b";".join(
            (
                str(similarity_id).encode(),
                query,
                str(query_length).encode(),
                subject,
                str(subject_length).encode(),
                evalue,
                str(percent_identity).encode(),
                b".".join(spans) + b".",
            )
        )
        + b"\n"
    )


def convert_blast(
    blast_path: Path,
    fasta_path: Path,
    output_path: Path,
    evalue_cutoff: float = 1e-5,
    progress_every: int = 1_000_000,
) -> int:
    """Stream BLAST m8 hits into the byte-compatible OrthoMCL BPO layout."""
    lengths = read_fasta_lengths(fasta_path)
    similarity_id = 0
    current_key: tuple[bytes, bytes] | None = None
    current_evalue = b""
    weighted_identity = 0.0
    aligned_subject_length = 0
    spans: list[bytes] = []

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with blast_path.open("rb") as blast, output_path.open("wb", buffering=8 << 20) as out:
        for line_number, line in enumerate(blast, start=1):
            fields = line.rstrip(b"\r\n").split(b"\t")
            if len(fields) != 12:
                raise ValueError(
                    f"Expected 12 BLAST columns at {blast_path}:{line_number}; "
                    f"found {len(fields)}"
                )
            query, subject = fields[0], fields[1]
            key = (query, subject)
            if key != current_key:
                if current_key is not None and _numeric_evalue(current_evalue) <= evalue_cutoff:
                    similarity_id += 1
                    _write_hit(
                        out,
                        similarity_id,
                        current_key[0],
                        lengths[current_key[0]],
                        current_key[1],
                        lengths[current_key[1]],
                        current_evalue,
                        weighted_identity,
                        aligned_subject_length,
                        spans,
                    )
                    if progress_every and similarity_id % progress_every == 0:
                        print(f"Wrote {similarity_id:,} BPO hits", file=sys.stderr)
                current_key = key
                current_evalue = _normalize_evalue(fields[10])
                weighted_identity = 0.0
                aligned_subject_length = 0
                spans = []

            alignment_length = int(fields[3])
            mismatches = int(fields[4])
            query_span = abs(int(fields[7]) - int(fields[6])) + 1
            subject_span = abs(int(fields[9]) - int(fields[8])) + 1
            total_gaps = (alignment_length - query_span) + (
                alignment_length - subject_span
            )
            identical = alignment_length - mismatches - total_gaps
            percent_identity = identical / alignment_length * 100
            weighted_identity += percent_identity * subject_span
            aligned_subject_length += subject_span
            span_id = len(spans) + 1
            spans.append(
                b":".join(
                    (
                        str(span_id).encode(),
                        fields[6] + b"-" + fields[7],
                        fields[8] + b"-" + fields[9],
                    )
                )
            )

        if current_key is not None and _numeric_evalue(current_evalue) <= evalue_cutoff:
            similarity_id += 1
            _write_hit(
                out,
                similarity_id,
                current_key[0],
                lengths[current_key[0]],
                current_key[1],
                lengths[current_key[1]],
                current_evalue,
                weighted_identity,
                aligned_subject_length,
                spans,
            )

    return similarity_id


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("blast", type=Path, help="BLAST -m 8 output")
    parser.add_argument("fasta", type=Path, help="Combined FASTA searched by BLAST")
    parser.add_argument("output", type=Path, help="Destination BPO file")
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--progress-every", type=int, default=1_000_000)
    args = parser.parse_args()
    count = convert_blast(
        args.blast,
        args.fasta,
        args.output,
        args.evalue,
        args.progress_every,
    )
    print(f"Wrote {count:,} BPO hits to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
