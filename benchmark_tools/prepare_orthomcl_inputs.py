#!/usr/bin/env python3
"""Build OrthoMCL 1.4 all.fa and all.gg files from benchmark proteomes."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


def prepare_inputs(
    input_dir: Path,
    fasta_output: Path,
    gg_output: Path,
    summary_output: Path | None = None,
) -> dict[str, object]:
    """Write native-compatible combined FASTA and genome-gene files."""
    sources = sorted(input_dir.glob("*.fasta"))
    if not sources:
        raise ValueError(f"No .fasta files found in {input_dir}")

    taxa = [source.stem for source in sources]
    if len(set(taxa)) != len(taxa):
        raise ValueError("Input FASTA names do not produce unique taxon names")

    fasta_output.parent.mkdir(parents=True, exist_ok=True)
    gg_output.parent.mkdir(parents=True, exist_ok=True)
    fasta_temporary = fasta_output.with_name(fasta_output.name + ".partial")
    gg_temporary = gg_output.with_name(gg_output.name + ".partial")
    fasta_temporary.unlink(missing_ok=True)
    gg_temporary.unlink(missing_ok=True)

    observed: set[bytes] = set()
    per_taxon: dict[str, int] = {}
    sequence_count = 0

    try:
        with fasta_temporary.open("wb") as combined, gg_temporary.open("wb") as gg:
            for source, taxon in zip(sources, taxa, strict=True):
                genes: list[bytes] = []
                current_id: bytes | None = None
                with source.open("rb") as handle:
                    for line_number, line in enumerate(handle, start=1):
                        normalized = line.replace(b"\r", b"").replace(b"\n", b"")
                        if normalized.startswith(b">"):
                            fields = normalized[1:].split(None, 1)
                            if not fields:
                                raise ValueError(
                                    f"Empty FASTA header at {source}:{line_number}"
                                )
                            current_id = fields[0]
                            if current_id in observed:
                                raise ValueError(
                                    "Duplicate FASTA ID: "
                                    + current_id.decode(errors="replace")
                                )
                            observed.add(current_id)
                            genes.append(current_id)
                            combined.write(b">" + current_id + b"\n")
                        elif current_id is None:
                            if normalized:
                                raise ValueError(
                                    f"Sequence before first header at {source}:{line_number}"
                                )
                        else:
                            combined.write(normalized + b"\n")

                if not genes:
                    raise ValueError(f"No sequences found in {source}")
                gg.write(
                    taxon.encode() + b":" + b"".join(b" " + gene for gene in genes)
                )
                gg.write(b"\n")
                per_taxon[taxon] = len(genes)
                sequence_count += len(genes)

        os.replace(fasta_temporary, fasta_output)
        os.replace(gg_temporary, gg_output)
    except BaseException:
        fasta_temporary.unlink(missing_ok=True)
        gg_temporary.unlink(missing_ok=True)
        raise

    summary: dict[str, object] = {
        "input_directory": str(input_dir.resolve()),
        "taxon_count": len(taxa),
        "sequence_count": sequence_count,
        "taxa": per_taxon,
    }
    if summary_output is not None:
        summary_output.parent.mkdir(parents=True, exist_ok=True)
        temporary = summary_output.with_name(summary_output.name + ".partial")
        temporary.write_text(json.dumps(summary, indent=2) + "\n")
        os.replace(temporary, summary_output)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("fasta_output", type=Path)
    parser.add_argument("gg_output", type=Path)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    summary = prepare_inputs(
        args.input_dir.resolve(),
        args.fasta_output.resolve(),
        args.gg_output.resolve(),
        args.summary.resolve() if args.summary else None,
    )
    print(
        f"Prepared {summary['sequence_count']:,} sequences from "
        f"{summary['taxon_count']} taxa"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
