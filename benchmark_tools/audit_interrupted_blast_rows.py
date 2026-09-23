"""Validate an observed BLAST prefix without authorizing reuse or editing it."""

import argparse
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_interrupted_blast_bytes import identity
from benchmark_tools.audit_orthomcl_search_table import audit_table
from benchmark_tools.convert_orthomcl_blast import read_fasta_lengths
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

PARTIAL_SHA = "41063f82a09dd6b8ef4b4bbec783d0643e09f44279011dfb2687a1d43fe76f6d"
PREFIX_END = 33141804090
INPUT_SHA = "11c03d6575e22c2f8bb718e59b7a06911637438c2e45181e86548647d49199c0"
LOG_SHA = "876f55697e53bf937ef3a0b4f5f021ab7350853752251fee814f7002ef8b0617"


class ObservedPrefix:
    def __init__(self, path, end, expected_sha, ordinals, blocks):
        self.path, self.end, self.expected_sha = path, end, expected_sha
        self.ordinals, self.blocks = ordinals, blocks
        self.summary = None

    @contextmanager
    def open(self, mode):
        if mode != "rb":
            raise ValueError("Read-only prefix")
        with self.path.open("rb") as stream:
            before = identity(os.fstat(stream.fileno()))
            yield self.lines(stream)
            if self.summary is None:
                raise ValueError("Prefix scan did not finish")
            if before != identity(os.fstat(stream.fileno())) or before != identity(self.path.stat()):
                raise ValueError("Partial input changed during scan")

    def lines(self, stream):
        digest = hashlib.sha256()
        block_hash = hashlib.sha256()
        offset = start = rows = blocks = 0
        current = None
        previous_ordinal = -1

        def emit(end, final):
            self.blocks.write(json.dumps({"query": current.decode(), "start": start, "end": end,
                "rows": rows, "sha256": block_hash.hexdigest(), "final_observed_query": final,
                "input_ordinal_0based": previous_ordinal, "reuse_authorized": False}) + "\n")

        while offset < self.end:
            line = stream.readline(min(self.end - offset, 65536))
            if not line or not line.endswith(b"\n") or b"\0" in line:
                raise ValueError("Malformed or oversized row before frozen prefix end")
            query = line.split(b"\t", 1)[0]
            if query not in self.ordinals:
                raise ValueError("Unknown query in observed prefix")
            if query != current:
                ordinal = self.ordinals[query]
                if ordinal <= previous_ordinal:
                    raise ValueError("Query blocks do not follow FASTA order")
                if current is not None:
                    emit(offset, False)
                current, previous_ordinal, start, rows = query, ordinal, offset, 0
                block_hash = hashlib.sha256()
                blocks += 1
            digest.update(line)
            block_hash.update(line)
            rows += 1
            offset += len(line)
            yield line
        if current is None:
            raise ValueError("Empty observed prefix")
        emit(offset, True)
        tail_bytes = 0
        while chunk := stream.read(8 * 1024 * 1024):
            digest.update(chunk)
            tail_bytes += len(chunk)
        if digest.hexdigest() != self.expected_sha:
            raise ValueError("Interrupted file checksum mismatch")
        self.summary = {"prefix_end": self.end, "block_count": blocks,
                        "last_query": current.decode(), "last_query_start": start,
                        "excluded_tail_bytes": tail_bytes, "sha256": digest.hexdigest()}


def run(root, output):
    base = root / "benchmarks/results/qfo_corrected_orthomcl_v1"
    fasta, log = base / "work/all.fa", base / "search_execution/blast.log"
    inputs = [record(fasta), record(log)]
    if [item["sha256"] for item in inputs] != [INPUT_SHA, LOG_SHA]:
        raise ValueError("Changed frozen FASTA or diagnostics")
    lengths = read_fasta_lengths(fasta)
    ordinals = {gene: i for i, gene in enumerate(lengths)}
    output.mkdir(parents=True, exist_ok=False)
    report = {"status": "running", "search_admitted": False, "reuse_authorized": False,
              "inputs": inputs, "source": record(__file__),
              "validator": record(Path(__file__).with_name("audit_orthomcl_search_table.py"))}
    try:
        with (output / "query_blocks.jsonl").open("x") as blocks:
            prefix = ObservedPrefix(base / "work/all.blast.partial", PREFIX_END, PARTIAL_SHA, ordinals, blocks)
            report["content"] = audit_table(prefix, fasta, log)
        for item in inputs:
            check(item)
        report.update(status="observed_prefix_rows_validated_not_admitted", boundary=prefix.summary,
                      blocks=record(output / "query_blocks.jsonl"))
        report["limitations"] = ["Final observed query remains incomplete and excluded from reuse.",
            "Absent query rows do not prove a completed no-hit search.",
            "Structural validity and order do not prove native replay equivalence or durable completion."]
    except BaseException as exc:
        report.update(status="failed_or_interrupted", error=repr(exc))
        raise
    finally:
        with (output / "status.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps({"status": run(args.root.resolve(), args.output.resolve())["status"]}))
