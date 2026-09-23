"""Recheck every proposed retained block and full interrupted-file identity."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import re

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen

AUDIT_SHA = "a0d024ff33544b64eafcb7ebe88ba3948ad0cffe8aabcb1e4ef2ce300cc71bc9"
PARTIAL_SHA = "41063f82a09dd6b8ef4b4bbec783d0643e09f44279011dfb2687a1d43fe76f6d"


def stat_key(stat):
    return (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)


def scan(stream, blocks, boundary, expected_digest):
    if stream.tell() != 0:
        raise ValueError("Require full-file scan from zero")
    digest = hashlib.sha256()
    end, previous, count, retained_rows, retained_blocks = 0, -1, 0, 0, 0
    final = None
    for block in blocks:
        if (final is not None or type(block["start"]) is not int or type(block["end"]) is not int
                or block["start"] != end or block["end"] <= end
                or type(block["input_ordinal_0based"]) is not int
                or block["input_ordinal_0based"] <= previous
                or type(block["rows"]) is not int or block["rows"] <= 0
                or block["reuse_authorized"] is not False
                or type(block["final_observed_query"]) is not bool
                or re.fullmatch(r"[0-9a-f]{64}", block["sha256"]) is None):
            raise ValueError("Invalid prefix block inventory")
        remaining = block["end"] - block["start"]
        local = hashlib.sha256()
        newlines, last = 0, b""
        while remaining:
            chunk = stream.read(min(remaining, 8*1024*1024))
            if not chunk or b"\0" in chunk:
                raise ValueError("Truncated or NUL-containing prefix block")
            digest.update(chunk)
            local.update(chunk)
            remaining -= len(chunk)
            newlines += chunk.count(b"\n")
            last = chunk[-1:]
        if local.hexdigest() != block["sha256"] or newlines != block["rows"] or last != b"\n":
            raise ValueError("Changed prefix block bytes or rows")
        count += 1
        end, previous = block["end"], block["input_ordinal_0based"]
        if block["final_observed_query"]:
            final = block
        else:
            retained_blocks += 1
            retained_rows += block["rows"]
    if (final is None or count != boundary["block_count"] or final["query"] != boundary["last_query"]
            or final["start"] != boundary["last_query_start"] or end != boundary["prefix_end"]):
        raise ValueError("Wrong incomplete-query boundary")
    tail = 0
    while chunk := stream.read(8*1024*1024):
        digest.update(chunk)
        tail += len(chunk)
    if tail != boundary["excluded_tail_bytes"] or digest.hexdigest() != expected_digest:
        raise ValueError("Changed complete interrupted file")
    return dict(retained_blocks=retained_blocks, retained_rows=retained_rows,
        retained_end=final["start"], excluded_final_query=final["query"],
        excluded_final_complete_rows=final["rows"], excluded_tail_bytes=tail,
        full_file_bytes=end+tail, full_file_sha256=digest.hexdigest())


def audit(root, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    path = root / "benchmarks/work/qfo_blast_prefix_audit_20260923/status.json"
    original = read_frozen(path, AUDIT_SHA)
    if (original["status"] != "observed_prefix_rows_validated_not_admitted"
            or original["reuse_authorized"] is not False or original["search_admitted"] is not False
            or original["boundary"]["sha256"] != PARTIAL_SHA):
        raise ValueError("Unexpected original prefix audit")
    inputs = [record(path), original["blocks"], original["source"], original["validator"], *original["inputs"]]
    for item in inputs:
        check(item)
    partial = root / "benchmarks/results/qfo_corrected_orthomcl_v1/work/all.blast.partial"
    with partial.open("rb") as stream, Path(original["blocks"]["path"]).open() as blocks:
        before = stat_key(os.fstat(stream.fileno()))
        result = scan(stream, (json.loads(line) for line in blocks), original["boundary"], PARTIAL_SHA)
        if before != stat_key(os.fstat(stream.fileno())) or before != stat_key(partial.stat()):
            raise ValueError("Interrupted file changed during scan")
    if result["retained_blocks"] != 885224 or result["retained_end"] != 33141800005:
        raise ValueError("Wrong production retention boundary")
    for item in inputs:
        check(item)
    report = dict(status="prefix_block_bytes_rechecked_not_admitted", source=record(__file__),
        checked_inputs=inputs, partial_path=str(partial), result=result,
        scheduler_job_id=os.environ.get("SLURM_JOB_ID"),
        search_admitted=False, reuse_authorized=False, publication_ready=False,
        limitations=["Byte identity and row-count recheck, not new numerical alignment validation.",
            "Original audit supplies structural validation; source-path and replay evidence require separate review.",
            "Does not prove pre-interruption filesystem durability or authorize whole-search admission.",
            "Entire final observed query remains excluded; absent queries still require completed replay."])
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit(args.root.resolve(), args.output.absolute())
