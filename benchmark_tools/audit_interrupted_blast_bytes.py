"""Read-only byte integrity inventory; never certify partial BLAST queries."""

import argparse
import hashlib
import json
import os
from pathlib import Path


def scan(stream, chunk_size=8 * 1024 * 1024):
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")
    digest = hashlib.sha256()
    size = zeros = newlines = 0
    first_zero = last_zero = last_newline = None
    last_byte = None
    while chunk := stream.read(chunk_size):
        digest.update(chunk)
        count = chunk.count(b"\0")
        zeros += count
        if count:
            if first_zero is None:
                first_zero = size + chunk.find(b"\0")
            last_zero = size + chunk.rfind(b"\0")
        newlines += chunk.count(b"\n")
        if b"\n" in chunk:
            last_newline = size + chunk.rfind(b"\n")
        last_byte = chunk[-1]
        size += len(chunk)
    return {"bytes": size, "sha256": digest.hexdigest(), "nul_bytes": zeros,
            "first_nul_offset": first_zero, "last_nul_offset": last_zero,
            "nuls_form_single_trailing_run": bool(zeros and last_zero == size - 1
                and last_zero - first_zero + 1 == zeros),
            "newline_count": newlines, "last_newline_offset": last_newline,
            "ends_with_newline": last_byte == 10}


def identity(st):
    return (st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns)


def audit(path, output):
    if output.exists():
        raise FileExistsError(output)
    with path.open("rb") as stream:
        before = identity(os.fstat(stream.fileno()))
        result = scan(stream)
        if before != identity(os.fstat(stream.fileno())) or before != identity(path.stat()):
            raise ValueError("Input changed during audit")
    report = {"status": "interrupted_blast_byte_inventory_only", "input": str(path),
              "content": result, "search_admitted": False, "reuse_authorized": False,
              "limitations": ["No BLAST row, query completion, or ordering validation.",
                              "No-hit queries cannot be recovered from hit rows alone.",
                              "NUL-free bytes do not establish valid or complete results."]}
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.input.resolve(), args.output.resolve())))
