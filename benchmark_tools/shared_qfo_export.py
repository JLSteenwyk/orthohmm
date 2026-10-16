"""Explicit, lossless removal of one known redundant native export record."""

import hashlib
from itertools import zip_longest
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

HEADER = b"gene_a\tspecies_a\tgene_b\tspecies_b\n"
CONFIDENCE_HEADER = HEADER[:-1] + b"\tconfidence\n"
CONFIDENCES = {b"high", b"medium", b"low"}
SUCCESS = "known_duplicate_removed_from_derived_exports_only"


def normalize(pairs, confidence, output_pairs, output_confidence, duplicate, rows, unique):
    """Copy literal rows; gene ownership and candidate membership need native validation."""
    pairs, confidence, output_pairs, output_confidence = map(Path,
        (pairs, confidence, output_pairs, output_confidence))
    if (type(rows) is not int or type(unique) is not int or unique <= 0 or rows != unique + 1
            or len(duplicate) != 4 or any(not isinstance(v, str) or not v for v in duplicate)):
        raise ValueError("Require one known duplicate and exact original/distinct row counts")
    paths = [p.resolve() for p in (pairs, confidence, output_pairs, output_confidence)]
    if len(set(paths)) != 4 or output_pairs.exists() or output_confidence.exists():
        raise ValueError("Require separate untouched derived export paths")
    original = {"pairs": record(pairs), "confidence": record(confidence)}
    previous = previous_row = previous_confidence = None
    offsets = {"pairs": len(HEADER), "confidence": len(CONFIDENCE_HEADER)}
    skips, total, written, repeated, repeated_confidence = {}, 0, 0, 0, None
    with pairs.open("rb") as left, confidence.open("rb") as right, \
            output_pairs.open("xb") as out_left, output_confidence.open("xb") as out_right:
        if next(left, None) != HEADER or next(right, None) != CONFIDENCE_HEADER:
            raise ValueError("Invalid pair/confidence export headers")
        out_left.write(HEADER)
        out_right.write(CONFIDENCE_HEADER)
        for line, (a, b) in enumerate(zip_longest(left, right), 2):
            if a is None or b is None or not a.endswith(b"\n") or not b.endswith(b"\n"):
                raise ValueError("Unaligned or unterminated export rows")
            fields, annotated = a[:-1].split(b"\t"), b[:-1].split(b"\t")
            if (len(fields) != 4 or len(annotated) != 5 or annotated[:4] != fields
                    or annotated[4] not in CONFIDENCES
                    or any(not field or field.strip() != field for field in fields)):
                raise ValueError("Malformed or conflicting pair/confidence row")
            decoded = tuple(field.decode("utf-8") for field in fields)
            key = decoded[0], decoded[2]
            if key[0] >= key[1] or decoded[1] == decoded[3] or (previous is not None and key < previous):
                raise ValueError("Self/reversed/same-species or unordered pair")
            total += 1
            if previous == key:
                repeated += 1
                if (repeated != 1 or decoded != tuple(duplicate) or a != previous_row
                        or b != previous_confidence):
                    raise ValueError("Unexpected duplicate identity/count or confidence conflict")
                repeated_confidence = annotated[4].decode()
                skips = {name: {"offset": offsets[name], "bytes": len(value), "line": line,
                    "sha256": hashlib.sha256(value).hexdigest()}
                    for name, value in (("pairs", a), ("confidence", b))}
            else:
                out_left.write(a)
                out_right.write(b)
                written += 1
            offsets["pairs"] += len(a)
            offsets["confidence"] += len(b)
            previous, previous_row, previous_confidence = key, a, b
    if total != rows or written != unique or repeated != 1:
        raise ValueError("Original/distinct export counts differ from the planned recovery")
    for pin in original.values():
        check(pin)
    return {"status": SUCCESS, "original": original,
        "derived": {"pairs": record(output_pairs), "confidence": record(output_confidence)},
        "original_rows": total, "distinct_rows": written, "removed_rows": repeated,
        "duplicate": list(duplicate), "duplicate_confidence": repeated_confidence, "skips": skips,
        "native_inference_rerun": False,
        "scope": "Derived serialization only; no distinct relation or confidence value changed. Original files retained."}


def _equal_bytes(left, right, length=None):
    remaining = length
    while remaining is None or remaining:
        size = 1024 * 1024 if remaining is None else min(1024 * 1024, remaining)
        a, b = left.read(size), right.read(size)
        if a != b or (remaining is not None and len(a) != size):
            raise ValueError("Derived export differs beyond the known redundant row")
        if not a:
            break
        if remaining is not None:
            remaining -= len(a)


def verify_projection(proof):
    """Verify by byte splicing, independently of the normalizer's row filter."""
    if (proof.get("status") != SUCCESS or type(proof.get("removed_rows")) is not int
            or proof.get("removed_rows") != 1 or type(proof.get("distinct_rows")) is not int
            or proof["distinct_rows"] <= 0 or type(proof.get("original_rows")) is not int
            or proof.get("original_rows") != proof.get("distinct_rows", 0) + 1
            or proof.get("native_inference_rerun") is not False
            or proof.get("duplicate_confidence", "").encode() not in CONFIDENCES):
        raise ValueError("Invalid lossless-export proof")
    pair_row = ("\t".join(proof["duplicate"]) + "\n").encode("utf-8")
    confidence_row = pair_row[:-1] + b"\t" + proof["duplicate_confidence"].encode() + b"\n"
    for name, expected in (("pairs", pair_row), ("confidence", confidence_row)):
        original, derived = proof["original"][name], proof["derived"][name]
        check(original)
        check(derived)
        skip = proof["skips"][name]
        if (Path(original["path"]).resolve() == Path(derived["path"]).resolve()
                or type(skip["offset"]) is not int or skip["offset"] < len(expected)
                or skip["bytes"] != len(expected)
                or skip["sha256"] != hashlib.sha256(expected).hexdigest()
                or original["bytes"] - derived["bytes"] != len(expected)):
            raise ValueError("Invalid byte projection bounds or duplicate identity")
        with Path(original["path"]).open("rb") as source, Path(derived["path"]).open("rb") as target:
            source.seek(skip["offset"] - len(expected))
            if source.read(len(expected)) != expected or source.read(len(expected)) != expected:
                raise ValueError("Projection does not omit an identical adjacent record")
            source.seek(0)
            _equal_bytes(source, target, skip["offset"])
            if source.read(len(expected)) != expected:
                raise ValueError("Duplicate changed during projection readback")
            _equal_bytes(source, target)
        check(original)
        check(derived)
    return {"status": "full_byte_projection_verified", "distinct_rows": proof["distinct_rows"],
            "removed_identical_rows_per_file": 1, "original_files_unchanged": True}
