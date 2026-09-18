"""Compare selected retained FASTAs to archive members without extracting files."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import sys
import tarfile

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def verify_members(archive, expected):
    observed = {}
    member_count = 0
    with gzip.open(archive, "rb") as compressed:
        with tarfile.open(fileobj=compressed, mode="r|") as stream:
            for member in stream:
                member_count += 1
                if member.name not in expected:
                    continue
                if member.name in observed or not member.isfile():
                    raise ValueError("Duplicate or nonregular selected archive member")
                digest, size = hashlib.sha256(), 0
                with stream.extractfile(member) as handle:
                    for block in iter(lambda: handle.read(1024 * 1024), b""):
                        digest.update(block)
                        size += len(block)
                actual = {"bytes": size, "sha256": digest.hexdigest()}
                if size != member.size or actual != {k: expected[member.name][k] for k in actual}:
                    raise ValueError("Selected archive member differs from retained file")
                observed[member.name] = actual
        # Consume through gzip EOF so a corrupt checksum/trailer cannot be ignored
        # merely because tar's end marker occurred earlier.
        for _ in iter(lambda: compressed.read(1024 * 1024), b""):
            pass
    if set(observed) != set(expected):
        raise ValueError("Missing selected archive member")
    return {"members": observed, "tar_members_visited": member_count,
            "gzip_read_to_eof": True, "selected_member_bytes_verified": True}


def audit(archive, canonical, additional):
    before = record(archive)
    expected = {"Eukaryota/UP000008143_8364.fasta": record(canonical),
                "Eukaryota/UP000008143_8364_additional.fasta": record(additional)}
    result = verify_members(archive, expected)
    for identity in [before, *expected.values()]:
        check(identity)
    return {"status": "selected_qfo_archive_members_verified", **result,
            "archive": before, "retained_files": expected, "source": record(__file__),
            "limitations": ["Establishes selected retained files match this local archive, not an independently authenticated publisher checksum.",
                            "Gzip EOF/checksum and tar iteration verified; no general validation of all biological records or archive contents.",
                            "No files extracted, overwritten, remapped or used for new benchmark scoring."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("archive", "canonical", "additional", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.archive, args.canonical, args.additional)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
