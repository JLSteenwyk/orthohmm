"""Validate a direct-review tar against external anchors before extracting regular files."""

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import tarfile


INDEX = "REVIEW_INDEX.json"


def identity(path):
    digest, size = hashlib.sha256(), 0
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
            size += len(block)
    return dict(bytes=size, sha256=digest.hexdigest())


def relative(name):
    path = PurePosixPath(name)
    if (not name or path.is_absolute() or ".." in path.parts or "\\" in name
            or str(path) != name or name == "."):
        raise ValueError("Unsafe archive member path")
    return name


def restore(archive, archive_sha256, index_sha256, output):
    archive, output = Path(archive).resolve(strict=True), Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if identity(archive)["sha256"] != archive_sha256:
        raise ValueError("Archive differs from external anchor")
    with tarfile.open(archive, "r:gz") as container:
        members = container.getmembers()
        names = [relative(m.name) for m in members]
        if (len(names) != len(set(names)) or INDEX not in names
                or any(not m.isfile() for m in members)):
            raise ValueError("Unexpected member type, missing index or duplicate")
        index_member = members[names.index(INDEX)]
        if index_member.size > 1024 * 1024 or index_member.mode not in (0o644, 0o664):
            raise ValueError("Unexpected index size or mode")
        raw_index = container.extractfile(index_member).read()
        if hashlib.sha256(raw_index).hexdigest() != index_sha256:
            raise ValueError("Index differs from external anchor")
        index = json.loads(raw_index)
        if (index["schema"] not in ("publication_direct_review_v1", "publication_direct_review_v2",
                                    "publication_direct_review_v3")
                or any(index[k] is not False for k in (
                    "publication_ready", "redistribution_clearance", "transitive_evidence_included"))):
            raise ValueError("Unexpected direct-review scope")
        rows = {relative(r["path"]): r for r in index["files"]}
        if INDEX in rows or len(rows) != len(index["files"]) or set(names) != set(rows) | {INDEX}:
            raise ValueError("Archive inventory differs from anchored index")
        # No extraction or copied-code execution occurs before all payload checks.
        for member in members:
            if member.name == INDEX:
                continue
            row = rows[member.name]
            if (row["mode"] not in (0o644, 0o755) or member.mode != row["mode"]
                    or member.size != row["bytes"]):
                raise ValueError("Payload mode or size differs")
            digest = hashlib.sha256()
            with container.extractfile(member) as stream:
                for block in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(block)
            if digest.hexdigest() != row["sha256"]:
                raise ValueError("Payload checksum differs")
        output.mkdir(parents=True, exist_ok=False)
        for member in members:
            path = output / member.name
            path.parent.mkdir(parents=True, exist_ok=True)
            with container.extractfile(member) as source, path.open("xb") as target:
                shutil.copyfileobj(source, target)
            path.chmod(member.mode)
    return dict(schema="direct_review_archive_restoration_v1", status="anchored_payloads_restored",
        archive=dict(path=str(archive), **identity(archive)), index_sha256=index_sha256,
        destination=str(output), members=len(members), payloads=len(rows),
        copied_code_executed=False, publication_ready=False,
        limitations=["Restoration checks direct payloads, not transitive evidence or redistribution rights.",
                     "Copied verifier must be byte-checked and actually executed separately."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--archive-sha256", required=True)
    parser.add_argument("--index-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--receipt", required=True, type=Path)
    args = parser.parse_args()
    if args.receipt.exists() or args.receipt.is_symlink():
        raise FileExistsError(args.receipt)
    result = restore(args.archive, args.archive_sha256, args.index_sha256, args.output)
    with args.receipt.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps(result))
