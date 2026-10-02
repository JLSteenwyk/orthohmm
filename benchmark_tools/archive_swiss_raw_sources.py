"""Stream private retained raw inputs into a bounded, checksum-pinned handoff."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import tarfile

from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.relocate_swiss_raw_sources import frozen_json, restore_inputs, validate_records

MAX_BYTES = 1024 ** 3
MAX_BINDING_BYTES = 2 * 1024 ** 2
MAX_RECORDS = 2000


def direct(path):
    path = Path(path)
    if not path.is_absolute() or path.resolve() != path or path.is_symlink():
        raise ValueError("Require direct absolute handoff path")
    return path


def declared_files(document):
    if (not isinstance(document, dict)
            or set(document) != {"schema_version", "source", "entries", "redistribution_authorized"}
            or type(document["schema_version"]) is not int or document["schema_version"] != 1
            or document["redistribution_authorized"] is not False
            or not isinstance(document["entries"], list)):
        raise ValueError("Malformed or wrongly scoped binding")
    validate_records([document["source"]])
    if not 0 < len(document["entries"]) <= MAX_RECORDS:
        raise ValueError("Binding exceeds record budget")
    entries = document["entries"]
    if any(not isinstance(e, dict) or set(e) != {"original", "artifact"} for e in entries):
        raise ValueError("Malformed original record")
    originals = [e["original"] for e in entries]
    validate_records(originals)
    files = {}
    for entry, item in zip(entries, originals):
        name = "inputs/" + item["sha256"]
        if entry["artifact"] != name:
            raise ValueError("Changed content-addressed artifact path")
        expected = {k: item[k] for k in ("bytes", "sha256")}
        if name in files and files[name] != expected:
            raise ValueError("Conflicting content address")
        files[name] = expected
    if sum(v["bytes"] for v in files.values()) > MAX_BYTES:
        raise ValueError("Raw payload exceeds size budget")
    return files, originals


def wire_size(files):
    size = sum(512 + ((v["bytes"] + 511) // 512) * 512 for v in files.values()) + 1024
    return ((size + tarfile.RECORDSIZE - 1) // tarfile.RECORDSIZE) * tarfile.RECORDSIZE


def inventory(bindings, digest):
    bindings = direct(bindings)
    if bindings.stat().st_size > MAX_BINDING_BYTES:
        raise ValueError("Binding exceeds size budget")
    document = frozen_json(bindings, digest)
    files, originals = declared_files(document)
    binding_record = record(bindings)
    files["bindings.json"] = {k: binding_record[k] for k in ("bytes", "sha256")}
    if wire_size(files) > MAX_BYTES:
        raise ValueError("Tar payload and padding exceed size budget")
    restore_inputs(originals, document.get("source"), bindings, digest)
    root = bindings.parent
    observed = {p.relative_to(root).as_posix() for p in root.rglob("*")}
    if observed != set(files) | {"inputs"}:
        raise ValueError("Unexpected raw handoff filesystem member")
    return files


def archive(bindings, binding_sha, output):
    bindings, output = direct(bindings), direct(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if output.is_relative_to(bindings.parent):
        raise ValueError("Archive must not mutate its input directory")
    files = inventory(bindings, binding_sha)
    with output.open("xb") as stream, gzip.GzipFile(fileobj=stream, mode="wb", filename="", mtime=0,
                                                   compresslevel=1) as compressed:
        with tarfile.open(fileobj=compressed, mode="w|", format=tarfile.USTAR_FORMAT) as handle:
            for name, expected in sorted(files.items()):
                info = tarfile.TarInfo(name)
                info.size, info.mode, info.mtime = expected["bytes"], 0o644, 0
                with (bindings.parent / name).open("rb") as source:
                    copied = DigestReader(source)
                    handle.addfile(info, copied)
                    if copied.digest.hexdigest() != expected["sha256"]:
                        raise ValueError("Raw payload changed while copying")
    if inventory(bindings, binding_sha) != files:
        raise ValueError("Raw handoff changed during archiving")
    if output.stat().st_size > MAX_BYTES:
        raise ValueError("Compressed archive exceeds size budget")
    return dict(status="private_swiss_raw_inputs_archived", archive=record(output),
                binding_sha256=binding_sha, members=len(files),
                payload_bytes=sum(v["bytes"] for v in files.values()),
                native_annotation_admission_rerun=False, redistribution_authorized=False,
                publication_ready=False)


class DigestReader:
    def __init__(self, stream):
        self.stream, self.digest = stream, hashlib.sha256()

    def read(self, size=-1):
        data = self.stream.read(size)
        self.digest.update(data)
        return data


class BoundedReader:
    def __init__(self, stream):
        self.stream, self.total = stream, 0

    def read(self, size=-1):
        allowance = MAX_BYTES - self.total
        request = allowance + 1 if size < 0 else min(size, allowance + 1)
        data = self.stream.read(request)
        self.total += len(data)
        if self.total > MAX_BYTES:
            raise ValueError("Actual decompressed stream exceeds size budget")
        return data


def restore(archive_path, archive_sha, binding_sha, output):
    archive_path, output = direct(archive_path), direct(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    if archive_path.stat().st_size > MAX_BYTES:
        raise ValueError("Compressed archive exceeds size budget")
    archive_record = record(archive_path)
    if archive_record["sha256"] != archive_sha:
        raise ValueError("Pinned archive changed")
    output.mkdir(parents=True)
    seen, files, logical_bytes = set(), None, 0
    with archive_path.open("rb") as stream, gzip.GzipFile(fileobj=stream, mode="rb") as compressed:
        bounded = BoundedReader(compressed)
        with tarfile.open(fileobj=bounded, mode="r|") as handle:
            for member in handle:
                name = member.name
                if (not member.isfile() or member.mode != 0o644 or member.size < 0
                        or member.pax_headers or name in seen or len(seen) > MAX_RECORDS):
                    raise ValueError("Unexpected or duplicate raw archive member")
                if files is None:
                    if name != "bindings.json" or member.size > MAX_BINDING_BYTES:
                        raise ValueError("Pinned binding must be the first bounded member")
                    expected = dict(bytes=member.size, sha256=binding_sha)
                else:
                    if name not in files:
                        raise ValueError("Unexpected raw archive path")
                    expected = files[name]
                if member.size != expected["bytes"]:
                    raise ValueError("Raw archive member size changed")
                logical_bytes += member.size
                if logical_bytes > MAX_BYTES:
                    raise ValueError("Raw payload exceeds size budget")
                target = output / name
                target.parent.mkdir(parents=True, exist_ok=True)
                digest, remaining = hashlib.sha256(), member.size
                with handle.extractfile(member) as source, target.open("xb") as destination:
                    while remaining:
                        data = source.read(min(1024 ** 2, remaining))
                        if not data:
                            raise ValueError("Truncated raw archive member")
                        remaining -= len(data)
                        digest.update(data)
                        destination.write(data)
                if digest.hexdigest() != expected["sha256"]:
                    raise ValueError("Raw archive member digest changed")
                seen.add(name)
                if files is None:
                    document = frozen_json(target, binding_sha)
                    files, _ = declared_files(document)
                    if wire_size({**files, "bindings.json": expected}) > MAX_BYTES:
                        raise ValueError("Tar payload and padding exceed size budget")
        # Consume the gzip footer and charge padding/trailing streams to the same bound.
        while bounded.read(1024 ** 2):
            pass
    if files is None or seen != set(files) | {"bindings.json"}:
        raise ValueError("Incomplete raw archive inventory")
    verified = inventory(output / "bindings.json", binding_sha)
    if bounded.total != wire_size(verified):
        raise ValueError("Unexpected decompressed tar footprint")
    if record(archive_path) != archive_record:
        raise ValueError("Pinned archive changed during restoration")
    return dict(status="private_swiss_raw_inputs_restored", archive=archive_record,
                binding=record(output / "bindings.json"), members=len(verified),
                payload_bytes=sum(v["bytes"] for v in verified.values()),
                decompressed_bytes=bounded.total, native_annotation_admission_rerun=False,
                redistribution_authorized=False, publication_ready=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    writer = commands.add_parser("archive")
    writer.add_argument("--bindings", type=Path, required=True)
    writer.add_argument("--binding-sha256", required=True)
    writer.add_argument("--output", type=Path, required=True)
    reader = commands.add_parser("restore")
    reader.add_argument("--archive", type=Path, required=True)
    reader.add_argument("--archive-sha256", required=True)
    reader.add_argument("--binding-sha256", required=True)
    reader.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "archive":
        result = archive(args.bindings, args.binding_sha256, args.output)
    else:
        result = restore(args.archive, args.archive_sha256, args.binding_sha256, args.output)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
