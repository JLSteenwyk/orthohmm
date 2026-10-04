"""Build, archive, restore and verify a pinned local publication candidate."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import tarfile


INDEX = "PACKAGE_INDEX.json"
SELECTION = "PACKAGE_SELECTION.json"
RUNNER = "bundle_publication_package.py"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def relative(value):
    require(isinstance(value, str), "Require relative path string")
    path = PurePosixPath(value)
    require(path.parts and not path.is_absolute() and ".." not in path.parts
            and str(path) == value and "\\" not in value, "Unsafe package path")
    return value


def identity(path):
    digest, size = hashlib.sha256(), 0
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
            size += len(block)
    return {"bytes": size, "sha256": digest.hexdigest()}


def pin(row):
    require(type(row["bytes"]) is int and row["bytes"] >= 0
            and isinstance(row["sha256"], str) and re.fullmatch(r"[0-9a-f]{64}", row["sha256"]),
            "Invalid package identity")
    return {key: row[key] for key in ("bytes", "sha256")}


def save(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write("\n")
    Path(path).chmod(0o644)


def fresh(path):
    path = Path(path).absolute()
    if path.exists() or path.is_symlink():
        raise FileExistsError("Refusing existing output: " + str(path))
    return path


def selected(value):
    require(value["schema"] == "publication_package_selection_v1", "Unknown selection schema")
    require(isinstance(value["version"], str)
            and re.fullmatch(r"orthohmm-study-[0-9]{4}\.[0-9]{2}\.[0-9]{2}-rc[1-9][0-9]*", value["version"]),
            "Require explicit study-candidate version")
    require(value["publication_ready"] is False and value["public_archive_uploaded"] is False,
            "Candidate cannot certify publication readiness or deposition")
    require(all(isinstance(value[key], str) and re.fullmatch(r"[0-9a-f]{40}", value[key])
                for key in ("scientific_revision", "workflow_revision")), "Require exact source revisions")
    files = {}
    for row in value["files"]:
        target, source = relative(row["target"]), relative(row["source"])
        require(target not in files and target not in {INDEX, SELECTION, RUNNER},
                "Duplicate or reserved package target")
        files[target] = {"source": source, **pin(row)}
    require(files and "README.md" in files, "Missing package instructions")
    require(all(not (set(parent.as_posix() for parent in PurePosixPath(name).parents)
                     & (set(files) | {INDEX, SELECTION, RUNNER})) for name in files),
            "File/directory target collision")
    return files


def verify(directory, manifest_sha256):
    directory = Path(directory).resolve(strict=True)
    index_path = directory / INDEX
    require(not index_path.is_symlink() and identity(index_path)["sha256"] == manifest_sha256,
            "Package index differs from external anchor")
    index = json.loads(index_path.read_text())
    require(index["schema"] == "publication_package_v1" and index["publication_ready"] is False
            and index["public_archive_uploaded"] is False and index["native_inference_repeated"] is False,
            "Package scope differs")
    rows = {}
    for row in index["files"]:
        name = relative(row["path"])
        path = directory / name
        require(name != INDEX and name not in rows and not path.is_symlink() and path.is_file()
                and path.resolve().is_relative_to(directory) and row["mode"] == 0o644
                and path.stat().st_mode & 0o777 == 0o644, "Invalid package inventory")
        require(identity(path) == pin(row), "Package payload differs: " + name)
        rows[name] = row
    actual = {path.relative_to(directory).as_posix() for path in directory.rglob("*")
              if path.is_file() or path.is_symlink()}
    require(actual == set(rows) | {INDEX}, "Extra or missing package payloads")
    require(SELECTION in rows and RUNNER in rows, "Missing package support")
    chosen = json.loads((directory / SELECTION).read_text())
    sources = selected(chosen)
    require(index["version"] == chosen["version"] and index["scientific_revision"] == chosen["scientific_revision"]
            and index["workflow_revision"] == chosen["workflow_revision"]
            and set(rows) == set(sources) | {SELECTION, RUNNER}, "Selection mapping differs")
    for name, row in sources.items():
        require(pin(rows[name]) == pin(row), "Selected payload identity differs")
    return {"status": "local_publication_package_verified", "version": index["version"],
            "files": len(rows), "payload_bytes": sum(row["bytes"] for row in rows.values()),
            "manifest": identity(index_path), "publication_ready": False,
            "public_archive_uploaded": False, "native_inference_repeated": False,
            "nested_component_execution_repeated": False}


def build(root, selection, selection_sha256, output):
    root, selection = Path(root).resolve(strict=True), Path(selection).resolve(strict=True)
    output = fresh(output)
    require(identity(selection)["sha256"] == selection_sha256, "Selection checksum mismatch")
    chosen = json.loads(selection.read_text())
    sources = selected(chosen)
    for row in sources.values():
        path = root / row["source"]
        require(path.resolve().is_relative_to(root) and path.is_file(), "Source escapes root or is absent")
        require(identity(path) == pin(row), "Source checksum mismatch: " + row["source"])
    output.mkdir(parents=True)
    for name, row in sources.items():
        path = output / name
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(root / row["source"], path)
        path.chmod(0o644)
    for path, destination in ((selection, output / SELECTION), (Path(__file__), output / RUNNER)):
        shutil.copyfile(path, destination)
        destination.chmod(0o644)
    rows = [{"path": path.relative_to(output).as_posix(), "mode": 0o644, **identity(path)}
            for path in sorted(output.rglob("*")) if path.is_file()]
    save(output / INDEX, {"schema": "publication_package_v1", "version": chosen["version"],
        "scientific_revision": chosen["scientific_revision"], "workflow_revision": chosen["workflow_revision"],
        "files": rows, "publication_ready": False, "public_archive_uploaded": False,
        "native_inference_repeated": False, "limitations": chosen["limitations"]})
    return verify(output, identity(output / INDEX)["sha256"])


def archive(directory, manifest_sha256, output):
    directory, output = Path(directory).resolve(strict=True), fresh(output)
    result = verify(directory, manifest_sha256)
    with output.open("xb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as compressed:
        with tarfile.open(fileobj=compressed, mode="w") as container:
            for path in sorted(directory.rglob("*")):
                if not path.is_file():
                    continue
                member = container.gettarinfo(str(path), path.relative_to(directory).as_posix())
                member.uid = member.gid = member.mtime = 0
                member.uname = member.gname = ""
                member.mode = 0o644
                with path.open("rb") as stream:
                    container.addfile(member, stream)
    return {"archive": {"path": str(output), **identity(output)}, "verified": result}


def restore(archive_path, archive_sha256, manifest_sha256, output):
    output = fresh(output)
    require(identity(archive_path)["sha256"] == archive_sha256, "Archive differs from external anchor")
    with tarfile.open(archive_path, "r:gz") as container:
        members = container.getmembers()
        names = [relative(member.name) for member in members]
        require(len(names) == len(set(names)) and all(member.isfile() and member.mode == 0o644 for member in members),
                "Unexpected archive member type, mode or duplicate")
        require(INDEX in names, "Missing archive index")
        contents = container.extractfile(INDEX).read()
        require(hashlib.sha256(contents).hexdigest() == manifest_sha256, "Archive index anchor differs")
        index = json.loads(contents)
        rows = {relative(row["path"]): row for row in index["files"]}
        require(len(rows) == len(index["files"]) and set(names) == set(rows) | {INDEX}, "Archive inventory differs")
        require(INDEX not in rows and SELECTION in rows and RUNNER in rows,
                "Invalid archive support inventory")
        chosen = json.loads(container.extractfile(SELECTION).read())
        sources = selected(chosen)
        require(set(rows) == set(sources) | {SELECTION, RUNNER}, "Archive selection mapping differs")
        # Validate payload streams before any extraction or execution.
        for member in members:
            if member.name == INDEX:
                continue
            digest, size = hashlib.sha256(), 0
            with container.extractfile(member) as stream:
                for block in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(block)
                    size += len(block)
            require({"bytes": size, "sha256": digest.hexdigest()} == pin(rows[member.name]),
                    "Archive member checksum mismatch")
        output.mkdir(parents=True)
        for member in members:
            path = output / member.name
            path.parent.mkdir(parents=True, exist_ok=True)
            with container.extractfile(member) as source, path.open("xb") as target:
                shutil.copyfileobj(source, target)
            path.chmod(0o644)
    return verify(output, manifest_sha256)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    builder = commands.add_parser("build")
    builder.add_argument("--root", type=Path, required=True)
    builder.add_argument("--selection", type=Path, required=True)
    builder.add_argument("--selection-sha256", required=True)
    builder.add_argument("--output", type=Path, required=True)
    verifier = commands.add_parser("verify")
    verifier.add_argument("directory", type=Path)
    verifier.add_argument("--manifest-sha256", required=True)
    archiver = commands.add_parser("archive")
    archiver.add_argument("directory", type=Path)
    archiver.add_argument("--manifest-sha256", required=True)
    archiver.add_argument("--output", type=Path, required=True)
    restorer = commands.add_parser("restore")
    restorer.add_argument("archive", type=Path)
    restorer.add_argument("--archive-sha256", required=True)
    restorer.add_argument("--manifest-sha256", required=True)
    restorer.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "build":
        result = build(args.root, args.selection, args.selection_sha256, args.output)
    elif args.command == "verify":
        result = verify(args.directory, args.manifest_sha256)
    elif args.command == "archive":
        result = archive(args.directory, args.manifest_sha256, args.output)
    else:
        result = restore(args.archive, args.archive_sha256, args.manifest_sha256, args.output)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
