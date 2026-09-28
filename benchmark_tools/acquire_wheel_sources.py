"""Acquire PyPI source candidates bound to exact retained wheels, without builds."""

import argparse
from datetime import datetime, timezone
from email.parser import BytesParser
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import tarfile
from urllib.parse import urlsplit
from urllib.request import urlopen

from benchmark_tools.export_dependency_notices import safe_member
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_publication_pipeline import save


def normalized(name):
    return re.sub(r"[-_.]+", "-", name).lower()


def select_source(metadata, wheel):
    if (normalized(metadata["info"]["name"]) != normalized(wheel["name"])
            or metadata["info"]["version"] != wheel["version"]):
        raise ValueError("PyPI release identity differs")
    matches = [r for r in metadata["urls"] if r["filename"] == Path(wheel["wheel"]["path"]).name]
    if len(matches) != 1 or matches[0]["packagetype"] != "bdist_wheel":
        raise ValueError("Exact wheel not uniquely published in release")
    match = matches[0]
    if match["digests"]["sha256"] != wheel["wheel"]["sha256"] or match["size"] != wheel["wheel"]["bytes"]:
        raise ValueError("Published wheel bytes differ")
    sources = [r for r in metadata["urls"] if r["packagetype"] == "sdist"]
    if len(sources) != 1:
        raise ValueError("Require exactly one source candidate")
    source = sources[0]
    safe_member(source["filename"])
    if PurePosixPath(source["filename"]).name != source["filename"] or not source["filename"].endswith(".tar.gz"):
        raise ValueError("Require plain tar.gz filename")
    url = urlsplit(source["url"])
    if url.scheme != "https" or url.netloc != "files.pythonhosted.org" or url.query or url.fragment:
        raise ValueError("Unexpected source download origin")
    if not re.fullmatch(r"[0-9a-f]{64}", source["digests"]["sha256"]):
        raise ValueError("Invalid source digest")
    if type(source["size"]) is not int or not 0 < source["size"] <= 100_000_000:
        raise ValueError("Unexpected source size")
    return source


def download(url, path, maximum):
    with urlopen(url, timeout=60) as response, path.open("xb") as output:
        if response.geturl() != url:
            raise ValueError("Unexpected redirect")
        size = 0
        while block := response.read(1024 * 1024):
            size += len(block)
            if size > maximum:
                raise ValueError("Download exceeds limit")
            output.write(block)
    return record(path)


def inspect_source(path, name, version):
    before = record(path)
    entries, notices, names, roots = [], [], set(), set()
    total = 0
    package = None
    with tarfile.open(path, "r:gz") as archive:
        for member in archive:
            member_name = safe_member(member.name)
            if member_name in names or len(names) >= 100_000:
                raise ValueError("Duplicate or excessive source members")
            names.add(member_name)
            parts = PurePosixPath(member_name).parts
            roots.add(parts[0])
            if member.isdir():
                continue
            if not member.isfile() or len(parts) < 2:
                raise ValueError("Unsupported source member")
            total += member.size
            if total > 1_000_000_000 or member.size > 100_000_000:
                raise ValueError("Source content exceeds limit")
            with archive.extractfile(member) as stream:
                data = stream.read()
            if len(data) != member.size:
                raise ValueError("Truncated source member")
            row = dict(member=member_name, bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
            entries.append(row)
            basename = parts[-1].lower()
            if basename.startswith(("license", "licence", "copying", "copyright", "notice", "authors")):
                notices.append(row)
            if len(parts) == 2 and parts[-1] == "PKG-INFO":
                package = BytesParser().parsebytes(data)
    if len(roots) != 1 or package is None:
        raise ValueError("Missing unique source root or root PKG-INFO")
    if normalized(package["Name"] or "") != normalized(name) or package["Version"] != version:
        raise ValueError("Source package identity differs")
    check(before)
    return dict(archive=before, package_name=package["Name"], package_version=package["Version"],
                files=entries, notice_candidates=notices, total_regular_bytes=total)


def acquire(inventory_path, packages, output):
    watched = record(inventory_path)
    inventory = json.loads(inventory_path.read_text())
    if inventory["status"] != "selected_wheel_elf_inventory":
        raise ValueError("Require selected wheel ELF inventory")
    requested = [normalized(p) for p in packages]
    if len(set(requested)) != len(requested) or not requested:
        raise ValueError("Require unique package names")
    wheels = []
    for name in requested:
        matches = [w for w in inventory["wheels"] if normalized(w["name"]) == name]
        if len(matches) != 1:
            raise ValueError("Package not unique in wheel inventory: " + name)
        wheels.extend(matches)
    output.mkdir(parents=True, exist_ok=False)
    rows = []
    for wheel in wheels:
        check(wheel["wheel"])
        name, version = normalized(wheel["name"]), wheel["version"]
        if not re.fullmatch(r"[a-z0-9-]+", name) or not re.fullmatch(r"[a-zA-Z0-9.+!-]+", version):
            raise ValueError("Unsupported release URL component")
        url = "https://pypi.org/pypi/" + name + "/" + version + "/json"
        directory = output / name
        directory.mkdir()
        metadata_record = download(url, directory / "pypi.json", 10_000_000)
        metadata = json.loads((directory / "pypi.json").read_text())
        source = select_source(metadata, wheel)
        artifact = download(source["url"], directory / source["filename"], source["size"])
        if artifact["bytes"] != source["size"] or artifact["sha256"] != source["digests"]["sha256"]:
            raise ValueError("Source download identity differs")
        inspection = inspect_source(Path(artifact["path"]), name, version)
        check(wheel["wheel"])
        rows.append(dict(name=name, version=version, wheel=wheel["wheel"], metadata_url=url,
                         metadata=metadata_record, source_publication=source, inspection=inspection))
    check(watched)
    return dict(status="wheel_bound_source_candidates_acquired", input=watched, source=record(__file__),
                retrieved_utc=datetime.now(timezone.utc).isoformat(), packages=rows,
                redistribution_clearance=False, publication_ready=False,
                limitations=["Exact wheel and source hashes share a PyPI release; this is not a reproducible-build attestation.",
                    "Source candidates may omit libraries bundled or statically linked into wheels.",
                    "Notice candidates are filename-based, not complete component attribution or legal clearance.",
                    "No source or wheel code was installed, imported, built or executed."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--package", action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    if args.receipt.exists():
        raise FileExistsError(args.receipt)
    save(args.receipt, acquire(args.inventory, args.package, args.output))
