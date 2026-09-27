"""Stage verified base archives for offline Conda installation; no installation."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
from urllib.parse import urlsplit, unquote


def identity(path):
    sha, md5 = hashlib.sha256(), hashlib.md5()
    size = 0
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            sha.update(block)
            md5.update(block)
            size += len(block)
    return dict(bytes=size, sha256=sha.hexdigest(), md5=md5.hexdigest())


def packages_from_receipt(data):
    packages = data["acquisition"]["packages"]
    if not isinstance(packages, list) or not packages:
        raise ValueError("Empty or invalid package list")
    names, filenames, result = set(), set(), []
    for entry in packages:
        package, archive = entry["package"], entry["archive"]
        name, filename = package["name"], package["fn"]
        if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", name):
            raise ValueError("Invalid package name")
        if (not isinstance(filename, str)
                or not re.fullmatch(r"[A-Za-z0-9_.+-]+\.(?:conda|tar\.bz2)", filename)):
            raise ValueError("Invalid archive filename")
        if name in names or filename in filenames:
            raise ValueError("Duplicate package name or filename")
        names.add(name)
        filenames.add(filename)
        sha, md5 = package["sha256"], package["md5"]
        if (not isinstance(sha, str) or not re.fullmatch(r"[0-9a-f]{64}", sha)
                or not isinstance(md5, str) or not re.fullmatch(r"[0-9a-f]{32}", md5)):
            raise ValueError("Invalid archive digest")
        if (type(archive["bytes"]) is not int or archive["bytes"] <= 0
                or archive["sha256"] != sha):
            raise ValueError("Inconsistent archive identity")
        url = urlsplit(package["url"])
        if (url.scheme != "https" or url.username or url.password
                or url.hostname not in {"repo.anaconda.com", "conda.anaconda.org"}
                or url.port not in {None, 443} or url.query or url.fragment
                or unquote(url.path.rsplit("/", 1)[-1]) != filename):
            raise ValueError("Invalid retained provider URL")
        if package["subdir"] not in {"linux-64", "noarch"}:
            raise ValueError("Unsupported package platform")
        result.append(dict(package=package, expected=dict(
            bytes=archive["bytes"], sha256=sha, md5=md5)))
    return result


def stage(receipt, cache, output):
    receipt, cache, output = Path(receipt), Path(cache), Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    source = receipt.read_bytes()
    packages = packages_from_receipt(json.loads(source))
    # Preflight the entire cache before creating any output; never trust old paths.
    for item in packages:
        path = cache / item["package"]["fn"]
        if path.is_symlink() or not path.is_file():
            raise ValueError("Missing or symlinked archive: " + str(path))
        if identity(path) != item["expected"]:
            raise ValueError("Changed archive: " + str(path))
    output.mkdir(parents=True, exist_ok=False)
    archives = output / "archives"
    archives.mkdir()
    records, lines = [], ["@EXPLICIT"]
    # A failed copy leaves an unadmitted directory for inspection, never a retry.
    for item in packages:
        package = item["package"]
        target = archives / package["fn"]
        with (cache / package["fn"]).open("rb") as src, target.open("xb") as dst:
            shutil.copyfileobj(src, dst, length=1024 * 1024)
        actual = identity(target)
        if actual != item["expected"]:
            raise ValueError("Copied archive differs: " + str(target))
        records.append(dict(package=package, archive=dict(
            path=str(target.relative_to(output)), **actual)))
        lines.append(target.as_uri() + "#" + actual["md5"])
    if receipt.read_bytes() != source:
        raise ValueError("Receipt changed during staging")
    explicit = output / "explicit.txt"
    with explicit.open("x") as stream:
        stream.write("\n".join(lines) + "\n")
    result = dict(status="base_archives_staged", installation_performed=False,
        publication_ready=False, source_receipt_sha256=hashlib.sha256(source).hexdigest(),
        packages=records, explicit=dict(path="explicit.txt", **identity(explicit)),
        limitations=["Conda installation and pip-wheel overlay are separate steps.",
            "Explicit file contains destination-specific file URLs; restage after relocation.",
            "Digests establish retained identity, not provider authenticity or security.",
            "No dependency solving, rights clearance or scientific validation is performed."])
    with (output / "staging.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("receipt", "cache", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    stage(args.receipt, args.cache, args.output)
