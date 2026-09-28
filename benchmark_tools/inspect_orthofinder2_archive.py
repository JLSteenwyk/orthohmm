"""Retrieve a public supplement and inventory it without extracting or executing files."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile
import time

from benchmark_tools.probe_dgx_step_separation import save

URL = "https://zenodo.org/records/1481147/files/OrthoFinder2_Zenodo.tar.gz?download=1"
SIZE = 1957180078
MD5 = "5907cdf22b211f2c4e6bead42eb17878"


def inventory(path, expected_size=SIZE, expected_md5=MD5):
    path = Path(path)
    if path.is_symlink() or not path.is_file() or path.stat().st_size != expected_size:
        raise ValueError("Archive size or file type differs")
    md5, sha = hashlib.md5(), hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024*1024), b""):
            md5.update(block)
            sha.update(block)
    if md5.hexdigest() != expected_md5:
        raise ValueError("Archive MD5 differs from public record")
    members, candidates = [], []
    with tarfile.open(path, "r|gz") as archive:
        for member in archive:
            row = dict(name=member.name, size=member.size,
                       kind="file" if member.isfile() else "directory" if member.isdir() else "other")
            members.append(row)
            name = member.name.lower()
            if any(token in name for token in ("treefam", "reference", ".nhx", "readme")):
                candidates.append(row)
    return dict(status="archive_checksum_verified_and_indexed", archive=str(path.absolute()),
                bytes=expected_size, md5=md5.hexdigest(), sha256=sha.hexdigest(),
                member_count=len(members), members=members, candidate_names=candidates,
                originals_recovered=False,
                limitations=["Candidate paths alone do not identify original release-7 trees or mapping.",
                             "No archive member extracted or executed; nested archives not yet inspected."])


def run(directory):
    directory = Path(directory).resolve(strict=True)
    archive = directory / "OrthoFinder2_Zenodo.tar.gz"
    if archive.is_symlink():
        raise ValueError("Refuse indirect archive path")
    start = dict(url=URL, expected_bytes=SIZE, expected_md5=MD5,
                 resume_bytes=archive.stat().st_size if archive.exists() else 0,
                 started_epoch=time.time())
    save(directory / "scheduled_download_started.json", start)
    with (directory / "download.log").open("x") as log:
        result = subprocess.run(["curl", "--fail", "--location", "--silent", "--show-error",
            "--connect-timeout", "30", "--max-time", "7000", "--continue-at", "-",
            "--output", str(archive), URL], stdout=log, stderr=subprocess.STDOUT, timeout=7030)
    save(directory / "scheduled_download_finished.json", dict(exit_code=result.returncode,
        finished_epoch=time.time(), retained_bytes=archive.stat().st_size if archive.exists() else 0))
    if result.returncode:
        raise RuntimeError("Download failed; partial bytes and log retained, no automatic retry")
    save(directory / "inventory.json", inventory(archive))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", required=True, type=Path)
    run(parser.parse_args().directory)
