"""Resume and inventory the public Broccoli supplement without extraction."""

import argparse
import hashlib
from pathlib import Path
import subprocess
import time
import zipfile

from benchmark_tools.probe_dgx_step_separation import save

URL = "https://zenodo.org/api/records/3710751/files/data_Zenodo.zip/content"
SIZE = 719576680
MD5 = "fd77f9ef7a5b82a87143b603902901d4"


def inventory(path, expected_size=SIZE, expected_md5=MD5):
    path = Path(path)
    if path.is_symlink() or not path.is_file() or path.stat().st_size != expected_size:
        raise ValueError("Archive size or file type differs")
    md5, sha = hashlib.md5(), hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            md5.update(block)
            sha.update(block)
    if md5.hexdigest() != expected_md5:
        raise ValueError("Archive MD5 differs from public record")
    with zipfile.ZipFile(path) as archive:
        members = [dict(name=m.filename, bytes=m.file_size,
                        compressed_bytes=m.compress_size, crc32=m.CRC,
                        directory=m.is_dir()) for m in archive.infolist()]
    candidates = [m for m in members if any(token in m["name"].lower()
                  for token in ("treefam", "reference", ".nhx", "readme"))]
    return dict(status="archive_checksum_verified_and_indexed", url=URL,
                archive=str(path.absolute()), bytes=expected_size,
                md5=md5.hexdigest(), sha256=sha.hexdigest(),
                member_count=len(members), members=members, candidate_names=candidates,
                originals_recovered=False,
                limitations=["ZIP member names are leads, not verified original reference content.",
                             "No extraction or execution; member CRCs and nested archives not checked."])


def run(directory):
    directory = Path(directory).resolve(strict=True)
    archive = directory / "data_Zenodo.zip"
    if archive.is_symlink():
        raise ValueError("Refuse indirect archive path")
    size = archive.stat().st_size if archive.exists() else 0
    if size > SIZE:
        raise ValueError("Partial archive exceeds expected size")
    save(directory / "scheduled_download_started.json", dict(url=URL,
        expected_bytes=SIZE, expected_md5=MD5, resume_bytes=size,
        started_epoch=time.time()))
    with (directory / "download.log").open("x") as log:
        result = subprocess.run(["curl", "--fail", "--location", "--silent", "--show-error",
            "--connect-timeout", "30", "--max-time", "7000", "--continue-at", "-",
            "--output", str(archive), URL], stdout=log, stderr=subprocess.STDOUT,
            timeout=7030)
    save(directory / "scheduled_download_finished.json", dict(exit_code=result.returncode,
        finished_epoch=time.time(), retained_bytes=archive.stat().st_size if archive.exists() else 0))
    if result.returncode:
        raise RuntimeError("Download failed; partial bytes retained, no automatic retry")
    save(directory / "inventory.json", inventory(archive))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    run(parser.parse_args().directory)
