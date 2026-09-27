"""Preserve all local Python source pins for the completed canonical assessment."""

import argparse
import gzip
import hashlib
import io
from pathlib import Path
import subprocess
import tarfile

from benchmark_tools.run_qfo_order_replay import record, save
from benchmark_tools.run_simulation_methods import read_frozen

PLAN_SHA = "41910a6cfc113c028bbcbe0bed2bb9d8ca20c05cf8bb530d1868e31c79ca499b"
REVISION = "97a277afe4b2e28a9d1a10bbc8dd905415daccb8"


def source_inventory(repo, records):
    selected = {}
    for item in records:
        path = Path(item["path"])
        if path.suffix != ".py" or not path.is_relative_to(repo):
            continue
        relative = path.relative_to(repo)
        if ".." in relative.parts:
            raise ValueError("Source path escapes repository")
        name = "sources/" + relative.as_posix()
        if name in selected and selected[name] != item:
            raise ValueError("Conflicting source identities")
        selected[name] = item
    if not selected:
        raise ValueError("Empty source inventory")
    return selected


def matches(payload, item):
    return len(payload) == item["bytes"] and hashlib.sha256(payload).hexdigest() == item["sha256"]


def verify_archive(path, expected):
    seen = set()
    with tarfile.open(path, "r:gz") as stream:
        for member in stream:
            if not member.isfile() or member.name in seen or member.name not in expected:
                raise ValueError("Unexpected, duplicate or nonregular source member")
            if not matches(stream.extractfile(member).read(), expected[member.name]):
                raise ValueError("Archived source does not match original pin")
            seen.add(member.name)
    if seen != set(expected):
        raise ValueError("Missing archived source")
    return len(seen)


def archive(repo, output):
    if output.exists():
        raise FileExistsError(output)
    plan_path = repo / "benchmarks/work/qfo_canonical_assessment_20260927/plan.json"
    plan = read_frozen(plan_path, PLAN_SHA)
    expected = source_inventory(repo, plan["checked_records"])
    output.mkdir(parents=True)
    target = output / "qfo-canonical-assessment-python-sources.tar.gz"
    rows = []
    with target.open("xb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as zipped:
        with tarfile.open(fileobj=zipped, mode="w") as stream:
            for name, item in sorted(expected.items()):
                path = Path(item["path"])
                payload = path.read_bytes() if path.is_file() else b""
                origin = "retained_exact_bytes"
                if not matches(payload, item):
                    relative = path.relative_to(repo).as_posix()
                    payload = subprocess.check_output(["git", "show", f"{REVISION}:{relative}"], cwd=repo)
                    origin = "frozen_git_blob"
                if not matches(payload, item):
                    raise ValueError("Cannot recover exact source: " + item["path"])
                member = tarfile.TarInfo(name)
                member.size, member.mode, member.mtime = len(payload), 0o644, 0
                stream.addfile(member, io.BytesIO(payload))
                rows.append(dict(member=name, original=item, origin=origin))
    count = verify_archive(target, expected)
    read_frozen(plan_path, PLAN_SHA)
    result = dict(status="canonical_assessment_python_sources_archived", plan=record(plan_path),
        revision=REVISION, archive=record(target), members=rows, independently_read_members=count,
        uncompressed_bytes=sum(x["bytes"] for x in expected.values()), source=record(__file__),
        original_files_modified=False, publication_ready=False, redistribution_clearance=False,
        limitations=["Python source bytes only; archive modes normalized to 0644, not historical permission evidence",
            "Does not archive data, references, containers, executables, wheels or complete runtime",
            "Current-source validators must not silently substitute these paths for original receipt records",
            "No public deposition, execution reproduction or redistribution rights clearance"])
    save(output / "manifest.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    archive(args.repo.resolve(), args.output.resolve())
