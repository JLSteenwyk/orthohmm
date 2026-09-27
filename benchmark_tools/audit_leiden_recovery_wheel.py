"""Verify wheel payload identity with the validated private Leiden snapshot."""

import argparse
import base64
import csv
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import zipfile

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

OVERLAY_SHA = "3cc5b1795d1d0b361a02e4d18fc98e0a1bc5307ab6c2c086fe28e7bc0e0c4e4a"
ROOTS = {"leidenalg", "leidenalg.libs", "leidenalg-0.11.0.dist-info"}
GENERATED = {"leidenalg-0.11.0.dist-info/" + n for n in ("INSTALLER", "REQUESTED", "RECORD")}


def payload(wheel):
    with zipfile.ZipFile(wheel) as archive:
        infos = archive.infolist()
        names = [i.filename for i in infos]
        if len(names) != len(set(names)):
            raise ValueError("Duplicate wheel members")
        files = {}
        for info in infos:
            path = PurePosixPath(info.filename)
            if (path.is_absolute() or ".." in path.parts or not path.parts
                    or path.parts[0] not in ROOTS or "\\" in info.filename
                    or (info.external_attr >> 16) & 0o170000 == 0o120000):
                raise ValueError("Unsafe or unexpected wheel member")
            if not info.is_dir():
                files[info.filename] = archive.read(info)
    name = "leidenalg-0.11.0.dist-info/RECORD"
    rows = list(csv.reader(io.StringIO(files[name].decode())))
    if any(len(row) != 3 for row in rows) or len(rows) != len({r[0] for r in rows}):
        raise ValueError("Malformed or duplicate RECORD rows")
    if {r[0] for r in rows} != files.keys():
        raise ValueError("RECORD inventory differs")
    for member, digest, size in rows:
        data = files[member]
        if member == name:
            if digest or size:
                raise ValueError("RECORD must not hash itself")
            continue
        expected = "sha256=" + base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b"=").decode()
        if digest != expected or size != str(len(data)):
            raise ValueError("Wheel RECORD hash or size differs")
    return files


def audit(repo, wheel):
    report_path = repo / "benchmark_tools/results/ob_leiden_overlay_probe_20260926.json"
    report_record = record(report_path)
    if report_record["sha256"] != OVERLAY_SHA:
        raise ValueError("Changed validated overlay evidence")
    report = json.loads(report_path.read_text())
    base = Path(report["plan"]["path"]).parent / "distribution"
    members = payload(wheel)
    rows, expected = [], set()
    for item in report["copied"]:
        copied = item["copy"]
        check(copied)
        name = Path(copied["path"]).relative_to(base).as_posix()
        if name in GENERATED:
            rows.append(dict(member=name, classification="installer_generated_metadata", original=copied))
            continue
        expected.add(name)
        if name not in members or hashlib.sha256(members[name]).hexdigest() != copied["sha256"] or len(members[name]) != copied["bytes"]:
            raise ValueError("Recovery wheel payload differs")
        rows.append(dict(member=name, classification="byte_identical_payload", original=copied))
    if set(members) - GENERATED != expected:
        raise ValueError("Unexpected or missing recovery payload")
    return dict(status="wheel_payload_matches_validated_leiden011_snapshot", wheel=record(wheel),
        overlay_report=report_record, members=rows, payload_files=len(expected),
        internal_record_verified=True, installed=False, source=record(__file__), publication_ready=False,
        limitations=["Payload identity, not equality of installer-generated metadata or an executed clean installation.",
            "Platform-specific historical runtime artifact; no current security or general portability claim.",
            "Other runtime dependencies, tools and production ordering-policy integration remain separate."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.repo.resolve(), args.wheel.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
