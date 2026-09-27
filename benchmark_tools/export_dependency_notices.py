"""Collect exact wheel-embedded notice candidates without exporting wheel code."""

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import stat
import zipfile

from benchmark_tools.inventory_dependency_notices import audit
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_publication_pipeline import save


def safe_member(name):
    path = PurePosixPath(name)
    if not path.parts or path.is_absolute() or ".." in path.parts or str(path) != name or "\\" in name:
        raise ValueError("Unsafe notice member path")
    return name


def select_wheels(inventories):
    selected = {}
    for inventory in inventories:
        if inventory["status"] != "local_wheel_notice_inventory":
            raise ValueError("Require a successful notice inventory")
        for wheel in inventory["wheels"]:
            digest = wheel["wheel"]["sha256"]
            if digest in selected:
                left = {k: v for k, v in selected[digest].items() if k != "wheel"}
                right = {k: v for k, v in wheel.items() if k != "wheel"}
                if left != right or selected[digest]["wheel"]["bytes"] != wheel["wheel"]["bytes"]:
                    raise ValueError("Conflicting inventories for identical wheel digest")
            else:
                selected[digest] = wheel
    if not selected:
        raise ValueError("Empty notice inventory")
    return [selected[k] for k in sorted(selected)]


def verify(directory):
    manifest = json.loads((directory / "NOTICE_INDEX.json").read_text())
    if manifest["scope"] != "wheel_notice_candidates_only" or manifest["redistribution_clearance"] is not False:
        raise ValueError("Wrong notice export scope")
    expected = {}
    digests = [w["wheel"]["sha256"] for w in manifest["wheels"]]
    if len(set(digests)) != len(digests):
        raise ValueError("Duplicate wheel inventory")
    for wheel in manifest["wheels"]:
        prefix = wheel["wheel"]["sha256"]
        for notice in wheel["notice_candidates"]:
            name = prefix + "/" + safe_member(notice["path"])
            if name in expected:
                raise ValueError("Duplicate notice candidate")
            expected[name] = dict(relative_path=name, wheel_sha256=prefix, member=notice["path"],
                                  bytes=notice["bytes"], sha256=notice["sha256"])
    names = set()
    for row in manifest["files"]:
        name = safe_member(row["relative_path"])
        if name in names:
            raise ValueError("Duplicate exported notice")
        if expected.get(name) != row:
            raise ValueError("Notice does not match embedded candidate inventory")
        names.add(name)
        path = directory / name
        if path.is_symlink() or not path.resolve().is_relative_to(directory.resolve()):
            raise ValueError("Indirect or escaping notice")
        actual = record(path)
        if any(actual[k] != row[k] for k in ("bytes", "sha256")):
            raise ValueError("Exported notice bytes differ")
    actual = {str(p.relative_to(directory)) for p in directory.rglob("*") if p.is_file() or p.is_symlink()}
    if names != set(expected) or actual != names | {"NOTICE_INDEX.json"}:
        raise ValueError("Unexpected notice export inventory")
    return dict(status="wheel_notice_candidates_export_verified", files=len(names),
                wheels=len(manifest["wheels"]), index=record(directory / "NOTICE_INDEX.json"),
                redistribution_clearance=False, publication_ready=False)


def export(paths, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    watched, inventories = [], []
    for path in paths:
        item = record(path)
        inventory = json.loads(path.read_text())
        report = inventory["install_report"]
        # Recompute the inventory before using its member names as extraction instructions.
        if audit(Path(report["path"]), report["sha256"]) != inventory:
            raise ValueError("Retained notice inventory no longer reproduces")
        watched.append(item)
        inventories.append(inventory)
    wheels = select_wheels(inventories)
    output.mkdir(parents=True)
    files = []
    for wheel in wheels:
        check(wheel["wheel"])
        prefix = wheel["wheel"]["sha256"]
        with zipfile.ZipFile(wheel["wheel"]["path"]) as archive:
            for notice in wheel["notice_candidates"]:
                name = safe_member(notice["path"])
                info = archive.getinfo(name)
                kind = stat.S_IFMT(info.external_attr >> 16)
                if info.is_dir() or kind not in (0, stat.S_IFREG):
                    raise ValueError("Notice is not a regular file")
                data = archive.read(name)
                if len(data) != notice["bytes"] or hashlib.sha256(data).hexdigest() != notice["sha256"]:
                    raise ValueError("Embedded notice changed")
                relative = prefix + "/" + name
                target = output / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                with target.open("xb") as stream:
                    stream.write(data)
                files.append(dict(relative_path=relative, wheel_sha256=prefix, member=name,
                                  bytes=notice["bytes"], sha256=notice["sha256"]))
        check(wheel["wheel"])
    for item in watched:
        check(item)
    save(output / "NOTICE_INDEX.json", dict(scope="wheel_notice_candidates_only", inputs=watched,
        wheels=wheels, files=files, source=record(__file__), redistribution_clearance=False,
        publication_ready=False, limitations=["Exact embedded candidate texts, not a license compatibility or clearance determination.",
            "Heuristic inventory can miss notices and does not map every static/native component to its obligations.",
            "Wheel code, binaries, external tools, OS libraries and datasets are not exported.",
            "No public deposition or redistribution was performed."]))
    return verify(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    if args.receipt.exists():
        raise FileExistsError(args.receipt)
    result = export([p.resolve() for p in args.inventory], args.output.absolute())
    save(args.receipt, result)
    print(json.dumps(result, indent=2))
