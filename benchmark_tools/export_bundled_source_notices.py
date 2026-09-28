"""Export the four recovered igraph component notices, not source code or binaries."""

import argparse
import hashlib
import json
from pathlib import Path
import tarfile

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_publication_pipeline import save


SELECTION = {
    "libxml2-2.9.7/Copyright": "libxml2/Copyright",
    "xz-5.2.4/COPYING": "xz/COPYING",
    "gcc-8.5.0-20210514/COPYING3": "gcc/COPYING3",
    "gcc-8.5.0-20210514/COPYING.RUNTIME": "gcc/COPYING.RUNTIME",
}


def verify(directory, index_sha256):
    directory = Path(directory)
    index = directory / "SOURCE_NOTICE_INDEX.json"
    if directory.is_symlink() or index.is_symlink() or record(index)["sha256"] != index_sha256:
        raise ValueError("Source notice index identity differs")
    data = json.loads(index.read_text())
    if (data.get("scope") != "selected_bundled_source_notices"
            or data.get("redistribution_clearance") is not False):
        raise ValueError("Wrong source notice scope")
    seen = set()
    for row in data["files"]:
        member, relative = row["member"], row["relative_path"]
        if member not in SELECTION or relative != SELECTION[member] or member in seen:
            raise ValueError("Unexpected or duplicate source notice")
        seen.add(member)
        path = directory / relative
        if any(p.is_symlink() for p in (path, path.parent)):
            raise ValueError("Indirect notice path")
        actual = record(path)
        if any(actual[k] != row[k] for k in ("bytes", "sha256")):
            raise ValueError("Source notice payload differs")
    actual = {str(p.relative_to(directory)) for p in directory.rglob("*")
              if p.is_file() or p.is_symlink()}
    if seen != set(SELECTION) or actual != set(SELECTION.values()) | {index.name}:
        raise ValueError("Incomplete or extra notice export")
    return dict(status="source_notice_supplement_verified", files=len(seen), index=record(index),
                redistribution_clearance=False, publication_ready=False)


def export(inventory_path, output):
    inventory_path, output = Path(inventory_path), Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    source = record(inventory_path)
    inventory = json.loads(inventory_path.read_text())
    if inventory.get("status") != "distribution_source_packages_acquired_and_inventoried":
        raise ValueError("Require retained source material inventory")
    watched = [source, inventory["signature_receipt"], inventory["source_archive_inventory"]]
    archived = {}
    for package in inventory["source_packages"]:
        watched.append(package["download"]["file"])
        for member in package["members"]:
            ref = member["file"]
            archived[(ref["sha256"], ref["bytes"], ref["path"])] = ref
    files = []
    output.mkdir(parents=True, exist_ok=False)
    found = set()
    for group in inventory["selected_source_material"]:
        archive_ref = group["archive"]
        if (archive_ref["sha256"], archive_ref["bytes"], archive_ref["path"]) not in archived:
            raise ValueError("Archive is not bound to inventoried source package")
        expected = {}
        for row in group["selected_members"]:
            if row["member"] in SELECTION:
                if row["member"] in expected:
                    raise ValueError("Duplicate selected source record")
                expected[row["member"]] = row["file"]
        if not expected:
            continue
        check(archive_ref)
        watched.append(archive_ref)
        with tarfile.open(archive_ref["path"], mode="r|*") as archive:
            for member in archive:
                if member.name not in expected:
                    continue
                if member.name in found or not member.isfile() or member.size > 1_000_000:
                    raise ValueError("Duplicate or nonregular notice member")
                payload = archive.extractfile(member).read()
                identity = dict(bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())
                if any(identity[k] != expected[member.name][k] for k in identity):
                    raise ValueError("Archive notice does not match source inventory")
                found.add(member.name)
                relative = SELECTION[member.name]
                path = output / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                with path.open("xb") as handle:
                    handle.write(payload)
                files.append(dict(identity, member=member.name, relative_path=relative,
                                  archive=archive_ref))
        if not set(expected) <= found:
            raise ValueError("Missing selected source notice")
    if found != set(SELECTION):
        raise ValueError("Require all four selected component notices")
    for ref in watched:
        check(ref)
    save(output / "SOURCE_NOTICE_INDEX.json", dict(scope="selected_bundled_source_notices",
        files=sorted(files, key=lambda r: r["relative_path"]), inputs=watched,
        source=record(__file__), redistribution_clearance=False, publication_ready=False,
        limitations=["Exact selected source texts, not full transitive coverage or legal compatibility.",
            "Inventory and signature-receipt bytes are checked; RPM signature verification is not rerun here.",
            "Only four notice texts are exported, not headers, binaries, source archives or dataset payloads.",
            "Verify the index against the separately retained export receipt when relocating."]))
    return verify(output, record(output / "SOURCE_NOTICE_INDEX.json")["sha256"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    if args.receipt.exists() or args.receipt.is_symlink():
        raise FileExistsError(args.receipt)
    save(args.receipt, export(args.inventory.resolve(), args.output.absolute()))
