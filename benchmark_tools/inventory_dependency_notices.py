"""Inventory declarations and notice bytes in report-pinned local wheels, not legal clearance."""

import argparse
from email import policy
from email.parser import BytesParser
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
from urllib.parse import unquote, urlsplit
import zipfile

from packaging.utils import canonicalize_name, parse_wheel_filename
from packaging.version import Version

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def notice_candidate(name):
    path = PurePosixPath(name)
    return "licenses" in [part.casefold() for part in path.parts] or bool(re.match(
        r"^(licen[cs]e|copying|notice|copyright|patents|authors)([._-]|$)", path.name, re.I))


def inventory(item):
    info = item["download_info"]
    url = urlsplit(info["url"])
    if url.scheme != "file" or url.netloc or url.query or url.fragment:
        raise ValueError("Require a local file URL without host/query/fragment")
    wheel = Path(unquote(url.path))
    if not wheel.is_absolute() or wheel.suffix != ".whl":
        raise ValueError("Require an absolute local wheel path")
    identity = record(wheel)
    hashes = info["archive_info"].get("hashes", {})
    digests = [hashes["sha256"]] if "sha256" in hashes else []
    legacy = info["archive_info"].get("hash")
    if legacy is not None:
        if not legacy.startswith("sha256="):
            raise ValueError("Require SHA-256 archive identity")
        digests.append(legacy.removeprefix("sha256="))
    if not digests or any(value != identity["sha256"] for value in digests):
        raise ValueError("Wheel differs from recorded archive identity")
    name, version, _, _ = parse_wheel_filename(wheel.name)
    expected = item["metadata"]
    if canonicalize_name(expected["name"]) != name or Version(expected["version"]) != version:
        raise ValueError("Wheel filename and install report differ")
    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise ValueError("Duplicate wheel members")
        # Vendored packages may retain nested metadata; identify the wheel's own distribution.
        metadata_names = [n for n in names if n.endswith(".dist-info/METADATA")
                          and len(PurePosixPath(n).parts) == 2]
        if len(metadata_names) != 1:
            raise ValueError("Require one wheel metadata member")
        metadata_path = metadata_names[0]
        metadata_raw = archive.read(metadata_path)
        metadata = BytesParser(policy=policy.default).parsebytes(metadata_raw)
        if (canonicalize_name(metadata["Name"]) != name
                or Version(metadata["Version"]) != version):
            raise ValueError("Embedded metadata differs from install report")
        prefix = metadata_path.rsplit("/", 1)[0]
        declared = []
        for value in metadata.get_all("License-File", []):
            value = str(value)
            relative = PurePosixPath(value)
            if relative.is_absolute() or ".." in relative.parts or "\\" in value:
                raise ValueError("Unexpected declared license path")
            matches = [p for p in (prefix + "/licenses/" + value, prefix + "/" + value)
                       if p in names and not p.endswith("/")]
            declared.append(dict(declaration=value, matching_members=matches))
        candidates = sorted({n for n in names if not n.endswith("/") and notice_candidate(n)} |
                            {n for entry in declared for n in entry["matching_members"]})
        notices = []
        for n in candidates:
            content = archive.read(n)
            notices.append(dict(path=n, bytes=len(content), sha256=hashlib.sha256(content).hexdigest()))
        license_text = str(metadata.get("License", ""))
        result = dict(name=str(metadata["Name"]), version=str(metadata["Version"]), wheel=identity,
            metadata_path=metadata_path, metadata_sha256=hashlib.sha256(metadata_raw).hexdigest(),
            license_expression=metadata.get("License-Expression"),
            license_field_short=license_text if len(license_text) <= 200 else None,
            license_field_characters=len(license_text),
            license_field_sha256=hashlib.sha256(license_text.encode()).hexdigest(),
            license_classifiers=[str(x) for x in metadata.get_all("Classifier", []) if str(x).startswith("License ::")],
            declared_license_files=declared, notice_candidates=notices,
            unresolved_declarations=[x for x in declared if len(x["matching_members"]) != 1],
            native_members=[n for n in names if re.search(r"\.(so(?:\.[0-9.]+)?|pyd|dll|dylib)$", n, re.I)],
            member_count=len(names), redistribution_clearance=False)
    check(identity)
    return result


def audit(install_report, expected_sha):
    evidence = record(install_report)
    if evidence["sha256"] != expected_sha:
        raise ValueError("Install report identity differs")
    report = json.loads(Path(install_report).read_text())
    rows = [inventory(item) for item in report["install"]]
    if not rows or len({canonicalize_name(row["name"]) for row in rows}) != len(rows):
        raise ValueError("Require nonempty, unique installed distribution inventory")
    for row in rows:
        check(row["wheel"])
    check(evidence)
    return dict(status="local_wheel_notice_inventory", install_report=evidence, wheels=rows,
        source=record(__file__), redistribution_clearance=False, publication_ready=False,
        limitations=["Reports provider metadata and notice-candidate bytes, not legal interpretation or compatibility.",
            "Filename heuristics can miss notices; embedded native/static dependencies need separate review.",
            "No archive extraction, installation, native execution or redistribution performed.",
            "Pinned wheel bytes verified against retained install report, not current environment or full runtime.",
            "Does not cover external executables, OS libraries, containers or benchmark datasets."])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--install-report", type=Path, required=True)
    parser.add_argument("--report-sha", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    result = audit(args.install_report, args.report_sha)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
