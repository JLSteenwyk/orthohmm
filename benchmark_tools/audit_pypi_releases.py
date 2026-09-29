"""Snapshot release-specific PyPI advisories without installing packages."""

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import urlopen

from packaging.requirements import Requirement
from packaging.utils import canonicalize_name
from packaging.version import Version

from benchmark_tools.audit_dependency_lock import identity


def read_pins(text):
    pins = []
    seen = set()
    for line in text.splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split()
        req = Requirement(parts[0])
        specs = list(req.specifier)
        if (req.url or req.marker or req.extras or len(specs) != 1
                or specs[0].operator != "==" or "*" in specs[0].version):
            raise ValueError("Expected an unconditional exact pin")
        name = canonicalize_name(req.name)
        if name in seen:
            raise ValueError("Duplicate package")
        seen.add(name)
        hashes = []
        for part in parts[1:]:
            prefix = "--hash=sha256:"
            value = part.removeprefix(prefix)
            if not part.startswith(prefix) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
                raise ValueError("Expected SHA256 wheel hash")
            hashes.append(value)
        if not hashes:
            raise ValueError("Missing artifact hash")
        pins.append(dict(name=name, version=str(Version(specs[0].version)), hashes=hashes))
    if not pins:
        raise ValueError("Empty lock")
    return pins


def review(pin, payload):
    info = payload["info"]
    if (canonicalize_name(info["name"]) != pin["name"]
            or Version(info["version"]) != Version(pin["version"])):
        raise ValueError("Response release identity mismatch")
    advisories = payload["vulnerabilities"]
    if not isinstance(advisories, list) or any(not isinstance(a, dict) or not a.get("id") for a in advisories):
        raise ValueError("Invalid advisory inventory")
    artifacts = payload["urls"]
    if not isinstance(artifacts, list):
        raise ValueError("Invalid artifact inventory")
    matches = [dict(filename=a["filename"], sha256=a["digests"]["sha256"],
                    yanked=a["yanked"]) for a in artifacts
               if a["digests"]["sha256"] in pin["hashes"]]
    return dict(status="queried", advisories=advisories,
                active_advisories=[a["id"] for a in advisories if not a.get("withdrawn")],
                matching_public_artifacts=matches,
                unmatched_locked_hashes=sorted(set(pin["hashes"]) - {a["sha256"] for a in matches}))


def audit(locks, output):
    output.mkdir(parents=True, exist_ok=False)
    records = [identity(p) for p in locks]
    inventories = [read_pins(p.read_text()) for p in locks]
    cache = {}
    rows = []
    for lock, pins in zip(records, inventories):
        for pin in pins:
            key = (pin["name"], pin["version"])
            url = "https://pypi.org/pypi/{}/{}/json".format(*(quote(v, safe="") for v in key))
            if key not in cache:
                timestamp = datetime.now(timezone.utc).isoformat()
                try:
                    with urlopen(url, timeout=45) as response:
                        raw = response.read()
                    path = output / ("{}-{}.json".format(*key))
                    path.write_bytes(raw)
                    cache[key] = dict(payload=json.loads(raw), response=identity(path), fetched_at=timestamp)
                except (HTTPError, URLError, TimeoutError) as exc:
                    cache[key] = dict(error=str(exc), fetched_at=timestamp)
            saved = cache[key]
            row = dict(lock=lock["path"], **pin, url=url, fetched_at=saved["fetched_at"])
            if "error" in saved:
                row.update(status="unresolved", error=saved["error"])
            else:
                row.update(review(pin, saved["payload"]), response=saved["response"])
            rows.append(row)
    assert records == [identity(p) for p in locks], "Locks changed during audit"
    report = dict(locks=records, auditor=identity(Path(__file__)), rows=rows,
                  unique_queries=len(cache), comprehensive_security_clearance=False,
                  limitations=["A missing release or unmatched artifact is not a clean result.",
                               "Advisories describe public package versions, not equivalence of custom builds.",
                               "No dependency resolution, live payload verification, reachability, native-library or OS scan.",
                               "No known advisory is not evidence of absence of vulnerabilities."])
    (output / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", type=Path, action="append", required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.lock, args.output_directory)
    print(json.dumps(dict(queries=result["unique_queries"], rows=len(result["rows"]),
                          unresolved=sum(r["status"] != "queried" for r in result["rows"]))))
