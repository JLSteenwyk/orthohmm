"""Check every locked version against a retained GitHub alert snapshot."""

import argparse
import hashlib
import json
from pathlib import Path
try:
    import tomllib
except ImportError:
    import tomli as tomllib

from packaging.specifiers import SpecifierSet
from packaging.utils import canonicalize_name


def evaluate(lock, alerts):
    versions = {}
    for package in lock["package"]:
        versions.setdefault(canonicalize_name(package["name"]), []).append(package["version"])
    rows = []
    for alert in alerts:
        name = canonicalize_name(alert["dependency"]["package"]["name"])
        affected = SpecifierSet(alert["vulnerable_version_range"])
        locked = versions.get(name, [])
        vulnerable = [version for version in locked if affected.contains(version, prereleases=True)]
        rows.append({"number": alert["number"], "ghsa_id": alert["ghsa_id"], "package": name,
                     "locked_versions": locked, "affected_locked_versions": vulnerable,
                     "status": "affected" if vulnerable else "not_in_lock" if not locked else "outside_reported_range"})
    return rows


def identity(path):
    data = path.read_bytes()
    return {"path": str(path.resolve()), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", type=Path, required=True)
    parser.add_argument("--alerts", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    alerts = json.loads(args.alerts.read_text())["alerts"]
    if not alerts or any(a["dependency"]["manifest_path"] != "docs/uv.lock" for a in alerts):
        raise ValueError("This audit requires a nonempty docs-only alert snapshot")
    with args.lock.open("rb") as handle:
        lock = tomllib.load(handle)
    rows = evaluate(lock, alerts)
    report = {"lock": identity(args.lock), "alert_snapshot": identity(args.alerts), "auditor": identity(Path(__file__)),
              "comparisons": rows, "affected_alerts": sum(r["status"] == "affected" for r in rows),
              "limitations": ["Tests only the retained advisory ranges, not all vulnerabilities or exploitability.",
                              "Checks every locked branch, even versions not selected on the current interpreter.",
                              "GitHub alert closure requires a separate post-push API check."]}
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({"alerts_checked": len(rows), "affected_alerts": report["affected_alerts"]}))
