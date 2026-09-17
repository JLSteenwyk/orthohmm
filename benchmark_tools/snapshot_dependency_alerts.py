"""Read-only GitHub dependency-alert snapshot without credentials or exploit text."""

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import urllib.request


def summarize(alert):
    advisory = alert["security_advisory"]
    vulnerability = alert["security_vulnerability"]
    return {"number": alert["number"], "state": alert["state"], "dependency": alert["dependency"],
            "ghsa_id": advisory["ghsa_id"], "cve_id": advisory["cve_id"],
            "severity": advisory["severity"], "summary": advisory["summary"],
            "vulnerable_version_range": vulnerability["vulnerable_version_range"],
            "first_patched_version": vulnerability["first_patched_version"],
            "url": alert["html_url"], "references": [r["url"] for r in advisory["references"]]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--git-credential", action="store_true")
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    token = os.environ.get("GH_TOKEN") or os.environ.get("GITHUB_TOKEN")
    if args.git_credential:
        result = subprocess.run(["git", "credential", "fill"], input="protocol=https\nhost=github.com\n\n",
                                text=True, capture_output=True, check=True)
        token = dict(line.split("=", 1) for line in result.stdout.splitlines() if "=" in line).get("password")
    headers = {"Accept": "application/vnd.github+json", "X-GitHub-Api-Version": "2022-11-28"}
    if token:
        headers["Authorization"] = "Bearer " + token
    base = "https://api.github.com/repos/JLSteenwyk/orthohmm/dependabot/alerts"
    alerts, page = [], 1
    while True:
        suffix = "" if page == 1 else f"&page={page}"
        request = urllib.request.Request(f"{base}?state=open&per_page=100{suffix}", headers=headers)
        with urllib.request.urlopen(request, timeout=30) as response:
            batch = json.load(response)
        alerts.extend(summarize(a) for a in batch)
        if len(batch) < 100:
            break
        page += 1
    report = {"repository": "JLSteenwyk/orthohmm", "endpoint": base, "state_filter": "open",
              "retrieved_utc": datetime.now(timezone.utc).isoformat(), "alerts": alerts,
              "limitations": "Repository alerts are not a full installed-environment or reachability audit."}
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({"open_alerts": len(alerts), "manifests": sorted({a["dependency"]["manifest_path"] for a in alerts}),
                      "packages": sorted({a["dependency"]["package"]["name"].lower() for a in alerts})}))


if __name__ == "__main__":
    main()
