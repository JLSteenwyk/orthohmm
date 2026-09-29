"""Inspect public QfO fork default-branch inventories for original TreeFam leads."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen


def candidates(tree):
    if tree.get("truncated") is not False or not isinstance(tree.get("tree"), list):
        raise ValueError("Incomplete recursive inventory")
    return [r for r in tree["tree"] if r.get("type") == "blob" and
            ("treefam" in r["path"].lower() or r["path"].lower().endswith(".nhx"))]


def inspect(output):
    output.mkdir(parents=True, exist_ok=False)
    requests, rows = [], []

    def fetch(url, name):
        entry = dict(url=url, started_utc=datetime.now(timezone.utc).isoformat())
        requests.append(entry)
        try:
            with urlopen(Request(url, headers={"User-Agent": "OrthoHMM-public-source-audit"}), timeout=30) as response:
                raw = response.read(16 * 1024 * 1024 + 1)
                if len(raw) > 16 * 1024 * 1024:
                    raise ValueError("Response exceeds bound")
                entry.update(status=response.status, final_url=response.url)
            path = output / name
            path.write_bytes(raw)
            entry.update(file=str(path.resolve()), bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())
            return json.loads(raw)
        except Exception as error:
            entry.update(error_type=type(error).__name__, error=str(error))
            return None

    forks = fetch("https://api.github.com/repos/qfo/benchmark-webservice/forks?per_page=100", "forks.json")
    if isinstance(forks, list):
        for i, fork in enumerate(forks):
            repo, branch = fork["full_name"], fork["default_branch"]
            url = f"https://api.github.com/repos/{repo}/git/trees/{quote(branch, safe='')}?recursive=1"
            tree = fetch(url, f"tree_{i:02d}.json")
            row = dict(repository=repo, default_branch=branch, status="unresolved")
            if tree is not None:
                try:
                    row.update(status="inventory_inspected", tree_sha=tree["sha"],
                               entries=len(tree["tree"]), candidates=candidates(tree))
                except (ValueError, KeyError, TypeError) as error:
                    row["error"] = str(error)
            rows.append(row)
    report = dict(status="public_fork_default_branch_search", requests=requests, forks=rows,
        pagination_complete=isinstance(forks, list) and len(forks) < 100,
        original_inputs_admitted=False, publication_ready=False,
        limitations=["Default-branch inventories only; not full fork history, other branches or deleted/private objects.",
            "Filename leads are not validation of contents or equivalence to the retained QfO reference.",
            "No contact, authenticated access, code execution or benchmark changes."])
    (output / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = inspect(args.output)
    print(json.dumps(dict(forks=len(result["forks"]), rows=result["forks"], pagination_complete=result["pagination_complete"])))
