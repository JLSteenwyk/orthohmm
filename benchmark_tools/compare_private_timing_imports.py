"""Compare a private candidate's imported files with retained native lookup bytes."""

import argparse
import hashlib
import json
from pathlib import Path


def record(path):
    return dict(path=str(path.resolve()), bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def compare(old, new, core_root):
    pinned = {r["path"]: r for r in old["files"]}
    selected = {name: path for name, path in old["modules"].items()
                if "site-packages" in Path(path).parts or Path(path).is_relative_to(core_root)}
    rows, absent = [], []
    for name, old_path in sorted(selected.items()):
        if name not in new["modules"]:
            absent.append(dict(module=name, prior_path=old_path))
            continue
        if old_path not in pinned:
            raise ValueError("Prior imported file lacks a hash")
        observed = record(Path(new["modules"][name]))
        rows.append(dict(module=name, prior=pinned[old_path], current=observed,
                         identical_sha256=pinned[old_path]["sha256"] == observed["sha256"]))
    return dict(compared=len(rows), identical=sum(r["identical_sha256"] for r in rows),
        changed=[r for r in rows if not r["identical_sha256"]], omitted=absent, comparisons=rows,
        additional_module_names=sorted(set(new["modules"])-set(old["modules"])))


def run(prior, expected_sha, installation, output):
    if output.exists():
        raise FileExistsError(output)
    if record(prior)["sha256"] != expected_sha:
        raise ValueError("Historical import report checksum differs")
    installed = json.loads(installation.read_text())
    source = Path(installed["import_report"]["path"])
    if record(source) != installed["import_report"]:
        raise ValueError("Candidate import report changed")
    baseline = Path(installed["baseline"]["path"])
    if record(baseline) != installed["baseline"]:
        raise ValueError("Scientific baseline changed")
    configuration = json.loads(baseline.read_text())
    checked_core = []
    for ref in configuration["core_sources"]:
        current = record(Path(ref["absolute_path"]))
        if current["sha256"] != ref["sha256"] or current["bytes"] != ref["bytes"]:
            raise ValueError("Frozen scientific source changed")
        checked_core.append(current)
    result = compare(json.loads(prior.read_text()), json.loads(source.read_text()), Path(configuration["core_root"]))
    result.update(status="private_timing_import_comparison", core_sources_checked=checked_core,
        sources=[record(prior), record(installation), record(source), record(baseline), record(Path(__file__))],
        scientific_execution_authorized=False,
        limitations=["Observed Python/module files only, not all package payloads, ELF dependency closure or every inference branch.",
                     "Omitted startup hooks and optional modules remain explicit; same versions do not establish same binaries.",
                     "No inference output parity, native collector integration or timing admission."])
    with output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prior", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--installation", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.prior.resolve(), args.sha256, args.installation.resolve(), args.output.resolve())
    print(json.dumps(dict(compared=result["compared"], identical=result["identical"],
                         changed=[r["module"] for r in result["changed"]], omitted=[r["module"] for r in result["omitted"]])))
