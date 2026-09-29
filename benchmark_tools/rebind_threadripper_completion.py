"""Reinventory the existing runtime for the native-completion collector change."""

import argparse
import json
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.snapshot_runtime_trees import inventory


def changes(old, new, allowed):
    before = {r["path"]: r for r in old}
    after = {r["path"]: r for r in new}
    if len(before) != len(old) or len(after) != len(new) or before.keys() != after.keys():
        raise ValueError("Runtime entry inventory changed")
    changed = [dict(path=p, before=before[p], after=after[p]) for p in sorted(before) if before[p] != after[p]]
    if any(row["path"] not in allowed for row in changed):
        raise ValueError("Unexpected runtime change: " + str([r["path"] for r in changed]))
    return changed


def write(path, value):
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def run(previous, expected_sha, work, output):
    if output.exists() or work.exists():
        raise FileExistsError("Require fresh work and output paths")
    previous = previous.resolve()
    prior = read_frozen(previous, expected_sha)
    source = Path(__file__).resolve().parent
    check(prior["baseline"])
    check(prior["command_plan"])
    check(prior["generator_source"])
    allowed = {str(source / name) for name in ("measure_threadripper_scaling.py", "replay_threadripper_scaling.py")}
    work.mkdir(parents=True)
    manifests, changed, old_records, new_records = [], [], [], []
    for i, (path, sha) in enumerate(prior["runtime_specs"]):
        old = read_frozen(Path(path), sha)
        new = inventory(old["roots"])
        target = work / f"{i}.json"
        write(target, new)
        try:
            changed.extend(changes(old["records"], new["records"], allowed))
        except ValueError as error:
            old_paths = {r["path"] for r in old["records"]}
            new_paths = {r["path"] for r in new["records"]}
            write(work / "rejection.json", dict(status="runtime_rebinding_rejected",
                reason=str(error), prior=record(Path(path)), observed=record(target),
                added=sorted(new_paths-old_paths), removed=sorted(old_paths-new_paths),
                scientific_execution_authorized=False))
            raise
        old_records.extend(old["records"])
        new_records.extend(new["records"])
        manifests.append(record(target))
    if {r["path"] for r in changed} != allowed:
        raise ValueError("Expected exactly the collector and replay change")
    helper = inventory([source / "native_completion.py"])
    if set(r["path"] for r in helper["records"]) & set(r["path"] for r in new_records):
        raise ValueError("Completion helper already covered")
    new_records.extend(helper["records"])
    target = work / "completion_helper.json"
    write(target, helper)
    manifests.append(record(target))
    for row in prior["baseline_paths"]:
        old = [r for r in old_records if r["path"] == row]
        new = [r for r in new_records if r["path"] == row]
        if len(old) != 1 or old != new:
            raise ValueError("Scientific baseline path changed or missing")
    result = dict(status="prospective_threadripper_completion_runtime_bound", supersedes=record(previous),
        baseline=prior["baseline"], command_plan=prior["command_plan"], baseline_paths=prior["baseline_paths"],
        generator_source=record(source / "snapshot_runtime_trees.py"), source=record(Path(__file__).resolve()),
        changed_records=changed, added_helper_records=helper["records"], records=len(new_records),
        runtime_manifests=manifests, runtime_specs=[[r["path"], r["sha256"]] for r in manifests],
        scientific_execution_authorized=False, publication_ready=False,
        limitations=["Fresh explicit inventories, not a hermetic OS snapshot or continuous mutation guard.",
                     "Only collector/replayer changed; scientific baseline remains unchanged. Native lookup and fixture validation still required."])
    output.parent.mkdir(parents=True, exist_ok=True)
    write(output, result)
    return record(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--previous", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.previous, args.sha256, args.work.resolve(), args.output.resolve())))
