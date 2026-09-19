"""Describe archived frontier changes without repairing or admitting measurements."""

import argparse
import json
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_cgroup_frontier import validate
from benchmark_tools.summarize_frontier_overhead import terminal_scheduler_rows


def changes(left, right):
    for point in (left, right):
        validate(point)
    a = left["inventory_before"]["identities"]
    b = right["inventory_before"]["identities"]
    return dict(
        boot_changed=left["boot_id"] != right["boot_id"],
        target_changed=left["target"] != right["target"],
        added={k: b[k] for k in sorted(b.keys() - a.keys())},
        removed={k: a[k] for k in sorted(a.keys() - b.keys())},
        replaced={k: dict(before=a[k], after=b[k])
                  for k in sorted(a.keys() & b.keys()) if a[k] != b[k]},
    )


def diagnose(archive):
    accounting = archive / "accounting.txt"
    scheduler = terminal_scheduler_rows(accounting.read_text(), 21889)
    evidence = [record(accounting)]
    runs = []
    for index in range(18):
        directory = archive / "pressure_frontier_overhead_v2" / f"run_{index:02d}"
        paths = sorted((directory / "measurement").glob("point_*.json"))
        points, invalid, transitions = [], [], []
        for path in paths:
            evidence.append(record(path))
            point = json.loads(path.read_text())["frontier"]
            try:
                validate(point)
            except (ValueError, KeyError, TypeError) as error:
                invalid.append(dict(point=path.name, error=str(error)))
                points.append(None)
            else:
                points.append(point)
        for i, (left, right) in enumerate(zip(points, points[1:])):
            if left is None or right is None:
                continue
            delta = changes(left, right)
            if any(delta.values()):
                transitions.append(dict(left=paths[i].name, right=paths[i+1].name, **delta))
        runs.append(dict(index=index, scheduler=scheduler[index], points=len(points),
                         invalid_points=invalid, transitions=transitions))
    for item in evidence:
        check(item)
    return dict(status="frontier_identity_diagnostic_only", runs=runs, evidence=evidence,
                source=record(__file__), scientific_timings_admitted=False,
                limitations=["Observed identities only, not service activity or causal timing interference.",
                             "Stable snapshots cannot exclude between-sample transient scopes.",
                             "Failures and missing observations are not repaired or replaced."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = diagnose(args.archive.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
