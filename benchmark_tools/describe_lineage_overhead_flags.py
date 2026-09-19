"""Describe the complete audited lineage panel without attributing CPU residuals."""

import argparse
from collections import Counter
import gzip
import json
from pathlib import Path
import re

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.summarize_pressure_panel_flags import distribution

AUDIT_SHA = "7e4cf99291450c883bf69abe0d984d07737953c2792dcb191bdb5a3cbc73251e"
PERIODIC = {1, 3, 5, 6, 8, 10, 13, 15, 17}
FIELDS = ("residual_cores", "native_cores", "overhang_s", "root_minus_job_cpu_usec",
          "root_minus_system_cpu_usec", "system_minus_slurm_cpu_usec", "slurm_minus_job_cpu_usec")


def summarize(rows):
    return dict(count=len(rows), distributions={key: distribution([r[key] for r in rows])
                for key in FIELDS} if rows else {},
                reasons=dict(Counter(reason for r in rows for reason in r["reasons"])))


def describe(audit):
    runs = audit["runs"]
    if (audit["validated_tasks"] != 18 or [r["index"] for r in runs] != list(range(18))
            or any(r["status"] != "validated" for r in runs)):
        raise ValueError("Require complete validated panel; no survivor-only description")
    rows = []
    boundary = []
    for run in runs:
        if run["index"] not in PERIODIC:
            if (run["mode"] != "boundary" or run["interval_screening_available"] is not False
                    or run["narrow_flagged_intervals"] is not None):
                raise ValueError("Boundary interval coverage changed")
            boundary.append(run["index"])
            continue
        if run["mode"] != "periodic" or run["interval_screening_available"] is not True:
            raise ValueError("Missing periodic interval coverage")
        screen = run["lineage_screening"]
        intervals = screen["intervals"]
        flags = [i for i, x in enumerate(intervals) if not x["narrow"]["screen_passed"]]
        if (flags != run["narrow_flagged_intervals"] or flags != screen["narrow_flagged_intervals"]
                or len(intervals) != run["observation_points"] - 1):
            raise ValueError("Incomplete intervals or contradictory flags")
        for i, interval in enumerate(intervals):
            narrow, lineage = interval["narrow"], interval["lineage"]
            complements = lineage["signed_complements"]
            if (len(complements) != 3 or complements[0]["ancestor"] != "/"
                    or complements[0]["excluded_child"] != "/system.slice"
                    or complements[1]["ancestor"] != complements[0]["excluded_child"]
                    or complements[2]["ancestor"] != complements[1]["excluded_child"]
                    or not re.fullmatch(r"/system\.slice/spark-7ff0_slurmstepd\.scope/job_[0-9]+",
                                        complements[2]["excluded_child"])):
                raise ValueError("Unexpected lineage path")
            values = [x["signed_cpu_usec"] for x in complements]
            if sum(values) != lineage["root_minus_target_cpu_usec"]:
                raise ValueError("Adjacent signed complements do not telescope")
            rows.append(dict(task=run["index"], method=run["method"], interval=i,
                flagged=i in flags, reasons=narrow["reasons"],
                residual_cores=narrow["signed_unassigned_average_cores"],
                native_cores=narrow["native_cpu_s"] / narrow["wall_s"],
                overhang_s=narrow["outer_read_overhang_s"],
                root_minus_job_cpu_usec=lineage["root_minus_target_cpu_usec"],
                **dict(zip(FIELDS[4:], values))))
    groups = []
    for method in sorted({r["method"] for r in rows}):
        selected = [r for r in rows if r["method"] == method]
        groups.append(dict(method=method, all=summarize(selected),
            flagged=summarize([r for r in selected if r["flagged"]]),
            unflagged=summarize([r for r in selected if not r["flagged"]])))
    return dict(status="lineage_overhead_flags_described_not_causally_attributed",
        periodic_tasks=sorted(PERIODIC), boundary_tasks_without_interval_coverage=boundary,
        all=summarize(rows), flagged=summarize([r for r in rows if r["flagged"]]),
        unflagged=summarize([r for r in rows if not r["flagged"]]), methods=groups,
        flagged_intervals=[r for r in rows if r["flagged"]],
        scientific_timings_admitted=False, environmental_validity_established=False,
        limitations=["Post-outcome description of all audited periodic intervals, not causal inference.",
            "Signed non-atomic complements are not bounds or identified outside-process CPU times.",
            "Host/native residuals and ancestor complements use different read windows; do not subtract them.",
            "All original flags remain in the source audit; no thresholds or timing eligibility changed.",
            "Boundary arms have no interval coverage; no flags are inferred for them."])


def report(path):
    evidence = record(path)
    if evidence["sha256"] != AUDIT_SHA:
        raise ValueError("Wrong retained complete-panel audit")
    result = describe(json.loads(gzip.decompress(path.read_bytes())))
    check(evidence)
    return dict(result, audit=evidence, source=record(__file__))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = report(args.audit.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
