"""Summarize retained GNU-time stage records without claiming controlled timing."""

import argparse
import json
from pathlib import Path
import statistics

from benchmark_tools.gnu_time_companion import FIELDS, parse
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


LABELS = ("Elapsed (wall clock) time (h:mm:ss or m:ss)", "User time (seconds)",
          "System time (seconds)", "Maximum resident set size (kbytes)", "Exit status")


def verbose_time(text):
    values = {}
    for line in text.splitlines():
        line = line.strip()
        for label in LABELS:
            if line.startswith(label + ": "):
                if label in values:
                    raise ValueError("Duplicate GNU-time field")
                values[label] = line[len(label) + 2:]
    if set(values) != set(LABELS):
        raise ValueError("Incomplete GNU-time fields")
    parts = values[LABELS[0]].split(":")
    if len(parts) not in (2, 3):
        raise ValueError("Unexpected wall-time format")
    if not all(p.isdigit() for p in parts[:-1]):
        raise ValueError("Invalid elapsed hours/minutes")
    seconds = float(parts[-1])
    if not 0 <= seconds < 60 or (len(parts) == 3 and not 0 <= int(parts[1]) < 60):
        raise ValueError("Invalid elapsed components")
    elapsed = seconds + 60 * int(parts[-2]) + (3600 * int(parts[0]) if len(parts) == 3 else 0)
    converted = [str(elapsed), *[values[label] for label in LABELS[1:]]]
    result = parse("\n".join(f"{key}\t{value}" for key, value in zip(FIELDS, converted)))
    if result["exit_status"]:
        raise ValueError("Failed stage cannot enter completed resource panel")
    return result


def stage_resources(execution, names):
    expected = {s["name"] for s in execution["stages"]}
    if not names or not set(names) <= expected:
        raise ValueError("Missing required stage")
    result = []
    for name in names:
        candidates = [r for r in execution["logs"] if Path(r["path"]).name == name + ".time.txt"]
        stages = [s for s in execution["stages"] if s["name"] == name]
        if len(candidates) != 1 or len(stages) != 1 or stages[0]["returncode"]:
            raise ValueError("Invalid stage resource inventory")
        item = candidates[0]
        check(item)
        result.append(dict(name=name, source=item, **verbose_time(Path(item["path"]).read_text())))
    return dict(stages=result, elapsed_seconds=sum(r["elapsed_seconds"] for r in result),
                cpu_seconds=sum(r["user_seconds"] + r["system_seconds"] for r in result),
                max_reported_process_rss_kib=max(r["max_process_rss_kib"] for r in result))


def run(readback_path, output):
    if output.exists():
        raise FileExistsError(output)
    evidence = record(readback_path)
    if evidence["sha256"] != "af7a9e0cdba97b7b84c9fc6925e577594802b7bf74c99ac19840baf2e42594f4":
        raise ValueError("Unexpected admitted graph panel")
    readback = json.loads(readback_path.read_text())
    submission = json.loads(Path(readback["submission"]["path"]).read_text())
    check(submission["manifest"])
    manifest = json.loads(Path(submission["manifest"]["path"]).read_text())
    check(manifest["search_result"])
    search = json.loads(Path(manifest["search_result"]["path"]).read_text())
    records = []
    for row in readback["cells"]:
        check(row["execution"])
        execution = json.loads(Path(row["execution"]["path"]).read_text())
        cell = manifest["cells"][row["index"]]
        source = next(r for r in search["receipts"] if Path(r["path"]).parent.name == f'cell_{cell["dataset_index"]}')
        check(source)
        native = json.loads(Path(source["path"]).read_text())
        if native["status"] != "native_completed_pending_independent_readback":
            raise ValueError("Search did not complete")
        if row["arm"] == "hmm":
            phases = {"search": stage_resources(native, ["hmm"])}
        else:
            for suffix in ("_0", "_1"):
                names = [s["name"] for s in native["stages"] if s["name"].startswith("diamond_") and s["name"].endswith(suffix)]
                if len(names) != len(row["inputs"]):
                    raise ValueError("Incomplete species-specific DIAMOND stages")
            phases = {label: stage_resources(native, [s["name"] for s in native["stages"] if s["name"].startswith("diamond_") and s["name"].endswith(suffix)])
                      for label, suffix in (("database_preparation", "_0"), ("search", "_1"))}
        phases["graph"] = stage_resources(execution, ["graph"])
        records.append(dict(condition=row["condition"], seed=row["seed"], arm=row["arm"],
                            search_execution=source, graph_execution=row["execution"], phases=phases))
    summary = {}
    for arm in ("hmm", "diamond"):
        rows = [r for r in records if r["arm"] == arm]
        if len(rows) != 35:
            raise ValueError("Incomplete reporting resource panel")
        summary[arm] = {}
        for phase in rows[0]["phases"]:
            summary[arm][phase] = {}
            for metric in ("elapsed_seconds", "cpu_seconds", "max_reported_process_rss_kib"):
                values = [r["phases"][phase][metric] for r in rows]
                summary[arm][phase][metric] = dict(median=statistics.median(values), minimum=min(values), maximum=max(values))
    result = dict(status="descriptive_reporting_stage_resources", source=record(__file__), readback=evidence,
                  records=records, summary=summary, controlled_timing=False, publication_ready=False,
                  missing_measurements=["Offline numeric-export preparation", "Separate scoring resource measurement"],
                  limitations=["Concurrent shared-host runs; not matched controlled timing or a speedup claim.",
                               "DIAMOND searched broadly at E=1 then postfiltered; effort was not matched at selected E=1e-40.",
                               "RSS is the reported maximum process value, not simultaneous process-tree/cgroup memory.",
                               "Phase elapsed/CPU sums combine sequential native commands; wrapper/hash/copy overhead is excluded.",
                               "Small timings are rounded by GNU time; no confidence intervals or statistical ranking."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--readback", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.readback.resolve(), args.output.absolute())
