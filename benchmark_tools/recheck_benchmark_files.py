"""Stream-check selected retained files without upgrading historical provenance."""

import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import stat
import time


def require(condition, message):
    if not condition:
        raise ValueError(message)


def identity(info):
    return [info.st_dev, info.st_ino, info.st_mode, info.st_size,
            info.st_mtime_ns, info.st_ctime_ns]


def inspect(ref):
    path = Path(ref["path"])
    try:
        before = path.lstat()
        if not stat.S_ISREG(before.st_mode):
            return {"status": "nonregular", "before": identity(before)}
        fd = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
        with os.fdopen(fd, "rb") as handle:
            opened = os.fstat(handle.fileno())
            digest = hashlib.sha256()
            count = 0
            while chunk := handle.read(4 * 1024 * 1024):
                digest.update(chunk)
                count += len(chunk)
            finished = os.fstat(handle.fileno())
        after = path.lstat()
        stable = identity(before) == identity(opened) == identity(finished) == identity(after)
        observed = {"path": ref["path"], "bytes": count, "sha256": digest.hexdigest()}
        status = "changed_during_check" if not stable else "matches" if observed == ref else "mismatch"
        return {"status": status, "observed": observed, "before": identity(before),
                "after": identity(after), "descriptor_stable": stable}
    except FileNotFoundError:
        return {"status": "missing"}
    except OSError as error:
        return {"status": "unreadable", "error": str(error), "errno": error.errno}


def selections(register):
    rows = register["rows"]
    datasets = {"OrthoBench", "QfO", "ThreeKingdoms"}
    require(len(rows) == 24 and {row["dataset"] for row in rows} == datasets,
            "Require the complete three-dataset panel")
    methods = []
    for dataset in sorted(datasets):
        names = [row["key"] for row in rows if row["dataset"] == dataset]
        require(len(names) == 8 and len(set(names)) == 8, "Require eight distinct methods")
        methods.append(set(names))
    require(methods[0] == methods[1] == methods[2], "Method panels differ")
    pins, associations = {}, []
    for row in rows:
        selected = {"dataset": row["dataset"], "key": row["key"], "files": []}
        for role in ("input_records", "output_records"):
            data = row[role]
            require(isinstance(data, (list, dict)), "Unknown record container")
            records = data.items() if isinstance(data, dict) else enumerate(data)
            for name, raw in records:
                require(isinstance(raw, dict) and {"path", "bytes", "sha256"} <= raw.keys(),
                        "Incomplete file identity")
                ref = {k: raw[k] for k in ("path", "bytes", "sha256")}
                require(isinstance(ref["path"], str) and Path(ref["path"]).is_absolute()
                        and type(ref["bytes"]) is int and ref["bytes"] >= 0
                        and isinstance(ref["sha256"], str) and len(ref["sha256"]) == 64
                        and all(c in "0123456789abcdef" for c in ref["sha256"]), "Invalid file identity")
                require(ref["path"] not in pins or pins[ref["path"]] == ref,
                        "Conflicting selected file identity")
                pins[ref["path"]] = ref
                selected["files"].append({"role": role, "entry": name, "path": ref["path"]})
        associations.append(selected)
    return pins, associations


def run(register_path, expected_sha):
    path = Path(register_path).absolute()
    require(path.is_file() and not path.is_symlink(), "Nonregular register")
    content = path.read_bytes()
    require(hashlib.sha256(content).hexdigest() == expected_sha, "Register checksum differs")
    ref = {"path": str(path), "bytes": len(content), "sha256": expected_sha}
    pins, associations = selections(json.loads(content))
    started = time.time_ns()
    observations = {name: {"expected": pin, **inspect(pin)} for name, pin in sorted(pins.items())}
    require(path.read_bytes() == content, "Register changed during inspection")
    counts = dict(Counter(row["status"] for row in observations.values()))
    return {"schema": "selected_benchmark_current_file_recheck_v1", "register": ref,
            "started_unix_ns": started, "finished_unix_ns": time.time_ns(),
            "observer_affinity": sorted(os.sched_getaffinity(0)),
            "selected_rows": associations, "observations": observations,
            "status_counts": counts, "selected_unique_bytes": sum(r["bytes"] for r in pins.values()),
            "all_selected_files_match": counts == {"matches": len(pins)},
            "publication_ready": False, "scores_recomputed": False,
            "historical_input_consumption_proven": False, "complete_transitive_provenance": False,
            "limitations": [
                "Current direct input/output bytes only; shared files are checked once, not independent repeats.",
                "Historical consumption, executables, environments, conversion and scoring are not re-admitted.",
                "Missing, mismatched, nonregular, unreadable and unstable files remain distinct outcomes.",
                "Stat/descriptor stability brackets each read; this is not continuous or cross-file atomic integrity.",
                "Register input roles include downstream artifacts such as OrthoMCL BPO indexes, not only FASTAs.",
                "This is a bookkeeping workload on the shared host, not comparative performance measurement."]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--register", required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    output = Path(args.output)
    require(not output.exists() and not output.is_symlink(), "Output already exists")
    result = run(args.register, args.sha256)
    require(output.resolve() not in {Path(p).resolve() for p in result["observations"]},
            "Output aliases selected evidence")
    source = Path(__file__).resolve()
    data = source.read_bytes()
    result["source"] = {"path": str(source), "bytes": len(data),
                        "sha256": hashlib.sha256(data).hexdigest()}
    with output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({"files": len(result["observations"]), "status_counts": result["status_counts"],
                      "all_selected_files_match": result["all_selected_files_match"]}))


if __name__ == "__main__":
    main()
