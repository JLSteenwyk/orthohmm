"""Replay immutable controller polls and verify first-terminal record retention."""

import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.capture_array_scheduler import terminal_records
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def replay(directory, array_id, count):
    if type(array_id) is not int or array_id <= 0 or type(count) is not int or count <= 0:
        raise ValueError("Require positive array identity and task count")
    source = record(__file__)
    capture_record = record(directory / "capture.json")
    capture = json.loads(Path(capture_record["path"]).read_text())
    if (capture["array_id"] != array_id or capture["expected_tasks"] != count
            or capture["status"] != "complete_controller_capture" or capture["missing_tasks"] != []
            or capture["scientific_timings_admitted"] is not False
            or type(capture["polls"]) is not int or capture["polls"] <= 0):
        raise ValueError("Require matching complete capture")
    expected = [directory / f"poll_{i:06d}.json" for i in range(capture["polls"])]
    if sorted(directory.glob("poll_*.json")) != expected:
        raise ValueError("Missing or extra raw polls")
    retained, evidence, errors = {}, [capture_record], 0
    tasks = set(range(count))
    for i, path in enumerate(expected):
        if path.is_symlink():
            raise ValueError("Symlinked raw poll")
        evidence.append(record(path))
        poll = json.loads(path.read_text())
        if (poll["argv"] != ["scontrol", "show", "job", str(array_id), "--oneliner"]
                or type(poll["observed_unix_ns"]) is not int or poll["observed_unix_ns"] <= 0):
            raise ValueError("Wrong poll command or timestamp")
        if "error_type" in poll or poll.get("returncode") != 0:
            errors += 1
            continue
        try:
            parsed = terminal_records(poll["stdout"], array_id, tasks)
        except ValueError as error:
            if poll.get("parse_error") != str(error):
                raise ValueError("Unrecorded or changed parse failure") from error
            errors += 1
            continue
        if "parse_error" in poll:
            raise ValueError("Spurious parse failure")
        for index, raw in parsed.items():
            if index in retained:
                continue
            path = directory / f"scheduler_{index}.txt"
            if path.is_symlink() or path.read_text() != raw:
                raise ValueError("Retained scheduler record differs from first terminal poll")
            evidence.append(record(path))
            retained[index] = {"path": path.name, "poll": i,
                               "sha256": hashlib.sha256(raw.encode()).hexdigest()}
        if set(retained) == tasks and i != len(expected) - 1:
            raise ValueError("Unexpected polls after completed capture")
    if ({str(k): v for k, v in retained.items()} != capture["retained"]
            or set(retained) != tasks or errors != capture["observation_errors"]):
        raise ValueError("Capture summary differs from replayed evidence")
    helper = record(Path(__file__).with_name("capture_array_scheduler.py"))
    for item in [source, helper, *evidence]:
        check(item)
    return {"status": "first_terminal_scheduler_capture_replayed", "array_id": array_id,
        "tasks": count, "polls": len(expected), "observation_errors": errors,
        "capture": capture_record, "retained": retained, "checked_records": evidence,
        "source": source, "helper": helper, "scientific_timings_admitted": False,
        "limitations": ["Uses the collector's tested terminal-record parser; replay is not an independent parser.",
                        "Does not establish allocation policy, output validity or uncontended timing."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--array-id", type=int, required=True)
    parser.add_argument("--tasks", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = replay(args.directory.resolve(), args.array_id, args.tasks)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
