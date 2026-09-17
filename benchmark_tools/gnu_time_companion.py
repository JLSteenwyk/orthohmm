"""Explicit GNU-time wrapper and strict companion accounting, not total RSS."""

import csv
import io
import math
from pathlib import Path

FIELDS = ("elapsed_seconds", "user_seconds", "system_seconds", "max_process_rss_kib", "exit_status")
FORMAT = "\n".join(f"{key}\t{specifier}" for key, specifier in zip(FIELDS, ("%e", "%U", "%S", "%M", "%x")))


def command(native, output, executable="/usr/bin/time"):
    if (not native or any(not isinstance(arg, str) or not arg or "\0" in arg for arg in native)
            or not Path(output).is_absolute() or not Path(executable).is_absolute()):
        raise ValueError("Require explicit command and absolute GNU-time paths")
    return [str(executable), "-q", "-f", FORMAT, "-o", str(output), "--", *native]


def parse(text):
    rows = list(csv.reader(io.StringIO(text), delimiter="\t"))
    if len(rows) != len(FIELDS) or any(len(row) != 2 for row in rows) or [r[0] for r in rows] != list(FIELDS):
        raise ValueError("Unexpected GNU-time companion fields")
    result = {}
    for i, (key, text_value) in enumerate(rows):
        try:
            value = float(text_value) if i < 3 else int(text_value)
        except ValueError as error:
            raise ValueError("Invalid GNU-time numeric value") from error
        if not math.isfinite(value) or value < 0:
            raise ValueError("Invalid GNU-time numeric value")
        result[key] = value
    if result["exit_status"] > 255:
        raise ValueError("Invalid GNU-time exit status")
    return {**result, "semantics": {
        "elapsed_seconds": "GNU-time native child wall time, rounded by GNU time; wrapper launch overhead is excluded",
        "user_seconds": "Native child and waited-for descendants, as reported by GNU time; collector CPU excluded",
        "system_seconds": "Native child and waited-for descendants, as reported by GNU time; collector CPU excluded",
        "max_process_rss_kib": "Reported maximum process RSS, not simultaneous sum across native processes or cgroup memory"}}
