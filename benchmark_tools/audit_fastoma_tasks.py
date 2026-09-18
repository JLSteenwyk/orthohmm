"""Strict fresh-run FastOMA task-trace checks for native-output admission."""

from collections import Counter
import csv
import io
import math
from pathlib import Path
import re
import shlex

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_qfo_corrected_fastoma_assets import IMAGE_DIGEST

FIXED = {"check_input", "infer_roothogs", "batch_roothogs", "collect_subhogs",
         "extract_pairwise_ortholog_relations", "fastoma_report"}
ALLOWED = FIXED | {"omamer_run", "hog_big", "hog_rest"}
PROGRAMS = {"check_input": "fastoma-check-input", "infer_roothogs": "fastoma-infer-roothogs",
            "batch_roothogs": "fastoma-batch-roothogs", "collect_subhogs": "fastoma-collect-subhogs",
            "extract_pairwise_ortholog_relations": "fastoma-helper", "fastoma_report": "papermill",
            "omamer_run": "omamer", "hog_big": "fastoma-infer-subhogs", "hog_rest": "fastoma-infer-subhogs"}


def option(argv, flag):
    if argv.count(flag) != 1 or argv.index(flag) + 1 == len(argv):
        raise ValueError("Missing or repeated option: " + flag)
    return argv[argv.index(flag) + 1]


def wrapper_limits(text):
    commands = [shlex.split(line.strip()) for line in text.splitlines() if line.strip().startswith("docker run ")]
    if len(commands) != 1 or commands[0].count(IMAGE_DIGEST) != 1:
        raise ValueError("Require one pinned Docker invocation")
    argv = commands[0]
    if "--privileged" in argv or option(argv, "--network") != "none":
        raise ValueError("Unexpected container privilege/network")
    cpus = float(option(argv, "--cpus"))
    match = re.fullmatch(r"([0-9]+)([kmg]?)", option(argv, "--memory"), re.IGNORECASE)
    if match is None or not math.isfinite(cpus) or not 0 < cpus <= 180:
        raise ValueError("Invalid container resource limit")
    memory = int(match[1]) * {"": 1, "k": 1024, "m": 1024**2, "g": 1024**3}[match[2].lower()]
    if not 0 < memory <= 700 * 1024**3:
        raise ValueError("Invalid container memory limit")
    return {"cpus": cpus, "memory_bytes": memory}


def task_directory(work, short_hash):
    if not re.fullmatch(r"[0-9a-f]{2}/[0-9a-f]{6}", short_hash):
        raise ValueError("Invalid Nextflow trace work hash")
    matches = list(work.glob(short_hash + "*"))
    if len(matches) != 1 or not matches[0].is_dir() or matches[0].is_symlink():
        raise ValueError("Missing or ambiguous native task directory")
    return matches[0]


def inspect_task(directory):
    paths = [directory / name for name in (".command.sh", ".command.run", ".exitcode", ".command.log")]
    files = [record(path) for path in paths]
    if (directory / ".exitcode").read_text().strip() != "0":
        raise ValueError("Nonzero task exit code")
    limits = wrapper_limits((directory / ".command.run").read_text())
    script = (directory / ".command.sh").read_text()
    for item in files:
        check(item)
    return {"directory": str(directory), "files": files, "limits": limits}, shlex.split(script, comments=True)


def validate_trace(text, species_files):
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if not rows or not {"task_id", "hash", "name", "status", "exit"}.issubset(rows[0]):
        raise ValueError("Missing Nextflow trace columns/tasks")
    if len(set(species_files)) != len(species_files) or not species_files:
        raise ValueError("Require unique proteome filenames")
    task_ids, hashes, counts = set(), set(), Counter()
    for row in rows:
        name = row["name"].split(" (", 1)[0]
        if name not in ALLOWED or row["status"] != "COMPLETED" or row["exit"] != "0":
            raise ValueError("Unrecognized, failed, cached or retried task requires explicit review")
        if not row["task_id"].isdigit() or row["task_id"] in task_ids or row["hash"] in hashes:
            raise ValueError("Duplicate or invalid native task identity")
        task_ids.add(row["task_id"])
        hashes.add(row["hash"])
        counts[name] += 1
    if (any(counts[name] != 1 for name in FIXED) or counts["omamer_run"] != len(species_files)
            or counts["hog_big"] + counts["hog_rest"] == 0):
        raise ValueError("Incomplete or duplicated FastOMA process coverage")
    return rows, counts


def audit_tasks(trace, work, species_files):
    trace_record = record(trace)
    rows, counts = validate_trace(trace.read_text(), species_files)
    tasks, queries, batches = [], set(), {"hog_big": set(), "hog_rest": set()}
    batch_directory = None
    for row in rows:
        name = row["name"].split(" (", 1)[0]
        directory = task_directory(work, row["hash"])
        info, argv = inspect_task(directory)
        if PROGRAMS[name] not in argv:
            raise ValueError("Task script does not contain its expected native program")
        if name == "omamer_run":
            query = option(argv, "--query")
            if query not in species_files or query in queries:
                raise ValueError("Duplicate or foreign OMAmer query")
            queries.add(query)
        if name in batches:
            batch = option(argv, "--input-rhog-folder")
            if Path(batch).name != batch or batch in batches[name]:
                raise ValueError("Duplicate or unexpected HOG batch")
            batches[name].add(batch)
        if name == "batch_roothogs":
            batch_directory = directory
        if name == "extract_pairwise_ortholog_relations" and option(argv, "--type") != "ortholog":
            raise ValueError("Native pair output is not orthology")
        tasks.append({"trace": row, **info})
    if queries != set(species_files):
        raise ValueError("OMAmer does not cover exact proteome inventory")
    for name, folder in (("hog_big", "rhogs_big"), ("hog_rest", "rhogs_rest")):
        expected = {p.name for p in (batch_directory / folder).glob("*")}
        if expected != batches[name]:
            raise ValueError("Native HOG tasks do not cover batching output")
    for item in [trace_record, *[r for task in tasks for r in task["files"]]]:
        check(item)
    return {"status": "fresh_fastoma_task_trace_verified", "trace": trace_record,
            "process_counts": dict(counts), "tasks": tasks,
            "limitations": ["Failed, cached and retried tasks require separate review; no override is inferred.",
                "Task coverage, script presence and wrapper limits, not validation of biological contents.",
                "CPU/memory flags are not measured peak use or aggregate resource accounting.",
                "Published outputs, staged-data identity and native pair semantics require separate checks."]}
