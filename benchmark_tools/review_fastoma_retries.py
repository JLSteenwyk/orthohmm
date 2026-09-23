"""Review explicitly selected FastOMA retries without admitting native outputs."""

import csv
import io
import shlex

from benchmark_tools.audit_fastoma_tasks import option, task_directory, validate_trace, wrapper_limits
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def partition_attempts(text, species_files, retry_pairs):
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if not rows or not {"task_id", "hash", "name", "status", "exit"}.issubset(rows[0]):
        raise ValueError("Missing trace columns or attempts")
    ids = [r["task_id"] for r in rows]
    hashes = [r["hash"] for r in rows]
    if (any(not value.isdigit() for value in ids) or len(set(ids)) != len(ids)
            or len(set(hashes)) != len(hashes)):
        raise ValueError("Duplicate or invalid attempt identity")
    if (not retry_pairs or len(set(retry_pairs.values())) != len(retry_pairs)
            or set(retry_pairs) & set(retry_pairs.values())):
        raise ValueError("Require distinct explicitly reviewed retry pairs")
    failed = {r["hash"]: r for r in rows if r["status"] == "FAILED" and r["exit"] == "1"}
    successful = [r for r in rows if r["status"] == "COMPLETED" and r["exit"] == "0"]
    if set(failed) != set(retry_pairs) or len(failed) + len(successful) != len(rows):
        raise ValueError("Unreviewed failure, cache, or nonterminal attempt")
    successes = {r["hash"]: r for r in successful}
    names = [r["name"] for r in successful]
    if len(set(names)) != len(names):
        raise ValueError("Duplicate successful logical task")
    for old, new in retry_pairs.items():
        before, after = failed[old], successes.get(new)
        if (after is None or before["name"] != after["name"]
                or before["name"].split(" (", 1)[0] not in {"hog_big", "hog_rest"}
                or int(before["task_id"]) >= int(after["task_id"])):
            raise ValueError("Retry is not the later successful matching HOG task")
    # Reuse the strict process-coverage validator; original rows remain in the report.
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(successful)
    _, counts = validate_trace(stream.getvalue(), species_files)
    return rows, failed, successes, counts


def review(trace, work, species_files, retry_pairs):
    trace_record = record(trace)
    rows, failed, successes, counts = partition_attempts(trace.read_text(), species_files, retry_pairs)
    reviewed, checked = [], [trace_record, record(__file__)]
    for old, new in retry_pairs.items():
        attempts = []
        scripts = []
        trees = []
        for row, exit_code in ((failed[old], "1"), (successes[new], "0")):
            directory = task_directory(work, row["hash"])
            files = [record(directory / name) for name in (
                ".command.sh", ".command.run", ".exitcode", ".command.log", ".command.err")]
            if (directory / ".exitcode").read_text().strip() != exit_code:
                raise ValueError("Trace and native exit code disagree")
            script = (directory / ".command.sh").read_bytes()
            argv = shlex.split(script.decode(), comments=True)
            if not argv or argv[0] != "fastoma-infer-subhogs":
                raise ValueError("Unexpected retry executable")
            tree = record(directory / option(argv, "--species-tree"))
            trees.append(tree)
            scripts.append(script)
            limits = wrapper_limits((directory / ".command.run").read_text())
            attempts.append({"trace": row, "directory": str(directory), "files": files,
                             "tree": tree, "limits": limits})
            checked.extend([*files, tree])
        if scripts[0] != scripts[1] or trees[0] != trees[1]:
            raise ValueError("Retry changed scientific command or checked species tree")
        first, second = (a["limits"] for a in attempts)
        if first["cpus"] != second["cpus"] or second["memory_bytes"] != 2 * first["memory_bytes"]:
            raise ValueError("Retry resource change differs from reviewed doubling")
        reviewed.append({"failed": attempts[0], "successful": attempts[1]})
    for item in checked:
        check(item)
    return {"status": "explicit_fastoma_retry_attempts_reviewed_not_admitted",
            "trace": trace_record, "attempts": rows, "successful_process_counts": dict(counts),
            "retry_pairs": reviewed, "checked_records": checked,
            "native_outputs_admitted": False, "accuracy_evaluated": False,
            "limitations": ["Current command/tree identity does not prove historical input-byte immutability.",
                "Batch input identity, container mount coverage, and collection require separate validation.",
                "Resource flags are limits, not measured use or proof of OOM.",
                "Frozen retry-policy provenance and published native outputs remain separate checks."]}
