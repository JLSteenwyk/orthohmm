"""Summarize decoded successful opens; retain other mutations for explicit review."""

import argparse
from collections import Counter
import json
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.probe_dgx_step_separation import save

MUTATIONS = {"mkdir", "mkdirat", "rmdir", "unlink", "unlinkat", "rename", "renameat",
             "renameat2", "link", "linkat", "symlink", "symlinkat", "truncate", "chmod",
             "fchmodat", "chown", "lchown", "fchownat", "utime", "utimes", "utimensat"}


def summarize(lines):
    pending, opened, mutations, executions, unparsed = {}, {}, [], [], []
    nonfilesystem_opens = []
    calls = Counter()
    for number, raw in enumerate(lines, 1):
        match = re.fullmatch(r"(\d+)\s+(.+)", raw.rstrip("\n"))
        if not match:
            unparsed.append(dict(line=number, text=raw.rstrip("\n")))
            continue
        pid, text = int(match[1]), match[2]
        if text.startswith(("+++", "---")):
            continue
        if text.endswith("<unfinished ...>"):
            if pid in pending:
                raise ValueError("Multiple unfinished calls for one PID")
            pending[pid] = (number, text.removesuffix("<unfinished ...>"))
            continue
        resumed = re.match(r"<\.\.\. (\w+) resumed>(.*)", text)
        if resumed:
            if pid not in pending:
                raise ValueError("Resumed syscall without a start")
            first, prefix = pending.pop(pid)
            if not prefix.startswith(resumed[1] + "("):
                raise ValueError("Resumed syscall name differs")
            text = prefix + resumed[2]
            number = first
        call = re.fullmatch(r"(\w+)\((.*)\)\s+=\s+(.+)", text)
        if not call:
            unparsed.append(dict(line=number, pid=pid, text=text))
            continue
        name, arguments, result = call.groups()
        calls[name] += 1
        if result.startswith("-1 "):
            continue
        if name in {"open", "openat", "openat2", "creat"}:
            fd = re.match(r"\d+<([^<>]+)(?:<(?:char|block) \d+:\d+>)?>", result)
            if fd and re.fullmatch(r"pipe:\[\d+\]", fd[1]):
                nonfilesystem_opens.append(dict(line=number, pid=pid, descriptor=fd[1], arguments=arguments))
                continue
            if not fd or not fd[1].startswith("/") or "\\" in fd[1]:
                unparsed.append(dict(line=number, pid=pid, text=text))
                continue
            path = fd[1]
            row = opened.setdefault(path, dict(path=path, successful_opens=0,
                write_capable_opens=0, first_line=number))
            row["successful_opens"] += 1
            row["write_capable_opens"] += int(name == "creat" or bool(re.search(
                r"\b(?:O_WRONLY|O_RDWR|O_CREAT|O_TRUNC|O_APPEND|O_TMPFILE)\b", arguments)))
        elif name in MUTATIONS and re.match(r"0(?:\s|$)", result):
            mutations.append(dict(line=number, pid=pid, syscall=name, arguments=arguments))
        elif name in {"execve", "execveat"} and re.match(r"0(?:\s|$)", result):
            executions.append(dict(line=number, pid=pid, syscall=name, arguments=arguments))
    return dict(status="bounded_file_trace_summary", syscall_counts=dict(calls),
        opened_paths=[opened[k] for k in sorted(opened)], successful_mutations=mutations,
        nonfilesystem_opens=nonfilesystem_opens,
        successful_executions=executions, unparsed=unparsed,
        unfinished={str(k): v for k, v in pending.items()}, scientific_timings_admitted=False,
        limitations=["Successful file opens, not proof of bytes read/written or module execution.",
                     "Other mutation arguments are retained without resolving relative paths.",
                     "Inherited descriptors, memory-only operations and unexercised branches are not covered.",
                     "Unparsed or unfinished records require review; no exhaustive-access claim."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    before = record(args.trace)
    with args.trace.open() as handle:
        result = summarize(handle)
    if record(args.trace) != before:
        raise ValueError("Trace changed during summarization")
    result.update(trace=before, source=record(__file__))
    save(args.output, result)
