"""Extract retained Memcheck stacks without rerunning or admitting inference."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def summarize(path, expected_sha):
    if digest(path) != expected_sha:
        raise ValueError("Changed Memcheck XML")
    counts, leak_sites, occurrences = Counter(), Counter(), {}
    details, states, seen = [], [], set()
    protocol = None
    fatal_signals = 0
    parser = ET.iterparse(path, events=("start", "end"))
    event, root = next(parser)
    if event != "start" or root.tag != "valgrindoutput":
        raise ValueError("Not Valgrind XML")
    depth = 1
    for event, element in parser:
        if event == "start":
            depth += 1
            continue
        if depth == 2:
            if element.tag == "protocoltool":
                protocol = element.text
            elif element.tag == "status":
                states.append(element.findtext("state"))
            elif element.tag == "fatal_signal":
                fatal_signals += 1
            elif element.tag == "error":
                unique, kind = element.findtext("unique"), element.findtext("kind")
                if not unique or not kind or unique in seen:
                    raise ValueError("Missing or duplicate error identity")
                seen.add(unique)
                counts[kind] += 1
                stacks = [[{child.tag: child.text for child in frame}
                           for frame in stack.findall("frame")]
                          for stack in element.findall("stack")]
                if kind.startswith("Leak_"):
                    first = stacks[0][0] if stacks and stacks[0] else {}
                    leak_sites[(kind, first.get("obj", "unknown"), first.get("fn", "unknown"))] += 1
                else:
                    details.append(dict(unique=unique, kind=kind, what=element.findtext("what"),
                                        auxiliary=[e.text for e in element.findall("auxwhat")], stacks=stacks))
            elif element.tag == "errorcounts":
                for pair in element.findall("pair"):
                    unique, count = pair.findtext("unique"), int(pair.findtext("count"))
                    if not unique or unique in occurrences or count < 1:
                        raise ValueError("Invalid occurrence inventory")
                    occurrences[unique] = count
            root.remove(element)
        depth -= 1
    if protocol != "memcheck" or states != ["RUNNING", "FINISHED"]:
        raise ValueError("Incomplete Memcheck run")
    if set(occurrences) - seen:
        raise ValueError("Occurrence references unknown error")
    if digest(path) != expected_sha:
        raise ValueError("XML changed during extraction")
    return dict(status="retained_memcheck_stacks_extracted", source_sha256=expected_sha,
                source_bytes=path.stat().st_size, states=states, fatal_signal_records=fatal_signals,
                error_records_by_kind=dict(sorted(counts.items())), reported_occurrences=occurrences,
                nonleak_errors=details, leak_top_frames=[dict(kind=k, object=o, function=f, records=n)
                    for (k, o, f), n in sorted(leak_sites.items())],
                accuracy_admitted=False, publication_ready=False,
                limitations=["XML record counts are not event counts or unique allocations.",
                    "Stack locations do not prove the cause of earlier crashes or exclude other defects.",
                    "Leak records are summarized, not dismissed; full stacks remain in pinned XML.",
                    "No native execution, input/runtime revalidation or scientific admission performed."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xml", required=True, type=Path)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    report = summarize(args.xml, args.sha256)
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
