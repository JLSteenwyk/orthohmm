"""Pin a human-reviewed requirement register without certifying completion."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re


STATUSES = {"supported", "partial", "unmet"}
HEADER = re.compile(r"^## ((?:scope|engineering|completion|[1-7])\.[1-9][0-9]*) \| (\w+)$", re.M)


def requirements(text):
    lines = []
    for line in text.splitlines():
        line = line.strip()
        line = line[1:].lstrip() if line.startswith(">") else line
        if re.fullmatch(r"[1-7]\. .+", line) or line in {
            "Engineering and reporting requirements", "Completion criteria"
        }:
            lines.extend(["", line, ""])
        elif line.startswith("- "):
            lines.extend(["", line])
        else:
            lines.append(line)
    paragraphs = [" ".join(p.splitlines()) for p in re.split(r"\n\s*\n", "\n".join(lines)) if p.strip()]
    if not paragraphs or not paragraphs[0].startswith("Goal:"):
        raise ValueError("Expected authoritative Goal title")
    rows, counts, section, headings = [], Counter(), "scope", []
    for paragraph in paragraphs[1:]:
        heading = re.fullmatch(r"([1-7])\. .+", paragraph)
        if heading:
            section = heading[1]
            headings.append(section)
            continue
        if paragraph in {"Engineering and reporting requirements", "Completion criteria"}:
            section = "engineering" if paragraph.startswith("Engineering") else "completion"
            headings.append(section)
            continue
        if section not in {"scope", "completion"} and not paragraph.startswith("- "):
            raise ValueError("Unrecognized requirement paragraph: " + paragraph)
        counts[section] += 1
        rows.append({"id": f"{section}.{counts[section]}",
                     "requirement": paragraph[2:] if paragraph.startswith("- ") else paragraph})
    if headings != [str(i) for i in range(1, 8)] + ["engineering", "completion"]:
        raise ValueError("Missing, duplicated or reordered goal sections")
    if any(counts[s] == 0 for s in ["scope", "engineering", "completion", *map(str, range(1, 8))]):
        raise ValueError("Empty goal section")
    return rows


def reviewed(text):
    matches = list(HEADER.finditer(text))
    rows = {}
    for i, match in enumerate(matches):
        identifier, status = match.groups()
        body = text[match.end():matches[i + 1].start() if i + 1 < len(matches) else len(text)].strip()
        if identifier in rows or status not in STATUSES or not body:
            raise ValueError("Invalid or duplicate reviewed requirement: " + identifier)
        paths = re.findall(r"\[[^\]]+\]\(([^\s)]+)\)", body)
        if not paths or any("://" in p or Path(p).is_absolute() for p in paths):
            raise ValueError("Each review needs relative retained evidence links")
        rows[identifier] = {"status": status, "assessment": body, "evidence_paths": sorted(set(paths))}
    return rows


def identity(path):
    data = path.read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def build(goal, review, output):
    goal, review, output = Path(goal), Path(review), Path(output)
    goal_bytes, review_bytes = goal.read_bytes(), review.read_bytes()
    original, decisions = requirements(goal_bytes.decode()), reviewed(review_bytes.decode())
    if set(decisions) != {r["id"] for r in original}:
        raise ValueError("Reviewed inventory does not exactly cover authoritative requirements")
    evidence = {}
    for decision in decisions.values():
        for name in decision["evidence_paths"]:
            path = (review.parent / name).resolve(strict=True)
            if not path.is_relative_to(review.parent.resolve()) or not path.is_file():
                raise ValueError("Evidence must be a retained file beneath the register directory")
            evidence[name] = identity(path)
    rows = [{**r, **decisions[r["id"]]} for r in original]
    result = {
        "schema_version": 1,
        "status": "publication_requirement_audit_incomplete",
        "goal_completion_proven": False,
        "scope": "Human-reviewed direct evidence; byte pins do not re-admit raw scientific results or certify biological/statistical assumptions.",
        "goal": {"snapshot": "goal.txt", "bytes": len(goal_bytes), "sha256": hashlib.sha256(goal_bytes).hexdigest()},
        "review": identity(review), "generator": identity(Path(__file__)),
        "requirement_count": len(rows), "status_counts": dict(sorted(Counter(r["status"] for r in rows).items())),
        "requirements": rows, "evidence": dict(sorted(evidence.items())),
    }
    # Check inputs again before emitting a dated audit; never alter retained evidence.
    if goal.read_bytes() != goal_bytes or review.read_bytes() != review_bytes:
        raise ValueError("Goal or review changed during audit")
    for name, pin in evidence.items():
        if identity(review.parent / name) != pin:
            raise ValueError("Evidence changed during audit: " + name)
    output.mkdir(parents=True, exist_ok=False)
    (output / "goal.txt").write_bytes(goal_bytes)
    (output / "audit.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--goal", type=Path, required=True)
    parser.add_argument("--review", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = build(args.goal, args.review, args.output)
    print(json.dumps({k: result[k] for k in ("status", "goal_completion_proven", "requirement_count", "status_counts")}))


if __name__ == "__main__":
    main()
