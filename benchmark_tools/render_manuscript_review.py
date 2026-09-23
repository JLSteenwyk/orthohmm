"""Render a dated review with checked local assets; not a portable publication archive."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
from urllib.parse import unquote, urlsplit

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def local_assets(document, manuscript, repo):
    occurrences, targets = [], {}

    def visit(node):
        if isinstance(node, list):
            for child in node:
                visit(child)
        elif isinstance(node, dict):
            if node.get("t") in {"Link", "Image"}:
                url = node["c"][2][0]
                parsed = urlsplit(url)
                if not parsed.scheme and not parsed.netloc and parsed.path:
                    relative = Path(unquote(parsed.path))
                    path = (manuscript.parent / relative).resolve()
                    if relative.is_absolute() or not path.is_relative_to(repo):
                        raise ValueError("Local asset escapes repository")
                    if not path.is_file():
                        raise FileNotFoundError(path)
                    name = str(path.relative_to(repo))
                    if name not in targets:
                        targets[name] = record(path)
                    occurrences.append(dict(kind=node["t"], url=url, path=name))
            for child in node.values():
                visit(child)

    visit(document)
    return occurrences, targets


def render(repo, manuscript, output, report):
    repo, manuscript = repo.resolve(), manuscript.resolve()
    output, report = output.absolute(), report.absolute()
    if output.parent.resolve() != manuscript.parent or output == report:
        raise ValueError("Use distinct review paths with HTML beside the manuscript")
    for path in (output, report):
        if path.exists() or path.is_symlink():
            raise FileExistsError(path)
    binary = shutil.which("pandoc")
    if binary is None:
        raise FileNotFoundError("pandoc")
    sources = [record(p) for p in (manuscript, binary, __file__)]
    parse_command = [binary, "--from=markdown", "--to=json", str(manuscript)]
    parsed = subprocess.run(parse_command, capture_output=True, text=True, check=True, timeout=60)
    occurrences, targets = local_assets(json.loads(parsed.stdout), manuscript, repo)
    tracked = set(subprocess.check_output(["git", "ls-files", "-z"], cwd=repo, text=True).split("\0"))
    command = [binary, "--from=markdown", "--to=html5", "--standalone",
               "--metadata=title:OrthoHMM publication working draft", str(manuscript)]
    rendered = subprocess.run(command, capture_output=True, text=True, check=True, timeout=60)
    for item in [*sources, *targets.values()]:
        check(item)
    with output.open("x") as stream:
        stream.write(rendered.stdout)
    result = dict(status="manuscript_review_rendered", publication_ready=False,
        sources=sources, html=record(output), parse_command=parse_command, render_command=command,
        pandoc_version=subprocess.check_output([binary, "--version"], text=True),
        stderr=dict(parse=parsed.stderr, render=rendered.stderr),
        local_occurrences=len(occurrences), unique_targets=len(targets),
        occurrences=occurrences, targets=list(targets.values()),
        untracked_targets=sorted(set(targets) - tracked),
        limitations=["Presence and hashes do not validate scientific contents, anchors, transitive provenance or rights.",
            "Review requires repository-relative assets; not a standalone archive or inference reproduction.",
            "No browser, PDF or visual review is implied by this render receipt.",
            "External URLs and embedded raw HTML are not audited."])
    with report.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "manuscript", "output", "report"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    render(args.repo, args.manuscript, args.output, args.report)
