"""Export the retained main-text review and direct assets; verify offline."""

import argparse
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path, PurePosixPath
import posixpath
import re
import subprocess
from urllib.parse import unquote, urlsplit

MAIN = "benchmark_tools/results/PUBLICATION_MAIN_TEXT_20260927.md"
RENDER = "benchmark_tools/results/publication_main_render_20260929_v3.json"
PRINT = "benchmark_tools/results/publication_main_print_20260929_v3/print.json"
REVIEW = "benchmark_tools/results/publication_main_pdf_review_20260929_v3/report.json"
LEDGER = "benchmark_tools/results/PUBLICATION_PROGRESS.md"
RUNNER = "benchmark_tools/bundle_publication_review.py"
GUIDE = "benchmark_tools/PUBLICATION_REVIEW_COMPONENT.md"
REVIEW_REVISION = "e2e96746b21d552ea46b9ce63634519f9bae1f82"
LEDGER_REVISION = "fdfe7fb70c583a6ad4ea5c4ae06dd2663d31b07f"


def identity(content):
    return dict(bytes=len(content), sha256=hashlib.sha256(content).hexdigest())


def relative(name):
    if not isinstance(name, str):
        raise ValueError("Require a relative path string")
    path = PurePosixPath(name)
    if (not path.parts or path.is_absolute() or ".." in path.parts
            or str(path) != name or "\\" in name):
        raise ValueError("Unsafe or noncanonical review path")
    return name


def manuscript_path(name):
    name = relative(name)
    if not name.endswith(".md") or name in {LEDGER, RUNNER, GUIDE, "LICENSE.md"}:
        raise ValueError("Require a distinct repository-relative Markdown manuscript")
    return name


def stage_paths(stages=None, main=MAIN):
    if stages is None:
        return dict(render=RENDER, print=PRINT, review=REVIEW)
    if not isinstance(stages, dict) or set(stages) != {"render", "print", "review"}:
        raise ValueError("Require exactly render/print/review stage paths")
    paths = [relative(stages[key]) for key in ("render", "print", "review")]
    if (len(set(paths)) != 3 or any(not path.endswith(".json") for path in paths)
            or set(paths) & {main, LEDGER, RUNNER, GUIDE, "LICENSE.md", "REVIEW_INDEX.json"}):
        raise ValueError("Stage paths must be distinct JSON receipts, not support files")
    return dict(zip(("render", "print", "review"), paths))


def pin(row):
    if (type(row["bytes"]) is not int or row["bytes"] < 0
            or not isinstance(row["sha256"], str)
            or not re.fullmatch(r"[0-9a-f]{64}", row["sha256"])):
        raise ValueError("Invalid review file identity")
    return {key: row[key] for key in ("bytes", "sha256")}


def historical_root(render, main=MAIN):
    source = render["sources"][0]["path"]
    suffix = "/" + main
    if not source.startswith("/") or not source.endswith(suffix):
        raise ValueError("Render source does not identify the main text")
    return source[:-len(suffix)]


def mapped(row, root):
    path = row["path"]
    if not isinstance(path, str) or not path.startswith("/"):
        raise ValueError("Require absolute historical provenance paths")
    return relative(path[len(root) + 1:]) if path.startswith(root + "/") else None


class Links(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.urls = []

    def handle_starttag(self, tag, attrs):
        for key, value in attrs:
            if key in {"href", "src", "data", "poster"} and value:
                self.urls.append(value)


def direct_links(content, html_name):
    parser = Links()
    parser.feed(content.decode("utf-8"))
    result = []
    for url in parser.urls:
        parts = urlsplit(url)
        if parts.scheme or parts.netloc:
            if parts.scheme not in {"http", "https", "mailto"} or parts.netloc and not parts.scheme:
                raise ValueError("Unsupported nonportable HTML URL")
            continue
        if not parts.path:
            continue
        decoded = unquote(parts.path)
        if decoded.startswith("/") or "\\" in decoded:
            raise ValueError("Absolute or platform-dependent local HTML link")
        result.append(relative(posixpath.normpath(posixpath.join(posixpath.dirname(html_name), decoded))))
    return result


def evidence(render, printed, reviewed, main=MAIN):
    if (render["status"] != "manuscript_review_rendered"
            or printed["status"] != "verified_html_printed"
            or reviewed["status"] != "pdf_bounds_checked"
            or any(v["publication_ready"] is not False for v in (render, printed, reviewed))
            or printed["returncode"] != 0 or reviewed["bounds_violations"]
            or type(printed["page_count"]) is not int or printed["page_count"] < 1
            or printed["page_count"] != reviewed["page_count"]):
        raise ValueError("Require retained successful review stages, not publication readiness")
    root = historical_root(render, main)
    records = {}
    external = {}
    for row in [render["html"], *render["sources"], *render["targets"],
                printed["pdf"], *printed["checked_records"],
                *reviewed["checked_records"], *reviewed["rendered_pages"]]:
        name = mapped(row, root)
        table, key = (external, row["path"]) if name is None else (records, name)
        value = pin(row)
        if key in table and table[key] != value:
            raise ValueError("Conflicting historical file identities")
        table[key] = value
    targets = [mapped(row, root) for row in render["targets"]]
    if None in targets or len(set(targets)) != len(targets) or len(targets) != render["unique_targets"]:
        raise ValueError("Invalid direct-target inventory")
    occurrences = render["occurrences"]
    if len(occurrences) != render["local_occurrences"] or {o["path"] for o in occurrences} != set(targets):
        raise ValueError("Direct-target occurrences differ")
    for occurrence in occurrences:
        if direct_links(('<a href="' + occurrence["url"].replace('"', '&quot;') + '">').encode(), main) != [occurrence["path"]]:
            raise ValueError("Historical local link mapping differs")
    return root, records, external, targets


def verify(directory, manifest_sha):
    directory = Path(directory).resolve(strict=True)
    index = directory / "REVIEW_INDEX.json"
    if index.is_symlink() or identity(index.read_bytes())["sha256"] != manifest_sha:
        raise ValueError("Review index differs from external digest")
    manifest = json.loads(index.read_bytes())
    if (manifest["schema"] not in {"publication_direct_review_v1", "publication_direct_review_v2", "publication_direct_review_v3"}
            or manifest["publication_ready"] is not False
            or manifest["redistribution_clearance"] is not False
            or manifest["transitive_evidence_included"] is not False):
        raise ValueError("Review component scope differs")
    main = manuscript_path(manifest.get("main_text")) if manifest["schema"] == "publication_direct_review_v3" else MAIN
    if manifest["schema"] != "publication_direct_review_v3" and "main_text" in manifest:
        raise ValueError("Historical schema cannot override manuscript path")
    if manifest["schema"] in {"publication_direct_review_v2", "publication_direct_review_v3"}:
        if not isinstance(manifest.get("stages"), dict):
            raise ValueError("Missing explicit stage paths")
        stages = stage_paths(manifest["stages"], main)
    else:
        if "stages" in manifest:
            raise ValueError("Historical schema cannot override stage paths")
        stages = stage_paths()
    payloads = {}
    for row in manifest["files"]:
        name = relative(row["path"])
        path = directory / name
        if (name in payloads or path.is_symlink() or not path.is_file()
                or not path.resolve().is_relative_to(directory) or row["mode"] not in (0o644, 0o755)
                or path.stat().st_mode & 0o777 != row["mode"]
                or not re.fullmatch(r"[0-9a-f]{40}", row["git_revision"])
                or not re.fullmatch(r"[0-9a-f]{40}", row["git_blob"])):
            raise ValueError("Invalid review payload mapping, type or mode")
        content = path.read_bytes()
        if identity(content) != pin(row):
            raise ValueError("Review payload identity differs")
        payloads[name] = content
    actual = {p.relative_to(directory).as_posix() for p in directory.rglob("*") if p.is_file() or p.is_symlink()}
    if actual != set(payloads) | {"REVIEW_INDEX.json"}:
        raise ValueError("Extra or missing review payloads")
    support = {RUNNER, GUIDE, "LICENSE.md", *stages.values()}
    if not support <= payloads.keys():
        raise ValueError("Missing review support files")
    render, printed, reviewed = [json.loads(payloads[stages[key]]) for key in ("render", "print", "review")]
    _, expected, external, targets = evidence(render, printed, reviewed, main)
    if set(payloads) != set(expected) | support:
        raise ValueError("Unexpected or missing historical review inventory")
    if manifest["external_provenance_not_included"] != external:
        raise ValueError("External provenance inventory differs")
    for name, row in expected.items():
        if identity(payloads[name]) != row:
            raise ValueError("Historical review identity differs")
    if manifest["schema"] in {"publication_direct_review_v2", "publication_direct_review_v3"}:
        render_name = stages["render"]
        if (render_name not in expected or identity(payloads[render_name]) != expected[render_name]
                or not all(any(mapped(row, historical_root(render, main)) == render_name
                               and pin(row) == identity(payloads[render_name]) for row in stage["checked_records"])
                           for stage in (printed, reviewed))
                or not any(row == printed["pdf"] for row in reviewed["checked_records"])):
            raise ValueError("Selected render/print/review chain differs")
    html_name = mapped(render["html"], historical_root(render, main))
    links = direct_links(payloads[html_name], html_name)
    if set(links) != set(targets):
        raise ValueError("Rendered direct links differ from recorded targets")
    entrypoints = dict(html=html_name, markdown=main, pdf=mapped(printed["pdf"], historical_root(render, main)))
    if manifest["entrypoints"] != entrypoints:
        raise ValueError("Review entrypoints differ")
    return dict(status="publication_direct_review_verified", files=len(payloads),
                payload_bytes=sum(map(len, payloads.values())), direct_targets=len(targets),
                local_html_occurrences=len(links), page_count=printed["page_count"],
                manifest=identity(index.read_bytes()), entrypoints=entrypoints,
                publication_ready=False, transitive_evidence_included=False,
                inference_reproduced=False, redistribution_clearance=False)


def committed(repo, revision, name):
    relative(name)
    entry = subprocess.check_output(["git", "-C", str(repo), "ls-tree", "-z", revision, "--", name])
    if not entry:
        raise FileNotFoundError(name)
    header, recorded = entry.rstrip(b"\0").split(b"\t", 1)
    mode, kind, blob = header.decode().split()
    if recorded.decode() != name or kind != "blob" or mode not in {"100644", "100755"}:
        raise ValueError("Require one regular committed review blob")
    content = subprocess.check_output(["git", "-C", str(repo), "cat-file", "blob", blob])
    return content, int(mode, 8) & 0o777, blob


def build(repo, review_revision, ledger_revision, workflow_revision, output, *, stages=None, main_text=None):
    main = MAIN if main_text is None else manuscript_path(main_text)
    repo, output = Path(repo).resolve(), Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    revisions = [subprocess.check_output(["git", "-C", str(repo), "rev-parse", "--verify", r + "^{commit}"], text=True).strip()
                 for r in (review_revision, ledger_revision, workflow_revision)]
    review_revision, ledger_revision, workflow_revision = revisions
    selected = stage_paths(stages, main)
    render, printed, reviewed = [json.loads(committed(repo, review_revision, selected[key])[0])
                                for key in ("render", "print", "review")]
    root, expected, external, _ = evidence(render, printed, reviewed, main)
    support = {RUNNER, GUIDE, "LICENSE.md", *selected.values()}
    payloads, rows = {}, []
    for name in sorted(set(expected) | support):
        revision = workflow_revision if name in {RUNNER, GUIDE, "LICENSE.md"} else ledger_revision if name == LEDGER else review_revision
        content, mode, blob = committed(repo, revision, name)
        if name in expected and identity(content) != expected[name]:
            raise ValueError("Committed review bytes differ from historical receipt: " + name)
        payloads[name] = content
        rows.append(dict(path=name, mode=mode, git_revision=revision, git_blob=blob, **identity(content)))
    manifest = dict(schema="publication_direct_review_v2" if stages is not None else "publication_direct_review_v1", files=rows,
        entrypoints=dict(html=mapped(render["html"], root), markdown=main, pdf=mapped(printed["pdf"], root)),
        external_provenance_not_included=external, publication_ready=False,
        redistribution_clearance=False, transitive_evidence_included=False,
        limitations=["Only main-text direct links are included; links inside linked documents may be unavailable.",
            "Original receipts retain absolute provenance paths; verification does not read those paths.",
            "The preserved PDF may contain workstation-specific link annotations; use the relative HTML links.",
            "External URLs, fragment anchors, scientific correctness and rights are not certified.",
            "No inference, scoring, plotting, rendering, installation or public release is performed."])
    if stages is not None:
        manifest["stages"] = selected
    if main_text is not None:
        manifest.update(schema="publication_direct_review_v3", main_text=main, stages=selected)
    output.mkdir(parents=True, exist_ok=False)
    for row in rows:
        path = output / row["path"]
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as stream:
            stream.write(payloads[row["path"]])
        path.chmod(row["mode"])
    with (output / "REVIEW_INDEX.json").open("x") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return verify(output, identity((output / "REVIEW_INDEX.json").read_bytes())["sha256"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    builder = commands.add_parser("build")
    builder.add_argument("--repo", type=Path, required=True)
    builder.add_argument("--review-revision", default=REVIEW_REVISION)
    builder.add_argument("--ledger-revision", default=LEDGER_REVISION)
    builder.add_argument("--workflow-revision", required=True)
    builder.add_argument("--output", type=Path, required=True)
    builder.add_argument("--main-text", help="Explicit committed repository-relative Markdown manuscript")
    for role in ("render", "print", "review"):
        builder.add_argument("--" + role + "-receipt", help="Committed repository-relative JSON path; supply all three")
    verifier = commands.add_parser("verify")
    verifier.add_argument("directory", type=Path)
    verifier.add_argument("--manifest-sha256", required=True)
    args = parser.parse_args()
    if args.command == "build":
        stages = {role: getattr(args, role + "_receipt") for role in ("render", "print", "review")}
        if all(value is None for value in stages.values()):
            stages = None
        elif any(value is None for value in stages.values()):
            parser.error("Explicit stage selection requires all three receipts")
        result = build(args.repo, args.review_revision, args.ledger_revision, args.workflow_revision, args.output,
                       stages=stages, main_text=args.main_text)
    else:
        result = verify(args.directory, args.manifest_sha256)
    print(json.dumps(result, indent=2, sort_keys=True))
