"""Render and inspect selected CSL metadata using Pandoc's citation processor."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import unicodedata

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def normalize(value):
    return " ".join(unicodedata.normalize("NFC", value).split())


def text(node):
    if isinstance(node, dict):
        if node.get("t") == "Str":
            return node["c"]
        if node.get("t") in ("Space", "SoftBreak", "LineBreak"):
            return " "
        if node.get("t") == "Link":
            return text(node["c"][1])
        if node.get("t") in ("Para", "Plain"):
            return text(node["c"]) + " "
        return text(node.get("c", []))
    if isinstance(node, list):
        return "".join(text(v) for v in node)
    return ""


def entries(document):
    found = {}
    def visit(node):
        if isinstance(node, dict):
            if node.get("t") == "Div":
                attributes, blocks = node["c"]
                identity, classes, _ = attributes
                if "csl-entry" in classes:
                    if not identity.startswith("ref-") or identity[4:] in found:
                        raise ValueError("Missing or duplicate rendered citation identity")
                    found[identity[4:]] = normalize(text(blocks))
            for value in node.values():
                visit(value)
        elif isinstance(node, list):
            for value in node:
                visit(value)
    visit(document["blocks"])
    return found


def audit(records, document):
    rendered = entries(document)
    if len({r["id"] for r in records}) != len(records) or set(rendered) != {r["id"] for r in records}:
        raise ValueError("Rendered bibliography does not match selected citation inventory")
    checks = []
    for row in records:
        value, expected = rendered[row["id"]], {}
        for key in ("title", "genre", "DOI", "URL", "container-title"):
            if row.get(key):
                expected[key] = normalize(row[key])
        if row.get("number"):
            expected["number"] = "Identifier: " + str(row["number"])
        for i, author in enumerate(row.get("author", [])):
            for key in ("literal", "family", "suffix"):
                if author.get(key):
                    expected[f"author_{i}_{key}"] = normalize(author[key])
        if "issued" in row:
            expected["issued_year"] = str(row["issued"]["date-parts"][0][0])
        if "accessed" in row:
            expected["accessed"] = "Accessed: " + "-".join(
                str(part) if i == 0 else f"{part:02d}" for i, part in enumerate(row["accessed"]["date-parts"][0]))
        missing = {key: wanted for key, wanted in expected.items() if wanted not in value}
        checks.append(dict(id=row["id"], text=value, expected=expected, missing=missing))
    return dict(status="bibliography_render_fields_checked", entries=len(rendered), checks=checks,
                all_checked_fields_present=all(not r["missing"] for r in checks))


def render(bibliography, style, output, executable="pandoc"):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    binary = shutil.which(executable)
    if binary is None:
        raise FileNotFoundError(executable)
    evidence = [record(bibliography), record(style), record(binary), record(__file__)]
    records = json.loads(bibliography.read_text())
    version = subprocess.check_output([binary, "--version"], text=True)
    command = [binary, "--from=markdown", "--to=json", "--citeproc",
               "--bibliography=" + str(bibliography), "--csl=" + str(style)]
    result = subprocess.run(command, input="---\nnocite: '@*'\n---\n", text=True,
                            capture_output=True, check=True, timeout=60)
    document = json.loads(result.stdout)
    checked = audit(records, document)
    html_command = [binary, "--from=json", "--to=html5", "--standalone", "--metadata=title:Bibliography review"]
    html = subprocess.run(html_command, input=result.stdout, text=True, capture_output=True, check=True, timeout=60)
    for item in evidence:
        check(item)
    output.mkdir()
    (output / "bibliography.html").write_text(html.stdout)
    (output / "bibliography.pandoc.json").write_text(result.stdout)
    manifest = dict(status="selected_bibliography_rendered_for_review", evidence=evidence, pandoc_version=version,
        commands=[command, html_command], stderr=[result.stderr, html.stderr], field_audit=checked,
        outputs=[record(output / "bibliography.html"), record(output / "bibliography.pandoc.json")],
        publication_ready=False,
        limitations=["Project review style, not a selected journal's submission style.",
            "Checks field visibility and entry inventory, not semantic citation correctness or visual layout.",
            "Page-range abbreviation and author given-name typography are not validated by substring checks.",
            "Citation coverage, raw-source rights, and scientific evidence are separate requirements."])
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("bibliography", "style", "output"):
        parser.add_argument("--" + key, type=Path, required=True)
    args = parser.parse_args()
    result = render(args.bibliography.resolve(), args.style.resolve(), args.output.absolute())
    raise SystemExit(0 if result["field_audit"]["all_checked_fields_present"] else 1)
