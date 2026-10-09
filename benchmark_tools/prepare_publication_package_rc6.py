"""Preserve immutable rc5 payloads and select the terminal-failure reporting revision."""

import argparse
import json
from pathlib import Path
import re

from benchmark_tools.bundle_publication_package import (
    identity, pin, relative, require, save, selected, verify,
)


VERSION = "orthohmm-study-2026.10.09-rc6"
PARENT = "benchmarks/work/orthohmm-study-2026.10.06-rc5"
PARENT_INDEX_SHA = "66bc02d7821d1097808a870c38198bded9c88c379c24052b1828867af8692922"
PARENT_SELECTION_SHA = "3825124e5a4397d713911d236766cb9074f29562ad986eb12c59473f6f5cc724"
DIRECT = "benchmarks/work/native_qfo_terminal_failure_direct_review_20261009_v1"
DIRECT_INDEX_SHA = "506c2a0b64d95282debb4eeae5b111fcaae39dbf8f5ebeb5e96f0319a227e396"
READER_SHA = "948efb04abfa84aea1c6822d9145592a15d54416838245f0c1bf569f5d81bbab"


def compose(root, baseline, parent_index, parent_dir, revision, additions):
    root, parent_dir = Path(root).resolve(), Path(parent_dir).resolve()
    require(parent_dir.is_relative_to(root), "Parent package escapes selection root")
    require(isinstance(revision, str) and re.fullmatch(r"[0-9a-f]{40}", revision),
            "Require exact workflow revision")
    old = selected(baseline)
    indexed = {relative(row["path"]): row for row in parent_index["files"]}
    require(len(indexed) == len(parent_index["files"])
            and set(indexed) == set(old) | {"PACKAGE_SELECTION.json", "bundle_publication_package.py"},
            "Parent indexed/selected inventories differ")
    require(parent_index["schema"] == "publication_package_v1"
            and parent_index["version"] == baseline["version"] == "orthohmm-study-2026.10.06-rc5"
            and all(parent_index[key] == baseline[key]
                    for key in ("scientific_revision", "workflow_revision"))
            and all(parent_index[key] is False for key in (
                "publication_ready", "public_archive_uploaded", "native_inference_repeated"))
            and all(row["mode"] == 0o644 for row in indexed.values()), "Parent package scope differs")
    for name, row in old.items():
        require(pin(row) == pin(indexed[name]), "Parent selected payload differs from index")
    files = {}

    def add(target, source, expected=None):
        target, source = relative(target), relative(source)
        require(target not in files, "Duplicate rc6 target")
        path = root / source
        require(path.is_file() and not path.is_symlink() and path.resolve().is_relative_to(root),
                "Missing, symlinked or escaped source")
        observed = identity(path)
        if expected is not None:
            require(observed == pin(expected), "Immutable parent/component payload changed")
        files[target] = dict(target=target, source=source, **observed)

    for name, row in indexed.items():
        target = "history/rc5/" + name if name in (
            "README.md", "PACKAGE_SELECTION.json", "bundle_publication_package.py") else name
        add(target, (parent_dir / name).relative_to(root).as_posix(), row)
    add("history/rc5/PACKAGE_INDEX.json", (parent_dir / "PACKAGE_INDEX.json").relative_to(root).as_posix())
    for target, source, expected in additions:
        add(target, source, expected)
    result = dict(baseline, version=VERSION, workflow_revision=revision, files=list(files.values()),
        parent_index_sha256=PARENT_INDEX_SHA, parent_selection_sha256=PARENT_SELECTION_SHA,
        preserved_parent_selected_payloads=len(old), preserved_parent_indexed_payloads=len(indexed),
        preserved_parent_payloads=len(indexed), terminal_direct_index_sha256=DIRECT_INDEX_SHA,
        publication_ready=False, public_archive_uploaded=False,
        limitations=[*baseline["limitations"],
            "Inherited rc5 payloads come from its immutable package, never mutable original sources.",
            "terminal-direct-review/ is current 22-page reporting; all other review routes and their status statements are historical.",
            "Four admitted native QfO cells are unchanged; three cells remain missing, including native11 scoring OOM and native12 inference SIGSEGV.",
            "Manual/content/citation and archive-execution receipts are separate evidence/rc6 records, not additions to the strict direct component.",
            "Byte verification does not reproduce inference/scoring, execute nested components, establish independent-family confirmation or clear redistribution rights.",
            "No full-study runtime restoration, public deposition, archival DOI or isolated timing advantage is claimed."])
    selected(result)
    return result


def direct_items(root, directory, index):
    root, directory = Path(root).resolve(), Path(directory).resolve()
    require(directory.is_relative_to(root), "Direct component escapes selection root")
    require(index["schema"] == "publication_direct_review_v3"
            and all(index[key] is False for key in (
                "publication_ready", "redistribution_clearance", "transitive_evidence_included")),
            "Wrong direct-review scope")
    rows = {relative(row["path"]): row for row in index["files"]}
    require(len(rows) == len(index["files"]) and "REVIEW_INDEX.json" not in rows
            and all(row["mode"] == 0o644 for row in rows.values()),
            "Outer package cannot silently change direct component modes")
    items = [("terminal-direct-review/" + name,
              (directory / name).relative_to(root).as_posix(), row) for name, row in rows.items()]
    items.append(("terminal-direct-review/REVIEW_INDEX.json",
                  (directory / "REVIEW_INDEX.json").relative_to(root).as_posix(), None))
    return items


def run(root, output, revision):
    root, output = Path(root).resolve(), Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    parent, direct = root / PARENT, root / DIRECT
    for path, sha in ((parent / "PACKAGE_INDEX.json", PARENT_INDEX_SHA),
                      (parent / "PACKAGE_SELECTION.json", PARENT_SELECTION_SHA),
                      (direct / "REVIEW_INDEX.json", DIRECT_INDEX_SHA),
                      (root / "benchmark_tools/bundle_publication_package.py", READER_SHA)):
        require(not path.is_symlink() and identity(path)["sha256"] == sha,
                "Changed immutable index, selection or reader")
    verify(parent, PARENT_INDEX_SHA)
    baseline = json.loads((parent / "PACKAGE_SELECTION.json").read_text())
    index = json.loads((parent / "PACKAGE_INDEX.json").read_text())
    additions = direct_items(root, direct, json.loads((direct / "REVIEW_INDEX.json").read_text()))
    for name in (
        "NATIVE_QFO_TERMINAL_FAILURE_REVIEW_20261009.md",
        "NATIVE_QFO_TERMINAL_FAILURE_COMPONENT_20261009.md",
        "native_qfo_terminal_failures_20261009_v1_citations.json",
        "native_qfo_terminal_failures_20261009_v1_content_readback.json",
        "native_qfo_terminal_failure_direct_index_20261009_v1.json",
        "native_qfo_terminal_failure_direct_archive_20261009_v1.json",
        "native_qfo_terminal_failure_direct_restore_20261009_v1.json",
        "native_qfo_terminal_failure_direct_verify_20261009_v1.json",
        "native_qfo_terminal_failures_20261009_v1/scores.tsv",
        "native_qfo_terminal_failures_20261009_v1/status.tsv",
        "native_qfo_terminal_failures_20261009_v1/addendum.md",
    ):
        additions.append(("evidence/rc6/results/" + name, "benchmark_tools/results/" + name, None))
    for name in ("export_native_qfo_terminal_failures.py", "review_native12_composed_attempt.py",
                 "render_manuscript_review.py", "print_manuscript_review.py",
                 "audit_manuscript_citations.py", "restore_direct_review_archive.py",
                 "prepare_publication_package_rc6.py"):
        additions.append(("evidence/rc6/workflows/" + name, "benchmark_tools/" + name, None))
    for name in ("test_export_native_qfo_terminal_failures.py", "test_review_native12_composed_attempt.py",
                 "test_render_manuscript_review.py", "test_print_manuscript_review.py",
                 "test_audit_manuscript_citations.py", "test_restore_direct_review_archive.py",
                 "test_prepare_publication_package_rc6.py"):
        additions.append(("evidence/rc6/tests/" + name, "tests/unit/" + name, None))
    additions.extend([
        ("evidence/rc6/PUBLICATION_GOAL_CURRENT.txt", "benchmark_tools/PUBLICATION_GOAL_CURRENT.txt", None),
        ("evidence/rc6/PUBLICATION_PROGRESS.md", "benchmark_tools/results/PUBLICATION_PROGRESS.md", None),
        ("README.md", "benchmark_tools/PUBLICATION_PACKAGE_RC6_20261009.md", None),
    ])
    result = compose(root, baseline, index, parent, revision, additions)
    save(output, result)
    return dict(path=str(output), **identity(output), files=len(result["files"]),
                preserved_parent_indexed_payloads=result["preserved_parent_indexed_payloads"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--revision", required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.root, args.output, args.revision)))
