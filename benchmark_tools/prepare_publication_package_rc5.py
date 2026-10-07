"""Inherit immutable rc4 payloads and add the provisional native review, not mutable history."""

import argparse
import json
from pathlib import Path
import re

from benchmark_tools.bundle_publication_package import identity, pin, relative, require, selected


PARENT = "benchmarks/work/orthohmm-study-2026.10.04-rc4"
PARENT_INDEX_SHA = "550d1b06875a26d5719376d5f8562055f68b87ba5d05619b1273146c4c51cf71"
PARENT_SELECTION_SHA = "77784822a38c59257050752a02aa87dc3df48a6c9a3f07a764e06dd4df097283"
DIRECT = "benchmarks/work/native_main_review_component_20261006_v2"
DIRECT_INDEX_SHA = "c8cc4899338689b9b6cd033e6ed7017f2bd3c2e152633a55869923cf7bb5d82a"
READER_SHA = "948efb04abfa84aea1c6822d9145592a15d54416838245f0c1bf569f5d81bbab"


def compose(root, baseline, parent_index, parent_dir, revision, additions):
    root, parent_dir = Path(root).resolve(), Path(parent_dir).resolve()
    require(parent_dir.is_relative_to(root), "Parent package must be within selection root")
    require(re.fullmatch(r"[0-9a-f]{40}", revision), "Require exact workflow revision")
    old = selected(baseline)
    indexed = {relative(r["path"]): r for r in parent_index["files"]}
    require(len(indexed) == len(parent_index["files"])
            and set(indexed) == set(old) | {"PACKAGE_SELECTION.json", "bundle_publication_package.py"},
            "Parent indexed/selected inventories differ")
    require(parent_index["version"] == baseline["version"]
            and parent_index["scientific_revision"] == baseline["scientific_revision"]
            and parent_index["publication_ready"] is False
            and parent_index["native_inference_repeated"] is False,
            "Parent package scope differs")
    for name, row in old.items():
        require(pin(row) == pin(indexed[name]), "Parent selected payload differs from index")
    files = {}

    def add(target, source, expected=None):
        target, source = relative(target), relative(source)
        require(target not in files, "Duplicate rc5 target")
        path = root / source
        require(path.is_file() and not path.is_symlink() and path.resolve().is_relative_to(root),
                "Missing, symlinked or escaped source")
        observed = identity(path)
        if expected is not None:
            require(observed == pin(expected), "Immutable parent/component payload changed")
        files[target] = dict(target=target, source=source, **observed)

    for name, row in indexed.items():
        target = "history/rc4/" + name if name in (
            "README.md", "PACKAGE_SELECTION.json", "bundle_publication_package.py") else name
        add(target, (parent_dir / name).relative_to(root).as_posix(), row)
    add("history/rc4/PACKAGE_INDEX.json", (parent_dir / "PACKAGE_INDEX.json").relative_to(root).as_posix())
    for target, source, expected in additions:
        add(target, source, expected)
    result = dict(baseline, version="orthohmm-study-2026.10.06-rc5", workflow_revision=revision,
        files=list(files.values()), parent_index_sha256=PARENT_INDEX_SHA,
        parent_selection_sha256=PARENT_SELECTION_SHA, preserved_parent_selected_payloads=len(old),
        preserved_parent_indexed_payloads=len(indexed), preserved_parent_payloads=len(indexed),
        native_direct_index_sha256=DIRECT_INDEX_SHA,
        publication_ready=False, public_archive_uploaded=False,
        limitations=[*baseline["limitations"],
            "All inherited rc4 payloads come from its immutable package, never mutable original worktree sources.",
            "native-review/ contains the newer 41-page provisional snapshot; inherited current-review/ remains historical rc4.",
            "native-direct-review/ preserves 99 direct local targets, not complete transitive scientific dependencies or raw data.",
            "The existing shared-host native job/review remain separate; this archive does not admit pending scores or repair failed timings.",
            "No hermetic full-study reproduction, independent-family confirmation, redistribution clearance or public deposition."])
    selected(result)
    return result


def direct_items(root, directory, index):
    root, directory = Path(root).resolve(), Path(directory).resolve()
    require(directory.is_relative_to(root), "Direct component escapes selection root")
    require(index["schema"] == "publication_direct_review_v3"
            and all(index[k] is False for k in (
                "publication_ready", "redistribution_clearance", "transitive_evidence_included")),
            "Wrong direct-review scope")
    rows = {relative(r["path"]): r for r in index["files"]}
    require(len(rows) == len(index["files"]) and "REVIEW_INDEX.json" not in rows
            and all(r["mode"] == 0o644 for r in rows.values()),
            "Outer 0644 package cannot silently change direct component modes")
    items = [("native-direct-review/" + name, (directory / name).relative_to(root).as_posix(), row)
             for name, row in rows.items()]
    items.append(("native-direct-review/REVIEW_INDEX.json",
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
        require(identity(path)["sha256"] == sha, "Changed immutable index, selection or reader")
    baseline = json.loads((parent / "PACKAGE_SELECTION.json").read_text())
    index = json.loads((parent / "PACKAGE_INDEX.json").read_text())
    direct_index = json.loads((direct / "REVIEW_INDEX.json").read_text())
    additions = direct_items(root, direct, direct_index)
    for name, source in (
        ("document.pdf", "publication_main_with_figures_20261006_v2/document.pdf"),
        ("assembly.json", "publication_main_with_figures_20261006_v2/assembly.json"),
        ("figure-selection.json", "publication_figure_selection_20261006_v2.json"),
        ("generation.json", "native_main_text_generation_20261006_v2.json"),
        ("main-visual.json", "native_main_visual_review_20261006_v2.json"),
        ("appendix-visual.json", "native_assembled_visual_review_20261006_v2.json"),
        ("direct-execution.json", "native_main_review_component_execution_20261006_v2.json"),
        ("direct-restoration.json", "native_main_review_component_restoration_20261006_v2.json"),
        ("direct-verifier.strace", "native_main_review_component_verify_20261006_v2.strace"),
    ):
        additions.append(("native-review/" + name, "benchmark_tools/results/" + source, None))
    for source in ("prepare_native_main_text_v2.py", "prepare_native_figure_selection_v2.py",
                   "restore_direct_review_archive.py", "prepare_publication_package_rc5.py"):
        additions.append(("evidence/rc5/workflows/" + source, "benchmark_tools/" + source, None))
    for source in ("test_prepare_native_main_text_v2.py", "test_prepare_native_figure_selection_v2.py",
                   "test_restore_direct_review_archive.py", "test_native_main_review_component_v2.py",
                   "test_prepare_publication_package_rc5.py"):
        additions.append(("evidence/rc5/tests/" + source, "tests/unit/" + source, None))
    additions.extend([
        ("evidence/rc5/PUBLICATION_GOAL_CURRENT.txt", "benchmark_tools/PUBLICATION_GOAL_CURRENT.txt", None),
        ("evidence/rc5/PUBLICATION_PROGRESS.md", "benchmark_tools/results/PUBLICATION_PROGRESS.md", None),
        ("README.md", "benchmark_tools/PUBLICATION_PACKAGE_RC5_20261006.md", None),
    ])
    result = compose(root, baseline, index, parent, revision, additions)
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return dict(path=str(output), **identity(output), files=len(result["files"]),
                preserved_parent_indexed_payloads=result["preserved_parent_indexed_payloads"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--revision", required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.root, args.output, args.revision)))
