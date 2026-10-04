"""Preserve rc3 payload bytes and add reviewed exposure/reproducibility reporting."""

import argparse
import json
from pathlib import Path
import re

from benchmark_tools.bundle_publication_package import identity, selected, require


BASELINE_SHA = "1d10dbf8f65c390b80a5470092aa8c72dec2db73a1f625001b53574c939d57a8"
READER_SHA = "948efb04abfa84aea1c6822d9145592a15d54416838245f0c1bf569f5d81bbab"
EVIDENCE = (
    "PUBLICATION_MAIN_TEXT_20261004_v3.md", "PUBLICATION_EVIDENCE_INTEGRATION_RESULT_20261004.md",
    "publication_evidence_integration_execution_20261004.json", "PUBLICATION_EVIDENCE_REVIEW_RESULT_20261004.md",
    "publication_main_evidence_review_execution_20261004.json", "publication_main_evidence_visual_review_20261004.json",
    "DEVELOPMENT_FAMILY_INVENTORY_PROTOCOL_20261004.md", "DEVELOPMENT_FAMILY_INVENTORY_RESULT_20261004.md",
    "development_family_inventory_execution_20261004.json", "development_family_local_manifest_20261004.json",
    "CANDIDATE_TRACE_VARIATION_RESULT_20261004.md", "candidate_trace_variation_execution_20261004.json",
    "FACTORIAL_NATIVE_RESOURCE_LINKAGE_RESULT_20261004.md", "factorial_native_resource_linkage_execution_20261004.json",
    "QFO_ORTHOHMM_STAGE_LINKAGE_RESULT_20261004.md", "qfo_orthohmm_stage_linkage_execution_20261004.json",
    "FACTORIAL_RESOURCE_RESULT_20261004.md", "ALL_BENCHMARK_METADATA_INTEGRATION_RESULT_20261004.md",
)
REPORTS = (
    ("inventory.json", "development_family_inventory_20261004/inventory.json"),
    ("diagnostic.json", "candidate_trace_variation_20261004/diagnostic.json"),
    ("linkage.json", "factorial_native_resource_linkage_20261004/linkage.json"),
    ("qfo_register.json", "qfo_orthohmm_stage_metadata_20261004/register.json"),
    ("factorial_resources.json", "factorial_retained_resources_20261004/resources.json"),
    ("all_tool_register.json", "all_benchmark_metadata_integrated_20261004_v3/register.json"),
    ("family_evidence.tsv", "development_family_inventory_20261004/family_evidence.tsv"),
)


def compose(root, baseline, revision, additions):
    require(re.fullmatch(r"[0-9a-f]{40}", revision), "Require exact workflow revision")
    previous = selected(baseline)
    files = {}
    def add(target, source, expected=None):
        require(target not in files, "Duplicate extension target")
        path = root / source
        require(path.is_file() and path.resolve().is_relative_to(root.resolve()), "Missing or escaped selected source")
        observed = identity(path)
        if expected is not None:
            require(observed == expected, "Inherited rc3 payload changed")
        files[target] = {"target": target, "source": source, **observed}
    for target, row in previous.items():
        kept = "evidence/PUBLICATION_PACKAGE_RC3_20261004.md" if target == "README.md" else target
        add(kept, row["source"], {k: row[k] for k in ("bytes", "sha256")})
    for target, source in additions:
        add(target, source)
    result = {**baseline, "version": "orthohmm-study-2026.10.04-rc4", "workflow_revision": revision,
              "files": list(files.values()), "parent_selection_sha256": BASELINE_SHA,
              "preserved_parent_payloads": len(previous),
              "preserved_parent_selected_payloads": len(previous),
              "limitations": [*baseline["limitations"],
                  "The current-review PDF renders the third source; inherited rc3 graphical-review-pending statements are historical, not current status.",
                  "Addenda replay validates retained metadata projections and arithmetic, not raw Git/local report scans, native traces, execution accounting or scientific independence.",
                  "The 114 local-only historical report inputs are pinned but not included; full exposure discovery requires them and the original Git snapshot."]}
    selected(result)
    return result


def run(root, output, revision):
    require(not output.exists() and not output.is_symlink(), "Refusing existing selection")
    path = root / "benchmark_tools/results/publication_package_rc3_selection_20261004.json"
    require(identity(path)["sha256"] == BASELINE_SHA, "Changed parent selection")
    require(identity(root / "benchmark_tools/bundle_publication_package.py")["sha256"] == READER_SHA, "Changed retained rc3 reader")
    baseline = json.loads(path.read_text())
    additions = [("evidence/" + name, "benchmark_tools/results/" + name) for name in EVIDENCE]
    additions += [("addenda-replay/" + target, "benchmark_tools/results/" + source) for target, source in REPORTS]
    additions += [("current-review/" + target, "benchmark_tools/results/" + source) for target, source in (
        ("document.pdf", "publication_main_with_figures_20261004_v3/document.pdf"),
        ("main.pdf", "publication_main_print_20261004_v3/document.pdf"),
        ("main.html", "publication_main_review_20261004_v3.html"),
        ("render.json", "publication_main_render_20261004_v3.json"),
        ("print.json", "publication_main_print_20261004_v3/print.json"),
        ("bounds.json", "publication_main_pdf_review_20261004_v3/report.json"),
        ("selection.json", "publication_figure_selection_20261004_v3.json"),
        ("assembly.json", "publication_main_with_figures_20261004_v3/assembly.json"),
    )]
    additions += [
        ("addenda-replay/replay_publication_addenda.py", "benchmark_tools/replay_publication_addenda.py"),
        ("evidence/prepare_publication_package_rc4.py", "benchmark_tools/prepare_publication_package_rc4.py"),
        ("evidence/test_replay_publication_addenda.py", "tests/unit/test_replay_publication_addenda.py"),
        ("evidence/test_prepare_publication_package_rc4.py", "tests/unit/test_prepare_publication_package_rc4.py"),
        ("evidence/test_publication_evidence_integration.py", "tests/unit/test_publication_evidence_integration.py"),
        ("evidence/test_evidence_main_review_artifacts.py", "tests/unit/test_evidence_main_review_artifacts.py"),
        ("history/rc3/PACKAGE_SELECTION.json", "benchmark_tools/results/publication_package_rc3_selection_20261004.json"),
        ("history/rc3/PACKAGE_INDEX.json", "benchmark_tools/results/publication_package_rc3_index_20261004.json"),
        ("history/rc3/bundle_publication_package.py", "benchmark_tools/bundle_publication_package.py"),
        ("evidence/addenda_prepared_tests.xml", "benchmarks/work/publication_addenda_prepared_tests_20261004.xml"),
        ("README.md", "benchmark_tools/PUBLICATION_PACKAGE_RC4_20261004.md"),
    ]
    result = compose(root, baseline, revision, additions)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as stream:
        stream.write(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return {"path": str(output), **identity(output), "files": len(result["files"]),
            "preserved_rc3_selected_payloads": len(baseline["files"])}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--revision", required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.root.resolve(), args.output.absolute(), args.revision), sort_keys=True))
