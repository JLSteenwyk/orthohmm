"""Select completed publication components for a local versioned candidate."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "benchmark_tools/results"
VERSION = "orthohmm-study-2026.10.04-rc1"


def record(path):
    path = Path(path)
    digest, size = hashlib.sha256(), 0
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
            size += len(block)
    return {"source": str(path.relative_to(ROOT)), "bytes": size, "sha256": digest.hexdigest()}


def build(output, revision):
    commit = subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "--verify", revision + "^{commit}"], text=True).strip()
    for name in ("benchmark_tools/bundle_publication_package.py", "benchmark_tools/results/prepare_publication_package_20261004.py",
                 "benchmark_tools/PUBLICATION_PACKAGE_20261004.md"):
        committed = subprocess.check_output(["git", "-C", str(ROOT), "show", commit + ":" + name])
        if committed != (ROOT / name).read_bytes():
            raise ValueError("Package helper/guide differs from selected commit")
    files = {}

    def add(target, path, expected=None):
        ref = record(path)
        if expected is not None and (ref["bytes"], ref["sha256"]) != (expected["bytes"], expected["sha256"]):
            raise ValueError("Retained component changed: " + str(path))
        if target in files and files[target] != dict(target=target, **ref):
            raise ValueError("Conflicting package target")
        files[target] = dict(target=target, **ref)

    handoff = json.loads((RESULTS / "final_publication_handoff_execution_20261004.json").read_text())
    add("archives/handoff.tar.gz", ROOT / handoff["archive"]["path"], handoff["archive"])
    add("anchors/HANDOFF_INDEX.json", RESULTS / "final_publication_handoff_index_20261004.json", handoff["public_index"])
    # The reporting execution records the public archive as a closed output.
    archive_path = RESULTS / "final_prepared_resource_reporting_component_20261004.tar.gz"
    archive_pin = {"bytes": 465108, "sha256": "f48b5a3946bc982f038d65a7a4b2222214c4ef49caacaa8190d48f2a4b11c8a8"}
    add("archives/resource-reporting.tar.gz", archive_path, archive_pin)
    add("anchors/RESOURCE_MANIFEST.json", RESULTS / "final_prepared_resource_reporting_manifest_20261004.json",
        {"bytes": 3997, "sha256": "e36c79c29ed47deec1b32cfded059d92c75726ecc0054fc87c84fc396b551bb8"})
    native = json.loads((RESULTS / "integrated_execution_archive_restoration_20260929.json").read_text())
    add("archives/native-execution-assets.tar", Path(native["archive"]["path"]), native["archive"])
    selection_path = RESULTS / "publication_figure_selection_20261004.json"
    selection = json.loads(selection_path.read_text())
    add("figure-replay/selection.json", selection_path)
    add("figure-replay/assemble.py", ROOT / "benchmark_tools/assemble_publication_review.py")
    for ref in [selection["main"], *selection["provenance"], *[item["pdf"] for item in selection["figures"]]]:
        add("figure-replay/" + ref["path"], ROOT / ref["path"], ref)
    assembly_path = RESULTS / "publication_main_with_figures_20261004/assembly.json"
    assembly = json.loads(assembly_path.read_text())
    add("review/document.pdf", Path(assembly["pdf"]["path"]), assembly["pdf"])
    add("review/assembly.json", assembly_path)
    for ref in assembly["rendered_pages"]:
        add("review/" + Path(ref["path"]).name, Path(ref["path"]), ref)
    evidence = [
        "PUBLICATION_MAIN_TEXT_20261004.md", "PUBLICATION_CLAIMS_20260916.md",
        "PUBLICATION_FINAL_FIGURE_REVIEW_20261004.md", "publication_final_figure_review_20261004.json",
        "PUBLICATION_MAIN_FINAL_REVIEW_20261004.md", "publication_main_final_visual_review_20261004.json",
        "FINAL_PUBLICATION_HANDOFF_20261004.md", "final_publication_handoff_execution_20261004.json",
        "FINAL_RESOURCE_REPORTING_COMPONENT_20261004.md", "final_prepared_resource_reporting_execution_20261004.json",
        "RESTORED_ARCHIVE_FULL_OB_RESULT_22377.md", "restored_archive_full_ob_result_22377.json",
        "RESTORED_ARCHIVE_EXECUTION_PROTOCOL_20260929.md", "INTEGRATED_ASSET_ARCHIVE_20260929.md",
        "integrated_execution_archive_restoration_20260929.json", "orthobench_public_input_chain_20260929.json",
        "THREADRIPPER_ASYNC_CALIBRATION_RESULT_22380.md",
        "PUBLICATION_BASE_ACQUISITION_20261002.md", "SOURCE_ORTHOBENCH_INPUT_SUPPORT_20261002.md",
    ]
    for name in evidence:
        add("evidence/" + name, RESULTS / name)
    add("README.md", ROOT / "benchmark_tools/PUBLICATION_PACKAGE_20261004.md")
    add("LICENSE.md", ROOT / "LICENSE.md")
    add("evidence/package_tests.xml", ROOT / "benchmarks/work/publication_package_tests_20261004.xml")
    add("evidence/package_test_source.py", ROOT / "tests/unit/test_bundle_publication_package.py")
    add("evidence/figure_assembly_tests.xml", ROOT / "benchmarks/work/final_figure_assembly_tests_20261004_v2.xml")
    add("evidence/prepare_publication_package_20261004.py", Path(__file__))
    chosen = {"schema": "publication_package_selection_v1", "version": VERSION,
              "scientific_revision": handoff["scientific_revision"], "workflow_revision": commit,
              "handoff_workflow_revision": handoff["workflow_revision"], "files": list(files.values()),
              "publication_ready": False, "public_archive_uploaded": False,
              "limitations": [
                  "Local versioned working candidate, not a completed scientific requirement audit or deposited release.",
                  "Native asset payloads are delivered, but raw inputs, base interpreter, bootstrap and OS require documented separate preparation.",
                  "Prior archive-to-results reproduction is retained evidence, not a new native execution from this outer package.",
                  "Shared-host timings have unknown, potentially method-dependent contention; no isolated speed ranking.",
                  "Not all transitive scientific raw outputs or challenge-specific uncertainty are bundled or resolved.",
                  "Historical nested guides and paths are preserved; the current outer README governs package routing.",
                  "No publication readiness, public upload, rights clearance, new DOI or journal submission is claimed."]}
    with output.open("x") as stream:
        json.dump(chosen, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps(dict(path=str(output.relative_to(ROOT)), **record(output), files=len(files)), indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--revision", required=True)
    args = parser.parse_args()
    build(args.output.resolve(), args.revision)
