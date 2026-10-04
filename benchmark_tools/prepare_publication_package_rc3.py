"""Preserve rc2 payload identities and add mechanism text/reporting to rc3."""

import argparse
import json
from pathlib import Path
import re

from benchmark_tools.bundle_publication_package import identity, selected, require


BASELINE_SHA = "593bb4f9298c99da9d3a7428e184a95e832fb73ee1c68759e453e30c144b888e"
EVIDENCE = (
    "PUBLICATION_MAIN_TEXT_20261004_v2.md",
    "SIMULATION_GENE_TREE_ORACLE_PROTOCOL_20261004.md",
    "SIMULATION_GENE_TREE_ORACLE_RESULTS_20261004.md",
    "simulation_gene_tree_oracle_execution_20261004.json",
    "simulation_gene_tree_oracle_readback_20261004.json",
    "SIMULATION_UPSTREAM_TRACE_PROTOCOL_20261004.md",
    "SIMULATION_UPSTREAM_TRACE_RESULTS_20261004.md",
    "SIMULATION_UPSTREAM_TRACE_REVIEW_20261004.md",
    "simulation_upstream_trace_execution_20261004.json",
    "simulation_upstream_trace_readback_20261004.json",
    "SIMULATION_ORACLE_RESIDUAL_PROTOCOL_20261004.md",
    "SIMULATION_ORACLE_RESIDUAL_RESULTS_20261004.md",
    "SIMULATION_ORACLE_RESIDUAL_REVIEW_20261004.md",
    "simulation_oracle_residual_execution_20261004.json",
    "simulation_oracle_residual_readback_20261004.json",
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
            require(observed == expected, "Inherited rc2 payload changed")
        files[target] = {"target": target, "source": source, **observed}
    for target, row in previous.items():
        preserved = "evidence/PUBLICATION_PACKAGE_RC2_20261004.md" if target == "README.md" else target
        add(preserved, row["source"], {k: row[k] for k in ("bytes", "sha256")})
    for target, source in additions:
        add(target, source)
    result = {**baseline, "version": "orthohmm-study-2026.10.04-rc3", "workflow_revision": revision,
              "files": list(files.values()), "parent_selection_sha256": BASELINE_SHA,
              "preserved_parent_payloads": len(previous),
              "limitations": [*baseline["limitations"],
                  "Mechanism replay uses checked summaries, not raw native inputs or independent biological validation.",
                  "Historical review PDF predates the new mechanism-integrated Markdown; new graphical review remains pending."]}
    selected(result)
    return result


def run(root, output, revision):
    require(not output.exists() and not output.is_symlink(), "Refusing existing selection")
    path = root / "benchmark_tools/results/publication_package_rc2_selection_20261004_v2.json"
    require(identity(path)["sha256"] == BASELINE_SHA, "Changed parent selection")
    baseline = json.loads(path.read_text())
    additions = [("evidence/" + name, "benchmark_tools/results/" + name) for name in EVIDENCE]
    additions.extend([
        ("evidence/replay_simulation_mechanisms.py", "benchmark_tools/replay_simulation_mechanisms.py"),
        ("evidence/prepare_publication_package_rc3.py", "benchmark_tools/prepare_publication_package_rc3.py"),
        ("evidence/test_replay_simulation_mechanisms.py", "tests/unit/test_replay_simulation_mechanisms.py"),
        ("evidence/test_prepare_publication_package_rc3.py", "tests/unit/test_prepare_publication_package_rc3.py"),
        ("evidence/mechanism_reporting_tests.xml", "benchmarks/work/simulation_mechanism_reporting_corrected_tests_20261004.xml"),
        ("evidence/mechanism_package_tests.xml", "benchmarks/work/simulation_mechanism_package_tests_20261004.xml"),
        ("README.md", "benchmark_tools/PUBLICATION_PACKAGE_RC3_20261004.md"),
    ])
    result = compose(root, baseline, revision, additions)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x") as handle:
        handle.write(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return {"path": str(output), **identity(output), "files": len(result["files"]),
            "preserved_rc2_payloads": result["preserved_parent_payloads"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--revision", required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.root.resolve(), args.output.absolute(), args.revision), sort_keys=True))
