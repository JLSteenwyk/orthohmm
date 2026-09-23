"""Inventory the original trees embedded in the frozen QfO SwissTree reference."""

import argparse
import ast
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess

from benchmark_tools.audit_qfo_swiss_counts import IMAGE_SHA, REFERENCE_SHA
from benchmark_tools.inventory_swisstree_duplications import FAMILIES
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def parse(stdout, expected=None):
    expected = set(FAMILIES.values()) - {None} if expected is None else set(expected)
    cases, slots, current = {}, 0, None
    for line in stdout.splitlines():
        if not line.startswith("RETAINED_"):
            continue
        fields = line.split("\t")
        if len(fields) != 3:
            raise ValueError("Malformed native inventory row")
        kind, family, value = fields
        if kind == "RETAINED_CASE":
            if slots or family in cases or family not in expected or not value.isdigit():
                raise ValueError("Invalid or incomplete family")
            current, slots = family, 1
            cases[family] = dict(mapped_proteins=int(value), leaves=[], node_annotations=[])
            continue
        if family != current or slots <= 0:
            raise ValueError("Invalid preorder traversal")
        if kind == "RETAINED_LEAF":
            if not value:
                raise ValueError("Empty leaf")
            cases[family]["leaves"].append(value)
            slots -= 1
        elif kind == "RETAINED_NODE":
            annotations = ast.literal_eval(value)
            if not isinstance(annotations, list) or any(not isinstance(x, str) for x in annotations):
                raise ValueError("Unexpected native annotation type")
            cases[family]["node_annotations"].append(annotations)
            slots += 1
        else:
            raise ValueError("Unknown native row")
    if slots or set(cases) != expected:
        raise ValueError("Incomplete native tree inventory")
    for case in cases.values():
        nodes = case["node_annotations"]
        case["leaf_count"] = len(case["leaves"])
        case["duplicate_leaf_labels"] = {label: count for label, count in
            sorted(Counter(case["leaves"]).items()) if count > 1}
        case["internal_node_count"] = len(nodes)
        case["explicit_D_Y_nodes"] = sum(any("D=Y" in value.split(":")
            for value in row) for row in nodes)
        case["empty_annotation_nodes"] = sum(not row for row in nodes)
        case["annotation_histogram"] = dict(sorted(Counter(x for row in nodes for x in row).items()))
    return cases


def inventory(repo, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    base = repo / "qfo_benchmark"
    reference = base / "benchmark-webservice/reference_data/2020/ReconciledTrees_SwissTrees.drw"
    image = base / "scoring/container_cache/qfobenchmark-darwin-2022.1.img"
    script = Path(__file__).with_suffix(".drw")
    paths = [reference, image, script,
             base / "benchmark-webservice/lib/RecTree",
             base / "benchmark-webservice/generateData/AddReconciledTree.drw"]
    identities = [record(path) for path in paths]
    if identities[0]["sha256"] != REFERENCE_SHA or identities[1]["sha256"] != IMAGE_SHA:
        raise ValueError("Changed frozen reference or Darwin image")
    if any(c in str(reference) for c in "'\n\r"):
        raise ValueError("Unsupported Darwin path quoting")
    command = ["singularity", "exec", str(image), "darwin", "-E"]
    result = subprocess.run(command, input=f"reference := '{reference}':\n" + script.read_text(),
                            text=True, capture_output=True, check=True, timeout=60)
    cases = parse(result.stdout)
    for case in cases.values():
        del case["leaves"]
        del case["node_annotations"]
    for identity in identities:
        check(identity)
    report = dict(status="retained_native_tree_inventory", source=record(__file__),
        checked_inputs=identities, command=command,
        stdout_sha256=hashlib.sha256(result.stdout.encode()).hexdigest(), stderr=result.stderr,
        families=cases, prediction_statistics_evaluated=False,
        leaf_mapping_validated=False, duplication_strata_admitted=False,
        limitations=["Embedded trees are the retained reference version, not current SwissTree downloads.",
            "MappedProts counts do not provide a leaf-to-benchmark-protein mapping.",
            "D_Y counts include all retained tree taxa, not only benchmark taxa.",
            "Empty annotations are not explicit speciation observations; the native generator defaults to S.",
            "No pair-relation reconstruction, outcome joins, new strata or confidence intervals performed."])
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    inventory(args.repo.resolve(), args.output)
