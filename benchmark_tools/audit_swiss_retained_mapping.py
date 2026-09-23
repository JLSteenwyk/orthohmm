"""Test retained-tree identifier mappings against native reference relations."""

import argparse
import ast
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess

from benchmark_tools.audit_qfo_swiss_counts import IMAGE_SHA, REFERENCE_SHA
from benchmark_tools.inventory_swiss_retained_trees import parse
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

MAPPING_SHA = "1c10f6ce5e53ebc3148dde02d16268b225c3c8817c952d1274daa41acbf9eb4d"
SPEC = (":D=F", ":D=N", "D=N", ":Ev=speciation", ":Ev=SPECIATION")
DUPL = (":D=T", ":D=Y", "D=Y", ":Ev=duplication", ":Ev=DUPLICATION")


def event(annotations):
    # The retained generator defaults absent/bootstrap-only annotations to S.
    for value in annotations:
        if any(token in value for token in SPEC):
            return "S"
        if any(token in value for token in DUPL):
            return "D"
        if re.fullmatch(r"[0-9.]+", value):
            continue
        raise ValueError("Unsupported annotation: " + value)
    return "S"


def reconstruct(rows, mapping):
    iterator = iter(rows)
    relations, collisions, mapped_labels, absent = {}, [], {}, []

    def walk():
        kind, value = next(iterator)
        if kind == "RETAINED_LEAF":
            candidates = [value, *value.split("_")]
            match = next((mapping[label] for label in candidates if label in mapping), None)
            if match is None:
                absent.append(value)
                return set()
            if type(match) is not int or match <= 0:
                raise ValueError("Invalid mapped protein number")
            mapped_labels[value] = match
            return {match}
        annotations = ast.literal_eval(value)
        left, right = walk(), walk()
        if not left:
            return right
        if not right:
            return left
        overlap = left & right
        if overlap:
            collisions.append(sorted(overlap))
            left -= right
        ev = event(annotations)
        for a in sorted(left):
            for b in sorted(right):
                relations[tuple(sorted((a, b)))] = ev
        return left | right

    members = walk()
    if next(iterator, None) is not None:
        raise ValueError("Unconsumed native tree nodes")
    return members, relations, dict(mapped_labels=mapped_labels,
        unmapped_leaf_occurrences=len(absent), intersecting_child_mappings=collisions)


def compare(stdout, mapping, expected=None):
    inventory = parse(stdout, expected)
    trees = {name: [] for name in inventory}
    members = {name: set() for name in inventory}
    truth = {name: {} for name in inventory}
    for line in stdout.splitlines():
        if line.startswith(("RETAINED_LEAF\t", "RETAINED_NODE\t")):
            kind, family, value = line.split("\t")
            trees[family].append((kind, value))
        elif line.startswith("SWISS_MEMBER\t"):
            _, family, number = line.split("\t")
            number = int(number)
            if number <= 0 or number in members[family]:
                raise ValueError("Invalid or duplicate reference member")
            members[family].add(number)
        elif line.startswith("SWISS_RELATION\t"):
            _, family, a, b, value = line.split("\t")
            key, ev = (int(a), int(b)), ast.literal_eval(value)
            if key[0] >= key[1] or key in truth[family] or ev not in {"D", "S"}:
                raise ValueError("Unsupported native reference relation")
            truth[family][key] = ev
    result = {}
    for family, rows in trees.items():
        if len(members[family]) != inventory[family]["mapped_proteins"] or not truth[family]:
            raise ValueError("Incomplete native member/relation export")
        if any(not set(pair) <= members[family] for pair in truth[family]):
            raise ValueError("Reference pair outside family")
        observed, relations, details = reconstruct(rows, mapping)
        missing, extra = truth[family].keys() - relations.keys(), relations.keys() - truth[family].keys()
        changed = [pair for pair in truth[family].keys() & relations.keys()
                   if truth[family][pair] != relations[pair]]
        result[family] = dict(reference_members=len(members[family]), mapped_members=len(observed),
            missing_members=sorted(members[family] - observed), extra_members=sorted(observed - members[family]),
            reference_pairs=len(truth[family]), reconstructed_pairs=len(relations),
            missing_pairs=len(missing), extra_pairs=len(extra), changed_events=len(changed),
            reference_event_counts=dict(Counter(truth[family].values())),
            reconstructed_event_counts=dict(Counter(relations.values())),
            exact_match=observed == members[family] and relations == truth[family], **details)
    return result


def audit(repo, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    refdir = repo / "qfo_benchmark/benchmark-webservice/reference_data/2020"
    image = repo / "qfo_benchmark/scoring/container_cache/qfobenchmark-darwin-2022.1.img"
    reference = refdir / "ReconciledTrees_SwissTrees.drw"
    mapping_path = refdir / "mapping.json.gz"
    scripts = [Path(__file__).with_name(name) for name in
               ("export_swiss_retained_relations.drw", "inventory_swiss_retained_trees.drw")]
    records = [record(path) for path in [reference, image, mapping_path, *scripts,
        Path(__file__).with_name("inventory_swiss_retained_trees.py"),
        repo / "qfo_benchmark/benchmark-webservice/generateData/AddReconciledTree.drw"]]
    if records[0]["sha256"] != REFERENCE_SHA or records[1]["sha256"] != IMAGE_SHA:
        raise ValueError("Changed frozen reference or container")
    if records[2]["sha256"] != MAPPING_SHA:
        raise ValueError("Changed frozen identifier mapping")
    if any(c in str(reference) for c in "'\n\r"):
        raise ValueError("Unsupported Darwin quoting")
    command = ["singularity", "exec", str(image), "darwin", "-E"]
    native = subprocess.run(command, input=f"reference := '{reference}':\n" +
        "\n".join(p.read_text() for p in scripts), text=True, capture_output=True, check=True, timeout=60)
    with gzip.open(mapping_path, "rt") as stream:
        mapping = json.load(stream)["mapping"]
    families = compare(native.stdout, mapping)
    for item in records:
        check(item)
    report = dict(status="retained_mapping_diagnostic", source=record(__file__), checked_inputs=records,
        command=command, stdout_sha256=hashlib.sha256(native.stdout.encode()).hexdigest(),
        stderr=native.stderr, families=families,
        exact_families=sum(row["exact_match"] for row in families.values()),
        prediction_statistics_evaluated=False, duplication_strata_admitted=False,
        limitations=["Reconstruction tests the retained mapping.json, not a recovered original IDIndex.db.",
            "Exact relation agreement supports mapping semantics but is not independent validation of curated biology.",
            "Native generator's default S and left-minus-right duplicate handling are preserved.",
            "No scoring joins or outcome-defined feature thresholds; mismatches remain explicit."])
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit(args.repo.resolve(), args.output)
