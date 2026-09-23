"""Extract frozen, prediction-independent mapped-tree duplication features."""

import argparse
import ast
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import statistics
import subprocess

from benchmark_tools.audit_swiss_retained_mapping import compare, event, SPEC, DUPL
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen

ADMISSION_SHA = "0d3e736a350782609c68764bc19387500ea64bd5045c5fcfb41930dd89c1ce9d"
HELPER_SHA = "871acee01b9e21441483230076cdadbc4d0d76df771f0fcfde5e3168683029ae"
PROTOCOL_SHA = "2bca0f2b2a0c3a5426be268484bab8302b404d66cfdb29c5872afefd41c79b2f"


def node_class(annotations):
    label = event(annotations)
    explicit = any(any(token in value for token in (*SPEC, *DUPL)) for value in annotations)
    return "explicit_duplication_nodes" if label == "D" else (
        "explicit_speciation_nodes" if explicit else "default_speciation_nodes")


def features(rows, mapping):
    iterator = iter(rows)
    counts = dict(explicit_duplication_nodes=0, explicit_speciation_nodes=0,
                  default_speciation_nodes=0, child_overlap_nodes=0)

    def walk():
        try:
            kind, value = next(iterator)
        except StopIteration as error:
            raise ValueError("Incomplete tree") from error
        if kind == "RETAINED_LEAF":
            match = next((mapping[x] for x in [value, *value.split("_")] if x in mapping), None)
            if match is not None and (type(match) is not int or match <= 0):
                raise ValueError("Invalid identifier mapping")
            return set() if match is None else {match}
        if kind != "RETAINED_NODE":
            raise ValueError("Unknown tree row")
        annotations = ast.literal_eval(value)
        left, right = walk(), walk()
        if left & right:
            counts["child_overlap_nodes"] += 1
            left -= right
        if left and right:
            counts[node_class(annotations)] += 1
        return left | right

    members = walk()
    if next(iterator, None) is not None:
        raise ValueError("Trailing tree nodes")
    n = sum(counts[k] for k in ("explicit_duplication_nodes", "explicit_speciation_nodes", "default_speciation_nodes"))
    fraction = None if not n else Fraction(counts["explicit_duplication_nodes"], n)
    return dict(**counts, informative_nodes=n, mapped_members=len(members),
        duplication_fraction=None if fraction is None else str(fraction))


def bins(families):
    values = {name: None if row["duplication_fraction"] is None else Fraction(row["duplication_fraction"])
              for name, row in families.items()}
    available = [value for value in values.values() if value is not None]
    if any(not 0 <= value <= 1 for value in available):
        raise ValueError("Invalid duplication fraction")
    median = statistics.median(available) if available else None
    result = {"lower_duplication_fraction": [], "upper_duplication_fraction": [], "missing_duplication_fraction": []}
    for name, value in sorted(values.items()):
        key = "missing_duplication_fraction" if value is None else (
            "lower_duplication_fraction" if value <= median else "upper_duplication_fraction")
        result[key].append(name)
    return None if median is None else str(median), result


def validate_native(stdout, families):
    seen = set()
    fields = ("mapped_members", "explicit_duplication_nodes", "explicit_speciation_nodes",
              "default_speciation_nodes", "child_overlap_nodes")
    for line in stdout.splitlines():
        if not line.startswith("DUPLICATION_FEATURE\t"):
            continue
        parts = line.split("\t")
        if len(parts) != 7 or parts[1] in seen or parts[1] not in families:
            raise ValueError("Invalid native feature row")
        name = parts[1]
        seen.add(name)
        if [int(x) for x in parts[2:]] != [families[name][key] for key in fields]:
            raise ValueError("Native feature counts disagree")
    if seen != set(families):
        raise ValueError("Incomplete native feature comparison")


def extract(repo, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    results = repo / "benchmark_tools/results"
    admission_path = results / "swiss_retained_mapping_v2_20260923.json"
    admission = read_frozen(admission_path, ADMISSION_SHA)
    if admission["exact_families"] != 18 or not all(r["exact_match"] for r in admission["families"].values()):
        raise ValueError("Require exact mapping of all families")
    helper = record(Path(__file__).with_name("audit_swiss_retained_mapping.py"))
    if helper["sha256"] != HELPER_SHA:
        raise ValueError("Changed frozen reconstruction helper")
    protocol = record(results / "SWISS_DUPLICATION_FEATURE_PROTOCOL_20260923.md")
    if protocol["sha256"] != PROTOCOL_SHA:
        raise ValueError("Changed frozen feature protocol")
    native_script = Path(__file__).with_name("check_swiss_duplication_features.drw")
    records = [record(admission_path), helper, protocol, record(native_script), *admission["checked_inputs"]]
    for item in records:
        check(item)
    reference, image, mapping_path = [Path(r["path"]) for r in admission["checked_inputs"][:3]]
    scripts = [Path(r["path"]) for r in admission["checked_inputs"][3:5]]
    command = ["singularity", "exec", str(image), "darwin", "-E"]
    native = subprocess.run(command, input=f"reference := '{reference}':\n" +
        "\n".join(p.read_text() for p in scripts), text=True, capture_output=True, check=True, timeout=60)
    with gzip.open(mapping_path, "rt") as stream:
        mapping = json.load(stream)["mapping"]
    if compare(native.stdout, mapping) != admission["families"]:
        raise ValueError("Fresh relationship reconstruction differs")
    trees = {name: [] for name in admission["families"]}
    for line in native.stdout.splitlines():
        if line.startswith(("RETAINED_LEAF\t", "RETAINED_NODE\t")):
            kind, name, value = line.split("\t")
            trees[name].append((kind, value))
    families = {name: features(rows, mapping) for name, rows in trees.items()}
    for name, row in families.items():
        if row["mapped_members"] != admission["families"][name]["reference_members"]:
            raise ValueError("Feature membership differs from admitted mapping")
    median, strata = bins(families)
    assignments = ["ReadProgram('/benchmark/lib/darwinit'):", "MappedLeaves := table():"]
    for name, family in admission["families"].items():
        for label, number in family["mapped_labels"].items():
            if any(c in name + label for c in "'\\\n\r|"):
                raise ValueError("Unsupported native label quoting")
            assignments.append(f"MappedLeaves['{name}|{label}'] := {number}:")
    independent = subprocess.run(command, input=f"reference := '{reference}':\n" +
        "\n".join(assignments) + "\n" + native_script.read_text(), text=True,
        capture_output=True, check=True, timeout=60)
    validate_native(independent.stdout, families)
    for item in records:
        check(item)
    report = dict(status="mapped_tree_duplication_features_unscored", families=families,
        median_fraction=median, primary_strata=strata, source=record(__file__),
        checked_inputs=records, native_stderr=native.stderr,
        independent_native_counts_checked=True,
        independent_native_stdout_sha256=hashlib.sha256(independent.stdout.encode()).hexdigest(),
        independent_native_stderr=independent.stderr,
        prediction_statistics_evaluated=False, publication_ready=False,
        limitations=["Annotation fraction on mapped informative nodes, not an evolutionary duplication rate.",
            "Unannotated nodes default to S in reference reconstruction but are separately counted here.",
            "Reference-derived and development-exposed; not independent of reference pair labels.",
            "Alias overlaps follow native semantics; informative nodes need not equal unique genes minus one."])
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    extract(args.repo.resolve(), args.output)
