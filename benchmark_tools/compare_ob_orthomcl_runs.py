"""Compare retained OrthoMCL partitions and per-family OrthoBench counts."""

import argparse
import json
from pathlib import Path

from benchmark_tools.audit_failed_recovery_refinement import record
from benchmark_tools.normalize_three_kingdoms_orthogroups import iter_orthomcl
from benchmark_tools.score_orthobench_partition import score_partition


def partition(rows):
    labels, seen, groups = set(), set(), set()
    for label, members in rows:
        group = frozenset(members)
        if label in labels or not group or len(group) != len(members) or seen.intersection(group):
            raise ValueError("Invalid partition membership or group labels")
        labels.add(label)
        seen.update(group)
        groups.add(group)
    return groups, seen


def differences(left, right):
    a, genes_a = partition(left)
    b, genes_b = partition(right)
    missing = genes_a - genes_b
    reduced = {group - missing for group in a if group - missing}
    return dict(april_groups=len(a), july_groups=len(b), identical_groups=len(a & b),
                april_only_genes=sorted(missing), july_only_genes=sorted(genes_b - genes_a),
                april_only_groups=sorted(sorted(g) for g in a - b),
                july_only_groups=sorted(sorted(g) for g in b - a),
                equal_after_removing_april_only_genes=reduced == b)


def audit(root):
    results = root / "benchmark_tools/results"
    inventory_path = results / "ob_orthomcl_run_gene_inventory_20260926.json"
    readback_path = results / "retained_ob_comparator_readback_20260926.json"
    inventory = json.loads(inventory_path.read_text())
    readback = json.loads(readback_path.read_text())
    if record(readback_path)["sha256"] != "5dcb4e65f277eb9debf2bbacb1d4e298c678ac13fc981e75c4ba78e4851eeb55":
        raise ValueError("Changed retained reference readback")
    expected_hashes = ["6234d4fd517826fdf4fd9c25cca14b8d0bf58aef089c407eae7a57324e53e32a",
                       "1add0db6640e0451184975b5cbdf085ff404bbf24dadb66402921cd4dfd8a227"]
    if [r["sha256"] for r in inventory["inputs"]] != expected_hashes:
        raise ValueError("Changed retained OrthoMCL inputs")
    ref_records = [r for r in readback["checked_records"] if Path(r["path"]).name.startswith("RefOG")]
    checked = [record(inventory_path), record(readback_path), *inventory["inputs"], *ref_records,
               *[record(Path(__file__).with_name(n).resolve()) for n in
                 ("compare_ob_orthomcl_runs.py", "score_orthobench_partition.py", "normalize_three_kingdoms_orthogroups.py")]]
    for item in checked:
        if record(item["path"]) != item:
            raise ValueError("Changed comparison evidence")
    reference, uncertain = {}, {}
    for item in ref_records:
        path = Path(item["path"])
        target = uncertain if path.parent.name == "low_certainty_assignments" else reference
        target[path.name] = set(path.read_text().splitlines())
    if len(reference) != 70:
        raise ValueError("Wrong reference inventory")
    raw = [list(iter_orthomcl(Path(item["path"]))) for item in inventory["inputs"]]
    delta = differences(*raw)
    scores = [score_partition([set(members) for _, members in rows], reference, uncertain) for rows in raw]
    changed_genes = {g for group in delta["april_only_groups"] + delta["july_only_groups"] for g in group}
    exposure = {name: sorted(genes & changed_genes) for name, genes in reference.items() if genes & changed_genes}
    if any(record(item["path"]) != item for item in checked):
        raise ValueError("Comparison inputs changed during readback")
    return dict(status="retained_orthomcl_partitions_compared", partition_difference=delta,
                changed_group_reference_exposure=exposure, april_score=scores[0], july_score=scores[1],
                per_family_counts_equal=scores[0]["refog_records"] == scores[1]["refog_records"],
                checked_records=checked, publication_ready=False,
                limitations=["Observed retained outputs, not proof of identical historical commands or inputs.",
                    "Reference scores can agree despite different non-reference groupings.",
                    "No inference rerun, counterfactual causal claim or general robustness conclusion."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = audit(args.root.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
