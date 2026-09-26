"""Bind retained OrthoBench comparator scores to current prediction-file readback."""

import argparse
import json
import math
from pathlib import Path

from benchmark_tools.audit_failed_recovery_refinement import record
from benchmark_tools.normalize_three_kingdoms_orthogroups import iter_fastoma, iter_orthomcl
from benchmark_tools.orthobench_stage_diagnostics import read_clusters
from benchmark_tools.score_orthobench_partition import score_partition

COMPARISON_SHA = "094f842ad8211450519d4238c0b3a5020a7969049d4b7306f5465d5acbf914a3"
REFERENCE_SHA = "660ead29c5b6ac0b8278cd1e62cdcdb0a513db81dda317e634b805e661d70ba9"


def agrees(score, expected):
    return (all(math.isclose(score[k], expected[k + "_percent"], rel_tol=0, abs_tol=1e-8)
                for k in ("f_score", "precision", "recall"))
            and score["exact_refogs"] == expected["exact_refogs"])


def audit(root):
    results = root / "benchmark_tools/results"
    comparison = results / "publication_comparison_orthomcl_complete_20260916.json"
    reference = results / "orthobench_paired_uncertainty_20260916.json"
    checked = [record(comparison), record(reference)]
    if [r["sha256"] for r in checked] != [COMPARISON_SHA, REFERENCE_SHA]:
        raise ValueError("Changed retained score/reference report")
    for name in ("audit_retained_ob_comparators.py", "normalize_three_kingdoms_orthogroups.py",
                 "orthobench_stage_diagnostics.py", "score_orthobench_partition.py"):
        checked.append(record(Path(__file__).with_name(name).resolve()))
    expected = {r["key"]: r["orthobench"] for r in json.loads(comparison.read_text())["methods"]}
    refs = json.loads(reference.read_text())["inputs"]
    checked.extend([*refs["references"], *refs["uncertain"]])
    for item in checked:
        if record(item["path"]) != item:
            raise ValueError("Changed reference bytes")
    reference_groups = {Path(r["path"]).name: set(Path(r["path"]).read_text().splitlines()) for r in refs["references"]}
    uncertain = {Path(r["path"]).name: set(Path(r["path"]).read_text().splitlines()) for r in refs["uncertain"]}
    if len(reference_groups) != 70:
        raise ValueError("Wrong reference inventory")
    base = root / "benchmarks/results"
    paths = [
        ("orthofinder_3_1_5_sequence_only", base / "orthofinder_v3_sequence_only_20260830/orthogroups.txt", None),
        ("sonicparanoid_2_0_9", base / "orthogroups_sonicparanoid.txt", None),
        ("proteinortho_6_3_6", base / "orthogroups_proteinortho.txt", None),
        ("fastoma_0_3_5", base / "fastoma_run/work/15/57adea05f3b88b17c43ee1e856acc3/OrthologousGroups.tsv", iter_fastoma),
        ("orthomcl_1_4", root.parents[1] / "SOFTWARE/ORTHOMCLV1.4/Jul_25/all_orthomcl.out", iter_orthomcl),
    ]
    rows = []
    for key, path, reader in paths:
        identity = record(path)
        groups = [set(genes) for _, genes in reader(path)] if reader else read_clusters(path)
        score = score_partition(groups, reference_groups, uncertain)
        matches = agrees(score, expected[key])
        rows.append(dict(key=key, prediction=identity, score=score, retained_score=expected[key],
                         agrees_with_retained_score=matches, parser=reader.__name__ if reader else "read_clusters"))
        checked.append(identity)
    if any(record(item["path"]) != item for item in checked):
        raise ValueError("Inputs changed during score readback")
    return dict(status="retained_ob_comparator_readback", all_scores_agree=all(r["agrees_with_retained_score"] for r in rows),
                rows=rows, checked_records=checked, publication_ready=False,
                limitations=["Current file readback matching retained statistics, not proof of historical consumption.",
                    "Does not establish native conversion, commands, versions, timing or input equivalence.",
                    "Generic orthogroups_fastoma.txt is a root-HOG diagnostic, not the final OG input.",
                    "No inference rerun or replacement of retained scores; mismatches remain explicit."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    report = audit(args.root.resolve())
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    raise SystemExit(0 if report["all_scores_agree"] else 1)
