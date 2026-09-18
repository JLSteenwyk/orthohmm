"""Independently reconstruct BUSCO pair counts using partition intersections."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path


def record(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def partition(path):
    membership = {}
    sizes = []
    with Path(path).open() as stream:
        for number, line in enumerate(stream, 1):
            fields = line.split()
            if not fields:
                continue
            if fields[0].endswith(":"):
                fields = fields[1:]
            if not fields:
                raise ValueError(f"Empty labeled group: {path}:{number}")
            group = len(sizes)
            for gene in fields:
                if gene in membership:
                    raise ValueError(f"Duplicate membership: {gene} in {path}:{number}")
                membership[gene] = group
            sizes.append(len(fields))
    return membership, sizes


def choose_two(n):
    return n * (n - 1) // 2


def count(reference, predictions):
    refs, ref_sizes = partition(reference)
    preds, pred_sizes = partition(predictions)
    intersections = Counter((r, preds[g]) for g, r in refs.items() if g in preds)
    restricted = Counter(preds[g] for g in refs if g in preds)
    tp = sum(choose_two(n) for n in intersections.values())
    truth = sum(choose_two(n) for n in ref_sizes)
    positive = sum(choose_two(n) for n in restricted.values())
    fp, fn = positive - tp, truth - tp
    precision = tp / positive if positive else 0.0
    recall = tp / truth if truth else 0.0
    f1 = 2 * tp / (positive + truth) if positive + truth else 0.0
    return {"reference_orthogroups": len(ref_sizes), "reference_genes": len(refs),
            "reference_genes_in_prediction": sum(restricted.values()),
            "predicted_orthogroups": len(pred_sizes), "true_positive_gene_pairs": tp,
            "false_positive_gene_pairs": fp, "false_negative_gene_pairs": fn,
            "precision": precision, "recall": recall, "f_score": f1,
            "reference_gene_coverage": sum(restricted.values()) / len(refs) if refs else 0.0}


def compare(actual, expected):
    for key, value in actual.items():
        reported = expected[key]
        if isinstance(value, int):
            matches = type(reported) is int and value == reported
        else:
            matches = isinstance(reported, (float, int)) and abs(value - reported) <= 1e-12
        if not matches:
            raise ValueError(f"Independent count mismatch for {key}: {value} != {reported}")


def audit_panel(repo, panel_path, reference):
    evidence = [record(panel_path), record(reference), record(__file__)]
    panel = json.loads(Path(panel_path).read_text())
    if evidence[1]["sha256"] != panel["dataset"]["reference_sha256"]:
        raise ValueError("Reference differs from historical panel")
    methods = []
    keys = set()
    for method in panel["methods"]:
        if method["key"] in keys:
            raise ValueError("Duplicate method key")
        keys.add(method["key"])
        directory = Path(repo) / method["result_dir"]
        groups = record(directory / "orthogroups.txt")
        score = record(directory / "score.txt")
        for item, key in ((groups, "orthogroups_sha256"), (score, "score_sha256")):
            if item["sha256"] != method["provenance"][key]:
                raise ValueError(f"Historical source changed: {item['path']}")
        evidence.extend([groups, score])
        counts = count(reference, groups["path"])
        compare(counts, method["score"])
        methods.append({"key": method["key"], "counts": counts,
                        "groups": groups, "historical_score_file": score})
    for item in evidence:
        if record(item["path"]) != item:
            raise ValueError(f"Source changed during audit: {item['path']}")
    return {"status": "historical_three_kingdoms_pair_counts_independently_verified",
            "method_count": len(methods), "methods": methods, "evidence": evidence,
            "accuracy_admitted": False,
            "limitations": ["Verifies normalized groups against historical summary; not native conversion.",
                            "Does not authenticate historical runtime or matched input consumption.",
                            "Only reference-gene co-membership is scored, including within-species pairs.",
                            "No penalty for predictions involving genes outside the reference universe."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--panel", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit_panel(args.repo, args.panel, args.reference)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"status": result["status"], "methods": result["method_count"]}))
