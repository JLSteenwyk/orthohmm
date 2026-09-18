"""Inventory retained SwissTrees domain annotations without evaluating predictions."""

import argparse
from collections import Counter
import json
from pathlib import Path
import statistics
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.snapshot_orthohmm_input_order import record

COUNTS_SHA = "2995868b0407ceda2e7422db6c3fc3b99523716e4bec8b1b43e15f296e40769b"


def features(annotation):
    length = annotation["length"]
    if type(length) is not int or length <= 0:
        raise ValueError("Invalid annotated sequence length")
    domains = annotation.get("pfam")
    if not isinstance(domains, dict):
        raise ValueError("Missing Pfam annotation namespace")
    instances = []
    for domain, value in sorted(domains.items()):
        if not domain.startswith("pfam_") or not value["instance"]:
            raise ValueError("Invalid Pfam annotation entry")
        for entry in value["instance"]:
            start, end = entry[:2]
            if type(start) is not int or type(end) is not int or not 0 <= start <= end <= length:
                raise ValueError("Invalid domain coordinates")
            instances.append((start, end, domain))
    if len(set(instances)) != len(instances):
        raise ValueError("Duplicate domain instance")
    return {"length": length, "pfam_types": sorted(domains), "pfam_type_count": len(domains),
            "pfam_instance_count": len(instances),
            "ordered_pfam_instances": [dict(start=a, end=b, domain=d) for a, b, d in sorted(instances)],
            "has_repeated_pfam_type": len(instances) > len(domains)}


def summarize(genes, annotated):
    available = [annotated[g] for g in genes if g in annotated]
    return {"reference_genes": len(genes), "annotated_genes": len(available),
            "missing_annotation_genes": sorted(set(genes) - set(annotated)),
            "annotated_zero_pfam_genes": sum(r["pfam_instance_count"] == 0 for r in available),
            "annotated_multiple_pfam_types": sum(r["pfam_type_count"] > 1 for r in available),
            "annotated_repeated_pfam_type": sum(r["has_repeated_pfam_type"] for r in available),
            "median_pfam_types_among_annotated": statistics.median(r["pfam_type_count"] for r in available) if available else None,
            "distinct_type_sets_among_annotated": len({tuple(r["pfam_types"]) for r in available}),
            "median_length_among_annotated": statistics.median(r["length"] for r in available) if available else None}


def inventory(count_path, directory):
    count_identity = record(count_path)
    if count_identity["sha256"] != COUNTS_SHA:
        raise ValueError("Changed reference-membership source")
    counts = json.loads(count_path.read_text())
    families = {r["family"]: r["represented_genes"] for r in counts["methods"][0]["families"]}
    genes = [g for members in families.values() for g in members]
    if len(families) != 18 or len(set(genes)) != len(genes):
        raise ValueError("Unexpected family inventory or shared reference gene")
    wanted = set(genes)
    found, sources, namespace_counts = {}, [], Counter()
    files = sorted(directory.glob("*.json"))
    if not files:
        raise ValueError("No annotation sources")
    for path in files:
        identity = record(path)
        data = json.loads(path.read_text())
        if not isinstance(data["feature"], dict):
            raise ValueError("Invalid feature map")
        selected = sorted(wanted.intersection(data["feature"]))
        for gene in selected:
            if gene in found:
                raise ValueError("Ambiguous accession across annotation files: " + gene)
            annotation = data["feature"][gene]
            namespace_counts.update(annotation.keys())
            found[gene] = {"source_file": path.name, **features(annotation)}
        if record(path) != identity:
            raise ValueError("Annotation changed during inventory")
        sources.append({**identity, "selected_accessions": selected})
    if record(count_path) != count_identity:
        raise ValueError("Reference membership changed")
    return {"status": "prediction_independent_swiss_annotation_inventory", "source": record(__file__),
            "reference_membership_source": count_identity, "annotation_sources": sources,
            "feature_namespaces_among_selected_records": dict(namespace_counts),
            "summary": summarize(genes, found),
            "families": {f: summarize(members, found) for f, members in families.items()},
            "genes": found, "prediction_statistics_evaluated": False,
            "limitations": ["Exact accession matches only; no fuzzy or silent alias mapping.",
                            "Missing annotation is distinct from a present record with no Pfam hits.",
                            "Pfam types, instances and repeats are annotation features, not validated domain loss or fragments.",
                            "Coordinates are retained as supplied; no coordinate-width or coverage-fraction calculation is made.",
                            "Annotations are the retained QfO FAS reference resource, not independent evidence for the FAS endpoint.",
                            "No error strata, uncertainty intervals or accuracy contrasts evaluated here."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("counts", "annotations", "output"):
        parser.add_argument("--" + flag, required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = inventory(args.counts, args.annotations)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
