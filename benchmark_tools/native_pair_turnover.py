"""Descriptive native-pair turnover, separate from benchmark OG accuracy."""

from collections import Counter
import csv
from pathlib import Path
import time

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def _pairs(path, owners):
    result = {}
    with Path(path).open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != ["gene_a", "species_a", "gene_b", "species_b", "confidence"]:
            raise ValueError("Unexpected native confidence header")
        for row in reader:
            a, b = row["gene_a"], row["gene_b"]
            sa, sb, confidence = row["species_a"], row["species_b"], row["confidence"]
            if (not a or not b or a == b or not sa or not sb or sa == sb
                    or confidence not in ("high", "medium", "low") or None in row):
                raise ValueError("Invalid native confidence row")
            for gene, species in ((a, sa), (b, sb)):
                if gene in owners and owners[gene] != species:
                    raise ValueError("Native species ownership changed")
                owners[gene] = species
            pair = tuple(sorted((a, b)))
            if pair in result:
                raise ValueError("Duplicate native ortholog pair")
            result[pair] = confidence
    return result


def native_pair_turnover(baseline, candidate):
    started = time.monotonic()
    pins = [record(Path(path)) for path in (baseline, candidate)]
    owners = {}
    old, new = (_pairs(path, owners) for path in (baseline, candidate))
    old_keys, new_keys = set(old), set(new)
    retained = old_keys & new_keys
    old_high = {pair for pair, confidence in old.items() if confidence == "high"}
    result = {"baseline_pairs": len(old), "candidate_pairs": len(new), "retained_pairs": len(retained),
        "added_pairs": len(new_keys - old_keys), "lost_pairs": len(old_keys - new_keys),
        "baseline_confidence_counts": dict(sorted(Counter(old.values()).items())),
        "candidate_confidence_counts": dict(sorted(Counter(new.values()).items())),
        "retained_confidence_transitions": dict(sorted(Counter(old[pair] + "->" + new[pair] for pair in retained).items())),
        "baseline_high_pairs": len(old_high), "baseline_high_pairs_absent": len(old_high - new_keys),
        "baseline_high_pairs_downgraded": sum(new[pair] != "high" for pair in old_high & new_keys),
        "files": pins, "wall_s": time.monotonic() - started,
        "scope": "Descriptive native cross-species pair/confidence turnover, not reference orthology accuracy or OG F1. Shared-host readback cost, not inference runtime."}
    for pin in pins:
        check(pin)
    return result
