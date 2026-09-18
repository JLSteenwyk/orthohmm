"""Count native SwissTrees relations incident to frozen-input missing identities."""

import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_qfo_swiss_counts import HEADER, LABELS, read_raw
from benchmark_tools.inventory_swiss_annotations import COUNTS_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

ALIASES_SHA = "b69f74e35aa3fa3af99a33b3b27f830d79df2abd3935b84e911f981186ac33ec"


def incident_counts(path, families, missing, expected_truth):
    all_counts, truth, members = read_raw(path, families)
    if truth != expected_truth:
        raise ValueError("Reference relation identities or truth differ")
    affected = {family: Counter({label: 0 for label in LABELS}) for family in families}
    affected_genes, both_missing, seen = set(), 0, set()
    with gzip.open(path, "rt", newline="") as stream:
        if stream.readline().rstrip("\r\n") != HEADER:
            raise ValueError("Changed raw header")
        for family, a, b, label in csv.reader(stream, delimiter="\t"):
            key = (family, *sorted((a, b)))
            if key in seen or key not in truth or label not in LABELS or truth[key] != (label in ("TP", "FN")):
                raise ValueError("Raw relation changed during counting")
            seen.add(key)
            if a in missing or b in missing:
                affected[family][label] += 1
                affected_genes.update({a, b} & missing)
                both_missing += int(a in missing and b in missing)
    if seen != set(truth):
        raise ValueError("Incomplete second raw pass")
    totals = {label: sum(row[label] for row in affected.values()) for label in LABELS}
    return {"all_counts": all_counts, "members": members, "affected_counts": affected,
            "affected_totals": totals, "affected_relations": sum(totals.values()),
            "both_endpoints_missing": both_missing, "missing_genes_incident": sorted(affected_genes)}


def audit(count_path, alias_path):
    inputs = [record(count_path), record(alias_path)]
    if [r["sha256"] for r in inputs] != [COUNTS_SHA, ALIASES_SHA]:
        raise ValueError("Changed frozen count or identity audit")
    counts, aliases = json.loads(count_path.read_text()), json.loads(alias_path.read_text())
    missing = {g for g, row in aliases["identities"].items() if row["status"] == "missing"}
    if missing != set(aliases["summary"]["missing_genes"]) or len(missing) != 14:
        raise ValueError("Unexpected missing identity inventory")
    families = counts["families"]
    raws = [row["raw_file"] for row in counts["methods"]]
    for identity in raws:
        check(identity)
    _, truth, members = read_raw(Path(raws[0]["path"]), families)
    if len(truth) != counts["reference_relation_count"]:
        raise ValueError("Reference relation count changed")
    methods = []
    for method, raw in zip(counts["methods"], raws):
        result = incident_counts(Path(raw["path"]), families, missing, truth)
        expected = {r["family"]: r for r in method["families"]}
        if result.pop("all_counts") != {f: r["counts_without_prior"] for f, r in expected.items()}:
            raise ValueError("Native confusion counts differ from admitted counts")
        if result.pop("members") != members or members != {f: set(r["represented_genes"]) for f, r in expected.items()}:
            raise ValueError("Reference membership changed")
        if set(result["missing_genes_incident"]) != missing:
            raise ValueError("Missing identities not represented in raw relations")
        methods.append({"method": method["method"], "raw_file": raw, **result})
    for identity in [*inputs, *raws]:
        check(identity)
    return {"status": "missing_input_reference_relations_counted", "source": record(__file__),
            "reader": record(Path(__file__).with_name("audit_qfo_swiss_counts.py")), "inputs": inputs,
            "missing_genes": sorted(missing), "reference_relation_count": len(truth),
            "families": families, "methods": methods, "scores_recomputed": False,
            "limitations": ["Missing identities are established for frozen factorial inputs, not independently reconstructed historical inputs.",
                            "Incident raw labels measure retained scoring behavior; they do not establish why input records are missing.",
                            "Counts are native raw relations before the half-count plus-one prior, not independent observations.",
                            "No benchmark denominator changes, available-input sensitivity scores or confidence intervals computed."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("counts", "aliases", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.counts, args.aliases)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
