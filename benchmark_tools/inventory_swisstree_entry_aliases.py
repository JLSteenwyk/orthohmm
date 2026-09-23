"""Inventory historical entry-name leaf aliases without admitting tree identities."""

import argparse
import json
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen

SOURCES = {
    "swiss_historical_fragment_admission_22117.json": "a480f32666bb96de3395213c7cad024c170318ca110c389f15e4756f8052adf8",
    "swisstree_model_feasibility_20260923.json": "6687e14f2398a602a990761378dc7623393c0654062e4f796de51f95ff60c654",
    "swisstree_multiline_model_diagnostic_20260923.json": "7e3461f619d1953da3e25f43c52ab3a5bb7af104232af1bf0288bdb73d95d2b7",
}


def aliases(leaves, genes, annotations):
    if not leaves or len(set(leaves)) != len(leaves):
        raise ValueError("Empty or duplicate leaf labels")
    rows = {}
    for gene in sorted(genes):
        annotation = annotations[gene]
        if annotation is None or annotation["accession"] != gene:
            raise ValueError("Require accession-matched historical annotation")
        entry = annotation["selection"]["name"]
        matches = {}
        for leaf in leaves:
            reasons = []
            if leaf == gene or leaf.rsplit("_", 1)[-1] == gene:
                reasons.append("accession_exact_or_suffix")
            if leaf == entry:
                reasons.append("historical_entry_name_exact")
            if reasons:
                matches[leaf] = reasons
        rows[gene] = dict(candidates=matches, historical_entry_name=entry,
            selection_class=annotation["selection_class"], sequence_sha256=annotation["sequence_sha256"],
            taxid=annotation["taxid"], status="unique_candidate" if len(matches) == 1 else "missing" if not matches else "ambiguous")
    reverse = {}
    for gene, row in rows.items():
        for leaf in row["candidates"]:
            reverse.setdefault(leaf, []).append(gene)
    collisions = {leaf: genes for leaf, genes in reverse.items() if len(genes) > 1}
    for genes in collisions.values():
        for gene in genes:
            rows[gene]["status"] = "shared_leaf_candidate"
    return dict(genes=rows, shared_leaf_candidates=collisions,
        unique_candidates=sum(r["status"] == "unique_candidate" for r in rows.values()),
        missing=sum(r["status"] == "missing" for r in rows.values()), mapping_admitted=False)


def inventory(root, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    data, records = [], [record(__file__)]
    for name, digest in SOURCES.items():
        path = root / "benchmark_tools/results" / name
        data.append(read_frozen(path, digest))
        records.append(record(path))
    annotations, models, multiline = data
    records.extend(annotations["records"] + models["inputs"] + multiline["records"])
    for item in records:
        check(item)
    families = {}
    for identifier, row in models["families"].items():
        family = row["qfo_family"]
        if identifier in multiline["families"]:
            leaf_sets = [c["leaves"] for c in multiline["families"][identifier]["candidates"]]
        else:
            leaf_sets = [row["leaf_names"]]
        families[family] = dict(identifier=identifier, chosen_candidate=None,
            candidates=[aliases(leaves, annotations["families"][family], annotations["annotations"])
                        for leaves in leaf_sets])
    if set(families) != set(annotations["families"]):
        raise ValueError("Incomplete reference family universe")
    for item in records:
        check(item)
    report = dict(status="historical_entry_alias_candidates_not_admitted", families=families, records=records,
        prediction_statistics_evaluated=False, mapping_admitted=False, publication_ready=False,
        limitations=["Entry names are from sequence-verified historical records, but current tree leaf sequence versions remain unverified.",
            "Only exact entry names are added; no species-name normalization, fuzzy matching or prediction-based choice.",
            "Multiple source candidates and ambiguous/shared aliases remain explicit; none is selected.",
            "Neither biological rooting nor ancestral event histories are validated by alias feasibility."])
    with output.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    inventory(args.root.resolve(), args.output.absolute())
