"""Descriptive native SwissTrees scores in unchanged duplication-annotation bins."""

import argparse
import csv
from fractions import Fraction
import gzip
import json
from pathlib import Path
from statistics import median
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import export_native_qfo_swiss_sequence_strata as sequence
from benchmark_tools.audit_qfo_swiss_counts import read_raw
from benchmark_tools.export_native_factorial_progress import require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

CELLS, METRICS = sequence.CELLS, sequence.METRICS
NAMES = ("all", "lower_duplication_fraction", "upper_duplication_fraction", "missing_duplication_fraction")
PINS = {
    "native_qfo_swiss_sequence_strata_20261006_v1/report.json":
        "4650f614d2fdeebd65cbcbf4999c61c83cc549560d2738bfbb5502f686f660e3",
    "native_qfo_swiss_sequence_strata_readback_20261006.json":
        "8d0356bc0e159d204939271578e23a6a5b6efcf6198b75a18c0f3d38e0e67c0a",
    "swiss_duplication_features_v2_20260923.json":
        "97b0c4755d6a9df258d5c3f60fc0d5d25f1e5c09c42216c754a245a67d1942ec",
    "swiss_retained_mapping_v2_20260923.json":
        "0d3e736a350782609c68764bc19387500ea64bd5045c5fcfb41930dd89c1ce9d",
    "SWISS_DUPLICATION_FEATURE_PROTOCOL_20260923.md":
        "2bca0f2b2a0c3a5426be268484bab8302b404d66cfdb29c5872afefd41c79b2f",
}
IDENTIFIERS_SHA = "1c10f6ce5e53ebc3148dde02d16268b225c3c8817c952d1274daa41acbf9eb4d"


def check_memberships(memberships, mapping, identifiers):
    require(set(mapping["families"]) == set(memberships), "Changed retained family inventory")
    for family, genes in memberships.items():
        original = mapping["families"][family]
        values = [identifiers.get(g) for g in genes]
        retained = list(original["mapped_labels"].values())
        require(all(type(v) is int and v > 0 for v in [*values, *retained])
                and len(set(values)) == len(genes) and set(values) == set(retained)
                and original["exact_match"] is True and original["mapped_members"] == len(genes),
                "Native and retained mapped-entry memberships differ")


def bins_for(memberships, feature):
    require(set(memberships) == set(feature["families"]) and memberships, "Changed feature family universe")
    universe = [g for genes in memberships.values() for g in genes]
    require(all(genes and len(set(genes)) == len(genes) for genes in memberships.values())
            and len(set(universe)) == len(universe), "Empty, duplicate or shared reference genes")
    fractions = {}
    for family, row in feature["families"].items():
        keys = ("explicit_duplication_nodes", "explicit_speciation_nodes", "default_speciation_nodes",
                "informative_nodes", "child_overlap_nodes", "mapped_members")
        require(all(type(row[k]) is int and row[k] >= 0 for k in keys), "Invalid feature counts")
        n = row["informative_nodes"]
        require(sum(row[k] for k in keys[:3]) == n and row["child_overlap_nodes"] <= n
                and row["mapped_members"] == len(memberships[family]), "Feature count or mapped-member disagreement")
        value = Fraction(row["explicit_duplication_nodes"], n) if n else None
        require(row["duplication_fraction"] == (None if value is None else str(value)), "Changed feature fraction")
        fractions[family] = value
    values = [v for v in fractions.values() if v is not None]
    midpoint = median(values) if values else None
    require(feature["median_fraction"] == (None if midpoint is None else str(midpoint)), "Changed exact median")
    bins = {name: [] for name in NAMES}
    for family, value in sorted(fractions.items()):
        bins["all"].append(family)
        bins[NAMES[3] if value is None else NAMES[1 + int(value > midpoint)]].append(family)
    require({k: bins[k] for k in NAMES[1:]} == feature["primary_strata"], "Changed frozen duplication bins")
    return bins


def project(native, memberships, feature):
    bins = bins_for(memberships, feature)
    families = sorted(memberships)
    require([(r["cell"], r["family"]) for r in native["family_rows"]] ==
            [(cell, family) for cell in CELLS for family in families], "Incomplete native family count inventory")
    values = {}
    for row in native["family_rows"]:
        value = sequence.family_statistics(row["counts_without_prior"])
        require(all(row[k] == value[k] for k in METRICS), "Changed native family statistic")
        values[row["cell"], row["family"]] = value
    rows = [dict(cell=cell, stratum=name, families=len(members), family_members=members,
                 status="descriptive" if members else "empty_bin",
                 prediction_semantics="group_clique" if cell == CELLS[0] else "resolved_native_pairs",
                 **(sequence.aggregate([values[cell, f] for f in members]) if members else dict.fromkeys(METRICS)))
            for cell in CELLS for name, members in bins.items()]
    differences = [dict(stratum=name, families=len(members), family_members=members,
                        status=rows[i]["status"],
                        **{k: None if not members else rows[i + len(bins)][k] - rows[i][k] for k in METRICS})
                   for i, (name, members) in enumerate(bins.items())]
    return bins, rows, differences


def export(repo, output):
    repo, output = Path(repo).resolve(), Path(output)
    require(not output.exists() and not output.is_symlink(), "Output already exists")
    root = repo / "benchmark_tools/results"
    documents, refs, checked = [], [], []
    for name, sha in PINS.items():
        ref = record(root / name)
        require(ref["sha256"] == sha, "Changed frozen source: " + name)
        refs.append(ref)
        checked.append(ref)
        documents.append(json.loads(Path(ref["path"]).read_text()) if name.endswith(".json") else None)
    native, reader, feature, mapping, _ = documents
    require(native["schema"] == "native_qfo_swiss_sequence_strata_v1"
            and native["source"] == record(sequence.__file__)
            and all(native[k] is False for k in ("new_uncertainty", "new_accuracy_or_resource_admission",
                                                "independent_confirmation", "publication_ready")), "Changed native projection scope")
    require(reader["report"] == refs[0] and reader["raw_rows_checked"] == 21530
            and reader["families_checked"] == 18 and reader["family_rows_checked"] == 36
            and reader["source"] == record(root.parent / "readback_native_qfo_swiss_sequence_strata.py"),
            "Changed native independent raw-count binding")
    require(feature["status"] == "mapped_tree_duplication_features_unscored"
            and feature["independent_native_counts_checked"] is True
            and feature["prediction_statistics_evaluated"] is False and feature["publication_ready"] is False
            and feature["source"] == record(root.parent / "extract_swiss_duplication_features.py")
            and refs[3] in feature["checked_inputs"] and refs[4] in feature["checked_inputs"], "Changed admitted feature scope")
    checked.extend([*native["checked_inputs"], native["source"], reader["source"], feature["source"],
                    record(read_raw.__code__.co_filename),
                    record(root / "NATIVE_QFO_SWISS_DUPLICATION_STRATA_PROTOCOL_20261006.md"),
                    record(root / "NATIVE_QFO_SWISS_DUPLICATION_ALIAS_AMENDMENT_20261006.md")])
    identifiers_ref = next(r for r in feature["checked_inputs"] if r["sha256"] == IDENTIFIERS_SHA)
    checked.append(identifiers_ref)
    for ref in checked:
        check(ref)
    memberships = json.loads(Path(native["strata"]["path"]).read_text())["family_memberships"]
    require(len(memberships) == 18 and sum(map(len, memberships.values())) == 563, "Changed native gene universe")
    with gzip.open(identifiers_ref["path"], "rt") as stream:
        identifiers = json.load(stream)["mapping"]
    check_memberships(memberships, mapping, identifiers)
    del identifiers
    binding = json.loads(Path(native["binding"]["path"]).read_text())
    truth_anchor, cells = None, []
    for cell in CELLS:
        ref = binding["bound_cells"][cell]["count_audit"]
        require(ref in checked, "Missing native count audit")
        audit = json.loads(Path(ref["path"]).read_text())
        selected = [r for r in audit["cells"] if r["cell"] == cell]
        require(len(selected) == 1 and selected[0]["raw_file"] in checked, "Missing native raw binding")
        selected = selected[0]
        counts, truth, members = read_raw(Path(selected["raw_file"]["path"]), sorted(memberships))
        require(len(truth) == 10765 and (truth_anchor is None or truth_anchor == truth)
                and all(sorted(members[f]) == sorted(memberships[f]) for f in memberships)
                and all(r["counts_without_prior"] == counts[r["family"]]
                        for r in native["family_rows"] if r["cell"] == cell), "Changed native raw counts/truth/members")
        truth_anchor = truth
        if cell == CELLS[1]:
            require(selected["resources"] is None and selected["timing_eligible"] is False
                    and selected["timing_admitted"] is False, "Failed R1 timing relabeled")
        cells.append(dict(cell=cell, raw_file=selected["raw_file"], admission=selected["admission"],
                          native_job_id=selected["native_job_id"],
                          **{k: selected[k] for k in ("timing_eligible", "timing_admitted") if k in selected}))
    bins, rows, differences = project(native, memberships, feature)
    require(feature["median_fraction"] == "7/48" and [len(v) for v in bins.values()] == [18, 9, 9, 0],
            "Changed fixed bin coverage")
    result = dict(schema="native_qfo_swiss_duplication_strata_v1", source=record(__file__),
                  native=refs[0], native_readback=refs[1], features=refs[2], mapping=refs[3], original_protocol=refs[4],
                  identifiers=identifiers_ref, membership_agreement="canonical accession/reference-entry sets; retained aliases deduplicated by entry ID",
                  checked_inputs=checked, inherited_feature_inputs=[r for r in feature["checked_inputs"] if r not in checked],
                  original_tree_traversal_repeated=False, memberships=memberships, bins=bins, rows=rows,
                  differences=differences, family_rows=native["family_rows"], cells=cells,
                  raw_relations_checked=21530, new_bootstrap_draws=0, new_uncertainty=False,
                  new_accuracy_or_resource_admission=False, independent_confirmation=False, publication_ready=False,
                  limitations=["Retrospective projection of unchanged reference-derived bins; no tuning or new hypothesis test.",
                      "Macro precision/recall then harmonic F1; empty-bin values are null, not zero.",
                      "Annotation fraction is related to reference labels, not true evolutionary duplication rate or causation.",
                      "Default-S nodes are not explicit speciation; informative nodes need not equal unique genes minus one.",
                      "Development exposure and family-size/composition confounding remain; no subgroup superiority.",
                      "Original dual traversal/acquisition inherited, not repeated; direct checks are not transitive admission.",
                      "No new bootstrap draws or historical interval transfer; initial HMM search remains on in both cells.",
                      "Failed R1 timing remains ineligible. Postprocessing is shared-host, not inference or isolated speed;",
                      "contention effects are unknown and potentially tool-dependent."])
    for ref in checked:
        check(ref)
    output.mkdir(parents=True)
    fields = ("cell", "stratum", "families", "status", *METRICS, "prediction_semantics")
    with (output / "scores.tsv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({k: "NA" if v is None else v for k, v in row.items()} for row in rows)
    lines = ["# Native SwissTrees Duplication-Annotation Bins", "",
             "Descriptive percentages; R1 minus R0 changes in percentage points. No new intervals or causal claims.", "",
             "| Bin | Families | R0 F1 | R1 F1 | F1 change | Precision change | Recall change |",
             "|---|---:|---:|---:|---:|---:|---:|"]
    for i, diff in enumerate(differences):
        values = [rows[i]["F1"], rows[i + len(bins)]["F1"], *[diff[k] for k in METRICS]]
        shown = ["NA" if v is None else format(100 * v, ".3f" if j < 2 else "+.3f") for j, v in enumerate(values)]
        lines.append("| " + " | ".join([diff["stratum"], str(diff["families"]), *shown]) + " |")
    (output / "TABLE.md").write_text("\n".join([*lines, "", *["- " + s for s in result["limitations"]], ""]))
    result["outputs"] = [record(output / name) for name in ("scores.tsv", "TABLE.md")]
    for ref in checked:
        check(ref)
    with (output / "report.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = export(args.repo, args.output)
    print(json.dumps(dict(rows=len(result["rows"]), differences=len(result["differences"]), bins=result["bins"])))
