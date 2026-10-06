"""Project completed native SwissTrees cells into unchanged Pfam annotation bins."""

import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import inventory_swiss_annotations as annotations
from benchmark_tools import export_native_qfo_swiss_sequence_strata as sequence
from benchmark_tools.audit_qfo_swiss_counts import read_raw
from benchmark_tools.export_native_factorial_progress import require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

CELLS = sequence.CELLS
METRICS = sequence.METRICS
PINS = {
    "native_qfo_swiss_sequence_strata_20261006_v1/report.json":
        "4650f614d2fdeebd65cbcbf4999c61c83cc549560d2738bfbb5502f686f660e3",
    "native_qfo_swiss_sequence_strata_readback_20261006.json":
        "8d0356bc0e159d204939271578e23a6a5b6efcf6198b75a18c0f3d38e0e67c0a",
    "swiss_domain_annotation_inventory_20260917.json":
        "d5e269158c2fb603acc0805a1c140a88758342d7b8e32b3094750ade22fbc06c",
    "SWISS_DOMAIN_STRATA_PROTOCOL_20260917.md":
        "8316866d11f988cef6a3e88f07e2802d1557a7ae581b100b85884f704b0afee8",
}
PRIMARY = ("median_pfam_types_below_two", "median_pfam_types_at_least_two")
REPEATS = ("repeated_type_fraction_below_quarter", "repeated_type_fraction_at_least_quarter")


def bins_for(memberships, inventory):
    require(len(memberships) == 18 and set(memberships) == set(inventory["families"]),
            "Changed reference family inventory")
    universe = [g for genes in memberships.values() for g in genes]
    require(len(universe) == len(set(universe)) and set(universe) == set(inventory["genes"]),
            "Missing, shared or extra annotated accession")
    bins = {"all": sorted(memberships), **{name: [] for name in (*PRIMARY, *REPEATS)}}
    for family in sorted(memberships):
        summary = annotations.summarize(memberships[family], inventory["genes"])
        require(summary == inventory["families"][family] and not summary["missing_annotation_genes"],
                "Changed annotation family summary")
        bins[PRIMARY[int(summary["median_pfam_types_among_annotated"] >= 2)]].append(family)
        fraction = summary["annotated_repeated_pfam_type"] / summary["reference_genes"]
        bins[REPEATS[int(fraction >= .25)]].append(family)
    require([len(bins[k]) for k in (*PRIMARY, *REPEATS)] == [12, 6, 15, 3]
            and set(bins[PRIMARY[1]]) == {"APP", "BAR", "HOX", "NOX", "TRFE", "VATB"}
            and set(bins[REPEATS[1]]) == {"MAPT", "PSEN", "TRFE"}, "Changed fixed domain bins")
    return bins


def project(native, memberships, inventory):
    bins = bins_for(memberships, inventory)
    families = sorted(memberships)
    family_rows = native["family_rows"]
    require([(r["cell"], r["family"]) for r in family_rows] ==
            [(cell, f) for cell in CELLS for f in families], "Wrong native family count inventory")
    values = {}
    for row in family_rows:
        value = sequence.family_statistics(row["counts_without_prior"])
        require(all(row[k] == value[k] for k in METRICS), "Native family statistic differs")
        values[row["cell"], row["family"]] = value
    rows = [dict(cell=cell, stratum=name, family_members=members, families=len(members),
                 prediction_semantics="group_clique" if cell == CELLS[0] else "resolved_native_pairs",
                 **sequence.aggregate([values[cell, f] for f in members]))
            for cell in CELLS for name, members in bins.items()]
    differences = [dict(stratum=name, family_members=members, families=len(members),
                        **{k: rows[len(bins) + i][k] - rows[i][k] for k in METRICS})
                   for i, (name, members) in enumerate(bins.items())]
    return bins, rows, differences


def reextract(inventory, checked):
    found, selected_sources = {}, []
    require(len(inventory["annotation_sources"]) == 78, "Wrong original annotation source count")
    for source in inventory["annotation_sources"]:
        if not source["selected_accessions"]:
            continue
        ref = {k: source[k] for k in ("path", "bytes", "sha256")}
        check(ref)
        data = json.loads(Path(ref["path"]).read_text())["feature"]
        selected = sorted(set(data) & set(inventory["genes"]))
        require(selected == source["selected_accessions"], "Changed annotation source membership")
        for gene in selected:
            require(gene not in found, "Ambiguous annotation accession")
            found[gene] = dict(source_file=Path(ref["path"]).name,
                               **annotations.features(data[gene]))
        check(ref)
        checked.append(ref)
        selected_sources.append(source)
    require(found == inventory["genes"], "Original annotation features differ")
    return selected_sources


def export(repo, output):
    repo, output = Path(repo).resolve(), Path(output)
    require(not output.exists() and not output.is_symlink(), "Output already exists")
    root = repo / "benchmark_tools/results"
    checked, documents, refs = [], {}, {}
    for name, sha in PINS.items():
        ref = record(root / name)
        require(ref["sha256"] == sha, "Changed frozen source: " + name)
        refs[name], documents[name] = ref, json.loads(Path(ref["path"]).read_text()) if name.endswith(".json") else None
        checked.append(ref)
    native, reader, inventory = [documents[name] for name in list(PINS)[:3]]
    require(native["schema"] == "native_qfo_swiss_sequence_strata_v1"
            and native["source"] == record(sequence.__file__)
            and all(native[k] is False for k in ("new_uncertainty", "new_accuracy_or_resource_admission",
                                                "independent_confirmation", "publication_ready")),
            "Changed native projection source/scope")
    require(reader["report"] == refs[list(PINS)[0]] and reader["raw_rows_checked"] == 21530
            and reader["families_checked"] == 18 and reader["family_rows_checked"] == 36
            and reader["source"] == record(root.parent / "readback_native_qfo_swiss_sequence_strata.py"),
            "Changed independent raw-count readback binding")
    require(inventory["status"] == "prediction_independent_swiss_annotation_inventory"
            and inventory["prediction_statistics_evaluated"] is False
            and inventory["source"] == record(annotations.__file__), "Changed annotation source/scope")
    checked.extend(native["checked_inputs"])
    checked.extend([native["source"], reader["source"], inventory["source"], record(read_raw.__code__.co_filename),
                    record(root / "NATIVE_QFO_SWISS_DOMAIN_STRATA_PROTOCOL_20261006.md")])
    for ref in checked:
        check(ref)
    feature = json.loads(Path(native["strata"]["path"]).read_text())
    memberships = feature["family_memberships"]
    require(sum(map(len, memberships.values())) == 563, "Changed native protein universe")
    selected_sources = reextract(inventory, checked)
    require(annotations.summarize(list(inventory["genes"]), inventory["genes"]) == inventory["summary"],
            "Original whole-inventory summary differs")
    binding = json.loads(Path(native["binding"]["path"]).read_text())
    cell_sources, truth_anchor = [], None
    for cell in CELLS:
        ref = binding["bound_cells"][cell]["count_audit"]
        require(ref in checked, "Missing native count audit")
        audit = json.loads(Path(ref["path"]).read_text())
        selected = [r for r in audit["cells"] if r["cell"] == cell]
        require(len(selected) == 1 and selected[0]["raw_file"] in checked, "Missing native raw binding")
        selected = selected[0]
        counts, truth, members = read_raw(Path(selected["raw_file"]["path"]), sorted(memberships))
        require(len(truth) == 10765 and (truth_anchor is None or truth == truth_anchor)
                and all(sorted(members[f]) == sorted(memberships[f]) for f in memberships)
                and all(r["counts_without_prior"] == counts[r["family"]]
                        for r in native["family_rows"] if r["cell"] == cell), "Native raw counts/truth differ")
        truth_anchor = truth
        if cell == CELLS[1]:
            require(selected["resources"] is None and selected["timing_eligible"] is False
                    and selected["timing_admitted"] is False, "Failed R1 timing relabeled")
        cell_sources.append(dict(cell=cell, raw_file=selected["raw_file"], admission=selected["admission"],
                                 native_job_id=selected["native_job_id"],
                                 **{k: selected[k] for k in ("timing_eligible", "timing_admitted") if k in selected}))
    bins, rows, differences = project(native, memberships, inventory)
    result = dict(schema="native_qfo_swiss_domain_strata_v1", source=record(__file__),
                  native=refs[list(PINS)[0]], native_readback=refs[list(PINS)[1]],
                  annotations=refs[list(PINS)[2]], original_protocol=refs[list(PINS)[3]],
                  checked_inputs=checked, selected_annotation_sources=selected_sources,
                  unselected_annotation_source_count=78 - len(selected_sources),
                  unselected_annotation_sources_rechecked=False, memberships=memberships,
                  bins=bins, rows=rows, differences=differences, family_rows=native["family_rows"],
                  cells=cell_sources, raw_relations_checked=21530, new_bootstrap_draws=0,
                  new_uncertainty=False, new_accuracy_or_resource_admission=False,
                  independent_confirmation=False, publication_ready=False,
                  limitations=["Retrospective fixed-bin projections, not a new inferential test or method tuning.",
                      "Macro precision/recall then harmonic F1, not pair-pooled or mean-family F1.",
                      "Pfam counts are annotations, not complete validated domain architectures or causal mechanisms.",
                      "Repeated-type high bin has three families; bins overlap; no new or historical intervals attached.",
                      "FAS reference annotations do not independently validate FAS; no fragment or domain-loss truth.",
                      "Direct source/raw/selected-annotation checks, not full transitive admission.",
                      "R1 failed timing remains ineligible. Shared-host postprocessing time is not inference cost;",
                      "contention effects are unknown and potentially tool-dependent, not isolated speed evidence."])
    for ref in checked:
        check(ref)
    output.mkdir(parents=True)
    fields = ("cell", "stratum", "families", *METRICS, "prediction_semantics")
    with (output / "scores.tsv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    lines = ["# Native SwissTrees Domain Bins", "", "Descriptive percentages; differences R1 minus R0 in percentage points.",
             "No new intervals, significance or causal domain effects. Initial HMM search on in both cells.", "",
             "| Bin | Families | R0 F1 | R1 F1 | F1 change | Precision change | Recall change |",
             "|---|---:|---:|---:|---:|---:|---:|"]
    for i, diff in enumerate(differences):
        lines.append(f"| {diff['stratum']} | {diff['families']} | {100 * rows[i]['F1']:.3f} | "
                     f"{100 * rows[i + len(bins)]['F1']:.3f} | " +
                     " | ".join(f"{100 * diff[k]:+.3f}" for k in METRICS) + " |")
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
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = export(args.repo, args.output)
    print(json.dumps(dict(rows=len(result["rows"]), bins=len(result["bins"]),
                          selected_annotation_sources=len(result["selected_annotation_sources"]))))
