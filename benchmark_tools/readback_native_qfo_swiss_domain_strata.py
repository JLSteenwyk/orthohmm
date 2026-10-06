"""Independent stdlib raw-count and annotation readback for native domain bins."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import statistics
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import readback_native_qfo_swiss_sequence_strata as raw

CELLS = ("p0_c0_r0", "p0_c0_r1")
NAMES = ("all", "median_pfam_types_below_two", "median_pfam_types_at_least_two",
         "repeated_type_fraction_below_quarter", "repeated_type_fraction_at_least_quarter")


def annotation_features(value, source):
    length = value["length"]
    raw.require(type(length) is int and length > 0 and isinstance(value.get("pfam"), dict),
                "Invalid annotation length/namespace")
    instances, types = [], sorted(value["pfam"])
    for domain in types:
        entries = value["pfam"][domain]["instance"]
        raw.require(domain.startswith("pfam_") and entries, "Invalid Pfam entry")
        for instance in entries:
            start, end = instance[:2]
            raw.require(type(start) is int and type(end) is int and 0 <= start <= end <= length,
                        "Invalid annotation coordinates")
            instances.append((start, end, domain))
    raw.require(len(instances) == len(set(instances)), "Duplicate annotation instance")
    frequency = Counter(d for _, _, d in instances)
    return dict(source_file=source, length=length, pfam_types=types, pfam_type_count=len(types),
                pfam_instance_count=len(instances), has_repeated_pfam_type=any(n > 1 for n in frequency.values()),
                ordered_pfam_instances=[dict(start=a, end=b, domain=d) for a, b, d in sorted(instances)])


def summary(genes, features):
    available = [features[g] for g in genes if g in features]
    return dict(reference_genes=len(genes), annotated_genes=len(available),
                missing_annotation_genes=sorted(set(genes) - set(features)),
                annotated_zero_pfam_genes=sum(r["pfam_instance_count"] == 0 for r in available),
                annotated_multiple_pfam_types=sum(r["pfam_type_count"] > 1 for r in available),
                annotated_repeated_pfam_type=sum(r["has_repeated_pfam_type"] for r in available),
                median_pfam_types_among_annotated=statistics.median(r["pfam_type_count"] for r in available) if available else None,
                distinct_type_sets_among_annotated=len({tuple(r["pfam_types"]) for r in available}),
                median_length_among_annotated=statistics.median(r["length"] for r in available) if available else None)


def verify_projection(report, counts, inventory):
    memberships = report["memberships"]
    families = sorted(memberships)
    genes = [g for f in families for g in memberships[f]]
    raw.require(len(families) == 18 and len(set(genes)) == len(genes)
                and set(genes) == set(inventory["genes"]) and set(families) == set(inventory["families"]),
                "Changed feature membership/coverage")
    bins = {name: [] for name in NAMES}
    bins["all"] = families
    for family in families:
        s = summary(memberships[family], inventory["genes"])
        raw.require(s == inventory["families"][family] and not s["missing_annotation_genes"],
                    "Changed annotation family summary")
        bins[NAMES[1 if s["median_pfam_types_among_annotated"] < 2 else 2]].append(family)
        bins[NAMES[3 if 4 * s["annotated_repeated_pfam_type"] < s["reference_genes"] else 4]].append(family)
    raw.require(bins == report["bins"] and [len(bins[k]) for k in NAMES] == [18, 12, 6, 15, 3]
                and bins[NAMES[2]] == ["APP", "BAR", "HOX", "NOX", "TRFE", "VATB"]
                and bins[NAMES[4]] == ["MAPT", "PSEN", "TRFE"], "Wrong frozen domain bins")
    raw.require([(r["cell"], r["stratum"]) for r in report["rows"]] ==
                [(cell, name) for cell in CELLS for name in NAMES], "Incomplete or duplicate score rows")
    raw.require([(r["cell"], r["family"]) for r in report["family_rows"]] ==
                [(cell, f) for cell in CELLS for f in families], "Incomplete or duplicate family rows")
    for row in report["family_rows"]:
        c = counts[row["cell"]]
        raw.require(row["counts_without_prior"] == dict(c[row["family"]])
                    and all(raw.equal(row[k], v) for k, v in raw.stats(c, [row["family"]]).items()),
                    "Wrong family counts/statistics")
    for row in report["rows"]:
        members = bins[row["stratum"]]
        raw.require(row["families"] == len(members) and row["family_members"] == members
                    and row["prediction_semantics"] == ("group_clique" if row["cell"] == CELLS[0] else "resolved_native_pairs")
                    and all(raw.equal(row[k], v) for k, v in raw.stats(counts[row["cell"]], members).items()),
                    "Wrong domain statistic/membership/semantics")
    raw.require([r["stratum"] for r in report["differences"]] == list(NAMES), "Incomplete domain differences")
    for diff in report["differences"]:
        members = bins[diff["stratum"]]
        left, right = [raw.stats(counts[cell], members) for cell in CELLS]
        raw.require(diff["family_members"] == members and diff["families"] == len(members)
                    and all(raw.equal(diff[k], right[k] - left[k]) for k in raw.METRICS), "Wrong domain difference")


def verify(path, digest):
    report_ref = raw.record(path)
    raw.require(report_ref["sha256"] == digest, "Changed domain report")
    report = json.loads(Path(path).read_text())
    raw.require(report["schema"] == "native_qfo_swiss_domain_strata_v1"
                and report["source"] == raw.record(Path(__file__).with_name("export_native_qfo_swiss_domain_strata.py"))
                and report["native"]["sha256"] == "4650f614d2fdeebd65cbcbf4999c61c83cc549560d2738bfbb5502f686f660e3"
                and report["native_readback"]["sha256"] == "8d0356bc0e159d204939271578e23a6a5b6efcf6198b75a18c0f3d38e0e67c0a"
                and report["annotations"]["sha256"] == "d5e269158c2fb603acc0805a1c140a88758342d7b8e32b3094750ade22fbc06c"
                and report["original_protocol"]["sha256"] == "8316866d11f988cef6a3e88f07e2802d1557a7ae581b100b85884f704b0afee8"
                and type(report["new_bootstrap_draws"]) is int and report["new_bootstrap_draws"] == 0
                and all(report[k] is False for k in ("new_uncertainty", "new_accuracy_or_resource_admission",
                                                    "independent_confirmation", "publication_ready",
                                                    "unselected_annotation_sources_rechecked")), "Wrong domain source/scope")
    checked = [report_ref, report["source"], raw.record(raw.__file__), report["native"],
               report["native_readback"], report["annotations"], report["original_protocol"],
               *report["checked_inputs"], *report["outputs"]]
    for ref in checked:
        raw.require(raw.record(ref["path"]) == ref, "Changed direct evidence")
    inventory = json.loads(Path(report["annotations"]["path"]).read_text())
    native = json.loads(Path(report["native"]["path"]).read_text())
    members = json.loads(Path(native["strata"]["path"]).read_text())["family_memberships"]
    raw.require(members == report["memberships"] and sum(map(len, members.values())) == 563
                and inventory["prediction_statistics_evaluated"] is False, "Changed original feature universe")
    selected = [r for r in inventory["annotation_sources"] if r["selected_accessions"]]
    raw.require(len(inventory["annotation_sources"]) == 78 and selected == report["selected_annotation_sources"]
                and report["unselected_annotation_source_count"] == 78 - len(selected), "Wrong selected annotation sources")
    features = {}
    for source in selected:
        ref = {k: source[k] for k in ("path", "bytes", "sha256")}
        raw.require(ref in checked, "Selected annotation not checked")
        data = json.loads(Path(ref["path"]).read_text())["feature"]
        wanted = sorted(set(data) & set(inventory["genes"]))
        raw.require(wanted == source["selected_accessions"], "Wrong annotation membership")
        for gene in wanted:
            raw.require(gene not in features, "Duplicate annotation accession")
            features[gene] = annotation_features(data[gene], Path(ref["path"]).name)
    raw.require(features == inventory["genes"] and summary(list(features), features) == inventory["summary"],
                "Independent annotation extraction differs")
    binding = json.loads(Path(native["binding"]["path"]).read_text())
    counts, truth_anchor, checked_rows = {}, None, 0
    raw.require([r["cell"] for r in report["cells"]] == list(CELLS), "Wrong native cells")
    for row in report["cells"]:
        cell = row["cell"]
        audit_ref = binding["bound_cells"][cell]["count_audit"]
        raw.require(audit_ref in checked, "Missing original count audit")
        audit = json.loads(Path(audit_ref["path"]).read_text())
        original = [r for r in audit["cells"] if r["cell"] == cell]
        raw.require(len(original) == 1 and all(row[k] == original[0][k] for k in row), "Changed native cell binding")
        raw.require(row["raw_file"] in checked, "Missing native raw file")
        counts[cell], truth = raw.raw_counts(row["raw_file"]["path"], members)
        raw.require(truth_anchor is None or truth_anchor == truth, "Changed native reference truth")
        truth_anchor, checked_rows = truth, checked_rows + len(truth)
        if cell == CELLS[1]:
            raw.require(original[0]["resources"] is None and row["timing_admitted"] is False
                        and row["timing_eligible"] is False, "Failed timing relabeled")
    raw.require(checked_rows == report["raw_relations_checked"] == 21530, "Wrong native raw count")
    verify_projection(report, counts, inventory)
    table = next(r for r in report["outputs"] if Path(r["path"]).name == "scores.tsv")
    with open(table["path"], newline="") as stream:
        parser = csv.DictReader(stream, delimiter="\t")
        raw.require(parser.fieldnames == ["cell", "stratum", "families", *raw.METRICS, "prediction_semantics"],
                    "Wrong exported table header")
        rows = list(parser)
    raw.require(len(rows) == len(report["rows"]), "Wrong table row count")
    for actual, expected in zip(rows, report["rows"]):
        raw.require(all(actual[k] == str(expected[k]) for k in ("cell", "stratum", "families", "prediction_semantics"))
                    and all(raw.equal(float(actual[k]), expected[k]) for k in raw.METRICS), "Wrong exported table values")
    for ref in checked:
        raw.require(raw.record(ref["path"]) == ref, "Evidence changed during readback")
    return dict(schema="native_qfo_swiss_domain_strata_readback_v1", source=raw.record(__file__), report=report_ref,
                checked_inputs=checked, proteins_checked=len(features), selected_annotation_sources_checked=len(selected),
                raw_rows_checked=checked_rows, score_rows_checked=len(rows), differences_checked=5,
                family_rows_checked=36, new_uncertainty=False, new_accuracy_or_resource_admission=False,
                independent_confirmation=False, publication_ready=False,
                limitations=["Independent annotation/raw parsing and alternative macro/harmonic arithmetic on original data.",
                             "Direct checks, not a new transitive admission, biological replicate or causal test.",
                             "No bootstrap draws, interval transfer, method tuning or inference timing repair."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    raw.require(not args.output.exists() and not args.output.is_symlink(), "Output already exists")
    result = verify(args.report, args.sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: result[k] for k in ("proteins_checked", "raw_rows_checked", "score_rows_checked")}))
