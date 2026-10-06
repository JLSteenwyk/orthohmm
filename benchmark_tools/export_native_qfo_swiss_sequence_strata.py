"""Descriptive native SwissTrees scores in unchanged input-only sequence bins."""

import argparse
import csv
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import audit_native_qfo_swiss_counts as ordinary
from benchmark_tools import audit_recovered_native_qfo_swiss_counts as recovered
from benchmark_tools import bind_native_qfo_swiss_uncertainty as binder
from benchmark_tools import export_native_qfo_scientific_scores as reporter
from benchmark_tools import prepare_corrected_swiss_sequence_strata as features
from benchmark_tools.audit_qfo_swiss_counts import read_raw
from benchmark_tools.export_native_factorial_progress import load, require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

CELLS = ("p0_c0_r0", "p0_c0_r1")
METRICS = ("F1", "PPV", "TPR")
STRATA_SHA = "912363276fa5456f11aeff9f398fce71aec8d377c8c8eda7b144105945e554f1"


def family_statistics(counts):
    require(set(counts) == {"TP", "FP", "FN", "TN"}
            and all(type(v) is int and v >= 0 for v in counts.values()), "Invalid native family counts")
    tp, fp, fn = (counts[k] / 2 + 1 for k in ("TP", "FP", "FN"))
    precision, recall = tp / (tp + fp), tp / (tp + fn)
    return dict(F1=2 * precision * recall / (precision + recall), PPV=precision, TPR=recall)


def aggregate(values):
    require(values, "Empty aggregate")
    precision = sum(row["PPV"] for row in values) / len(values)
    recall = sum(row["TPR"] for row in values) / len(values)
    return dict(F1=2 * recall * precision / (recall + precision), PPV=precision, TPR=recall)


def project(cells, strata):
    require([row["cell"] for row in cells] == list(CELLS), "Wrong native cell inventory")
    memberships = strata["family_memberships"]
    require(len(memberships) == 18, "Wrong reference family count")
    recomputed = features.define_strata(memberships, strata["genes"])
    require(all(recomputed[k] == strata[k] for k in recomputed), "Changed frozen feature bins/descriptors")
    families = sorted(memberships)
    bins = {"all": families, **strata["primary_strata"], **strata["secondary_strata"]}
    require(len(bins) == 11 and all(len(set(v)) == len(v) and set(v) <= set(families)
                                  for v in bins.values()), "Wrong or invalid bin inventory")
    rows, family_rows = [], []
    for cell in cells:
        require([r["family"] for r in cell["families"]] == families, "Changed native family order/coverage")
        values = {}
        for row in cell["families"]:
            family = row["family"]
            require(row["represented_genes"] == sorted(memberships[family]), "Feature/native member universe differs")
            value = family_statistics(row["counts_without_prior"])
            require(value == row["statistics_with_prior"], "Changed native family statistic")
            values[family] = value
            family_rows.append(dict(cell=cell["cell"], family=family, **value,
                                    counts_without_prior=row["counts_without_prior"]))
        require(aggregate(list(values.values())) == cell["aggregate"], "Native aggregate does not reproduce")
        require(math.isclose(cell["aggregate"]["F1"], cell["native_endpoint_f1"], rel_tol=0, abs_tol=5e-8),
                "Native rounded endpoint differs")
        for name, members in bins.items():
            rows.append(dict(cell=cell["cell"], stratum=name, families=len(members), family_members=members,
                status="descriptive" if members else "empty_bin",
                prediction_semantics="group_clique" if cell["cell"] == CELLS[0] else "resolved_native_pairs",
                **(aggregate([values[f] for f in members]) if members else dict.fromkeys(METRICS))))
    left = {r["stratum"]: r for r in rows if r["cell"] == CELLS[0]}
    right = {r["stratum"]: r for r in rows if r["cell"] == CELLS[1]}
    differences = [dict(stratum=name, families=left[name]["families"], family_members=left[name]["family_members"],
        status=left[name]["status"], **{metric: None if left[name][metric] is None else
        right[name][metric] - left[name][metric] for metric in METRICS}) for name in bins]
    return rows, differences, family_rows


def export(binding_path, binding_sha, strata_path, output):
    require(not output.exists() and not output.is_symlink(), "Output already exists")
    evidence = [record(Path(__file__).parent / "results/NATIVE_QFO_SWISS_SEQUENCE_STRATA_PROTOCOL_20261006.md"),
                record(read_raw.__code__.co_filename)]
    binding, binding_ref = load(binding_path, binding_sha, evidence)
    require(binding["schema"] == "native_qfo_retained_swiss_uncertainty_binding_v1"
            and binding["source"] == record(binder.__file__) and set(binding["bound_cells"]) == set(CELLS)
            and all(binding[k] is False for k in ("new_accuracy_or_resource_admission",
                    "independent_confirmation", "publication_ready")), "Wrong native binding source/scope")
    strata, strata_ref = load(strata_path, STRATA_SHA, evidence)
    require(strata["status"] == "corrected_swiss_sequence_strata_prepared_unscored"
            and strata["prediction_statistics_evaluated"] is False and strata["publication_ready"] is False
            and strata["source"] == record(features.__file__)
            and strata["helper"]["sha256"] == features.HELPER_SHA
            and strata["protocol"]["sha256"] == features.PROTOCOL_SHA, "Wrong frozen feature source/scope")
    evidence.extend([binding["source"], strata["source"], strata["helper"], strata["protocol"], *strata["inputs"]])
    snapshot, snapshot_ref = load(binding["snapshot"]["path"], binding["snapshot"]["sha256"], evidence)
    require(snapshot_ref == binding["snapshot"] and snapshot["source"] == record(reporter.__file__)
            and snapshot["schema"] == "native_qfo_scientific_reporting_snapshot_v1", "Wrong scientific snapshot")
    evidence.append(snapshot["source"])
    plan, plan_ref = load(snapshot["plan"]["path"], snapshot["plan"]["sha256"], evidence)
    require(plan_ref == snapshot["plan"], "Changed native plan binding")
    cells, truth_anchor, member_anchor = [], None, None
    snapshot_rows = {r["cell"]: r for r in snapshot["rows"]}
    for index, name in enumerate(CELLS):
        bound = binding["bound_cells"][name]
        audit, audit_ref = load(bound["count_audit"]["path"], bound["count_audit"]["sha256"], evidence)
        module = ordinary if index == 0 else recovered
        status = ("supplied_native_swiss_family_counts_verified" if index == 0 else
                  "supplied_recovered_native_swiss_family_counts_verified")
        require(audit_ref == bound["count_audit"] and audit["source"] == record(module.__file__)
                and audit["schema"] == ("native_qfo_swiss_family_count_audit_v1" if index == 0 else
                                        "recovered_native_qfo_swiss_family_count_audit_v1")
                and audit["status"] == status and audit["families"] == binding["families"]
                and audit["retained_counts"] == binding["retained_counts"]
                and all(audit[k] is False for k in ("historical_intervals_attached",
                    "new_accuracy_or_resource_admission", "independent_confirmation", "publication_ready")),
                "Changed count audit source/scope")
        selected = [r for r in audit["cells"] if r["cell"] == name]
        require(len(selected) == 1, "Missing or duplicate native audited cell")
        cell = selected[0]
        require(all(cell[k] == bound[k] for k in ("admission", "index", "native_job_id", "native_endpoint_f1"))
                and cell["aggregate"] == bound["count_aggregate"]
                and cell["admission"] == snapshot_rows[name]["admission"]
                and cell["native_endpoint_f1"] == snapshot_rows[name]["scores"]["SwissTrees"]
                and snapshot_rows[name]["status"] == ("supplied_native_admission" if index == 0 else
                                                     "supplied_recovered_scientific_admission")
                and cell["index"] == index + 6 and cell["raw_file"] in audit["checked_inputs"],
                "Native count/snapshot identity differs")
        admission, admission_ref = load(cell["admission"]["path"], cell["admission"]["sha256"], evidence)
        require(admission_ref == cell["admission"] and admission["accuracy_admitted"] is True
                and admission["conversion"]["cell"] == name and admission["conversion"]["plan"] == plan_ref
                and admission["conversion"]["input_fastas"] == strata["fasta_inputs"]
                and plan["runs"][index + 6]["inputs"] == strata["fasta_inputs"],
                "Feature/native input identities or admission differ")
        if index:
            require(admission["resources"] is None and admission["scientific_timings_admitted"] is False
                    and admission["eligible_for_timing_comparison"] is False and cell["resources"] is None
                    and cell["timing_admitted"] is False and cell["timing_eligible"] is False,
                    "Recovered cell relabels failed timing")
        check(cell["raw_file"])
        evidence.extend([audit["source"], cell["raw_file"]])
        counts, truth, members = read_raw(Path(cell["raw_file"]["path"]), audit["families"])
        require(len(truth) == audit["reference_relation_count"] and (truth_anchor is None or truth == truth_anchor)
                and (member_anchor is None or members == member_anchor)
                and all(dict(counts[r["family"]]) == r["counts_without_prior"]
                        and sorted(members[r["family"]]) == r["represented_genes"] for r in cell["families"]),
                "Actual raw family counts/truth/members differ")
        truth_anchor, member_anchor = truth, members
        cells.append(cell)
    rows, differences, family_rows = project(cells, strata)
    for ref in evidence:
        check(ref)
    output.mkdir(parents=True)
    fields = ("cell", "stratum", "families", "status", *METRICS, "prediction_semantics")
    with (output / "scores.tsv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({k: "NA" if v is None else v for k, v in row.items()} for row in rows)
    lines = ["# Native SwissTrees Frozen Sequence Strata", "",
        "Descriptive percentages; F1 is the harmonic mean of macro precision and recall. Differences are R1 minus R0 in percentage points.",
        "Empty bins are NA, not zero. These are development-exposed overlapping bins, with no new intervals or significance claims.", ""]
    for diff in differences:
        name = diff["stratum"]
        lines.extend([f"## {name} ({diff['families']} families)", "",
            "| Native cell | F1 | Precision | Recall |", "|---|---:|---:|---:|"])
        for row in (r for r in rows if r["stratum"] == name):
            values = ["NA" if row[k] is None else f"{100 * row[k]:.3f}" for k in METRICS]
            lines.append("| " + " | ".join([row["cell"], *values]) + " |")
        delta = ["NA" if diff[k] is None else f"{100 * diff[k]:+.3f}" for k in METRICS]
        lines.extend(["| R1 minus R0 | " + " | ".join(delta) + " |", ""])
    lines.append("Initial HMM search remains on in both cells. R0 uses group-clique pairs; R1 uses resolved native pairs. Global entropy is not divergence, relative shortness is not fragmentation, and absent fragment text does not prove completeness. Failed R1 timing stays ineligible.")
    with (output / "TABLE.md").open("x") as stream:
        stream.write("\n".join(lines) + "\n")
    for ref in evidence:
        check(ref)
    result = dict(schema="native_qfo_swiss_sequence_strata_v1", source=record(__file__), binding=binding_ref,
        strata=strata_ref, snapshot=snapshot_ref, plan=plan_ref, checked_inputs=evidence,
        rows=rows, differences=differences, family_rows=family_rows, raw_relations_checked=2 * len(truth_anchor),
        input_identity_agreement="78 recorded FASTA refs in original feature/plan/conversions; no fresh proteome rehash",
        outputs=[record(output / name) for name in ("scores.tsv", "TABLE.md")],
        new_uncertainty=False, new_accuracy_or_resource_admission=False, independent_confirmation=False,
        publication_ready=False, limitations=[
            "Unchanged input-only bins, native raw counts re-read; descriptive development-exposed projections only.",
            "No new resampling/intervals; differences and overlapping bins do not establish causal subgroup effects.",
            "Direct report/raw/source checks do not repeat transitive admission or whole-proteome extraction.",
            "Feature meanings retain original limits; no total-HMM control or calibrated evolutionary distance.",
            "Failed R1 timing remains ineligible; export resources are shared-host postprocessing, not inference timing."])
    with (output / "report.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("binding", "strata", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--binding-sha256", required=True)
    args = parser.parse_args()
    result = export(args.binding, args.binding_sha256, args.strata, args.output)
    print(json.dumps({k: len(result[k]) for k in ("rows", "differences", "family_rows")}))
