"""Independent pair/rational readback of the allocated profile-refinement effect."""

import argparse
import json
import math
from pathlib import Path

from benchmark_tools import readback_native_qfo_candidate_swiss as reader

CELLS = ("p0_c0_r1", "p1_c0_r1")
RELATIONS = 10765
METRICS = reader.METRICS
record, require = reader.record, reader.require


def verify(audit_path, audit_sha, binding_path, binding_sha, prior_path, prior_sha):
    refs = [record(p) for p in (audit_path, binding_path, prior_path)]
    require([r["sha256"] for r in refs] == [audit_sha, binding_sha, prior_sha], "Changed supplied report")
    audit, binding, prior = [json.loads(Path(r["path"]).read_text()) for r in refs]
    require(audit["schema"] == "allocated_native_qfo_swiss_family_count_audit_v1"
        and audit["status"] == "supplied_allocated_native_swiss_family_counts_verified"
        and audit["source"] == record(Path(__file__).with_name("audit_allocated_native_qfo_swiss_counts.py"))
        and type(audit["selected_index"]) is int and audit["selected_index"] == 10
        and [r["cell"] for r in audit["cells"]] == [CELLS[1]], "Changed allocated count scope/source")
    require(binding["schema"] == "allocated_native_qfo_retained_swiss_uncertainty_binding_v1"
        and binding["source"] == record(Path(__file__).with_name("bind_allocated_native_qfo_swiss_uncertainty.py"))
        and binding["replicates_reused"] == 100000 and binding["seed_reused"] == 20260922
        and binding["multiplicity_endpoints"] == 42 and binding["alpha"] == .05,
        "Changed allocated uncertainty scope")
    require(prior["schema"] == "native_qfo_retained_swiss_uncertainty_binding_v1"
        and prior["source"] == record(Path(__file__).with_name("bind_native_qfo_swiss_uncertainty.py"))
        and prior["bootstrap"] == binding["bootstrap"] and prior["retained_counts"] == binding["retained_counts"],
        "Changed prior binding provenance")
    for report in (audit, binding, prior):
        require(type(report["new_bootstrap_draws"]) is int and report["new_bootstrap_draws"] == 0
            and all(report[k] is False for k in (
                "new_accuracy_or_resource_admission", "independent_confirmation", "publication_ready")),
            "Unjustified admission or readiness claim")
    require(audit["historical_intervals_attached"] is False
        and binding["retained_counts"] == audit["retained_counts"]
        and binding["bound_cells"][CELLS[1]]["count_audit"] == refs[0], "Allocated audit binding differs")
    checked = [*refs, audit["source"], binding["source"], prior["source"], *audit["helpers"],
        *audit["checked_inputs"], *binding["helpers"], *binding["evidence"],
        binding["retained_counts"], binding["bootstrap"], record(reader.__file__),
        record(Path(__file__).with_name("readback_native_qfo_swiss_transitions.py"))]
    for ref in checked:
        require(record(ref["path"]) == ref, "Changed direct evidence")
    retained = json.loads(Path(audit["retained_counts"]["path"]).read_text())
    frozen = json.loads(Path(binding["bootstrap"]["path"]).read_text())
    reader.check_count_copy(frozen["counts"], audit["retained_counts"], retained)
    checked.append(frozen["counts"])
    require(frozen["replicates"] == 100000 and frozen["seed"] == 20260922
        and frozen["multiplicity_endpoints"] == 42, "Changed bootstrap provenance")
    families = audit["families"]
    require(families == retained["families"] == binding["families"] and len(families) == 18
        and len(set(families)) == 18 and audit["reference_relation_count"] == RELATIONS,
        "Changed reference universe")
    baseline_ref = binding["bound_cells"][CELLS[0]]["count_audit"]
    require(baseline_ref in checked, "Unchecked reference audit")
    baseline_audit = json.loads(Path(baseline_ref["path"]).read_text())
    require(baseline_audit["schema"] == "recovered_native_qfo_swiss_family_count_audit_v1"
        and baseline_audit["source"] == record(Path(__file__).with_name("audit_recovered_native_qfo_swiss_counts.py")),
        "Changed recovered reference source")
    selected = [r for r in baseline_audit["cells"] if r["cell"] == CELLS[0]]
    require(len(selected) == 1, "Missing reference cell")
    reference = selected[0]
    require(reference["resources"] is None and reference["timing_admitted"] is False
        and reference["timing_eligible"] is False, "Readback relabels failed timing")
    historical = {r["cell"]: r for r in retained["cells"]}
    points, truth, raw_rows = {}, None, 0
    for index, cell in zip((7, 10), (reference, audit["cells"][0])):
        name = cell["cell"]
        require(type(cell["index"]) is int and cell["index"] == index
            and type(cell["native_job_id"]) is int and cell["native_job_id"] > 0
            and binding["bound_cells"][name]["status"] == "native_records_matched"
            and cell["retained_family_records_identical"] is True
            and cell["retained_aggregate_identical"] is True, "Unmatched or wrong native cell")
        require(cell["raw_file"] in checked, "Unchecked native raw file")
        counts, members, decisions = reader.raw(cell["raw_file"]["path"], families)
        family_points, macro = reader.check_cell(cell, counts, members, families)
        require(cell["families"] == historical[name]["families"]
            and cell["aggregate"] == historical[name]["aggregate"], "Changed retained family records")
        current_truth = {key: label in ("TP", "FN") for key, label in decisions.items()}
        require(len(decisions) == RELATIONS and (truth is None or truth == current_truth),
            "Changed raw truth/coverage")
        truth = current_truth
        require(math.isclose(float(macro["F1"]), cell["native_endpoint_f1"], rel_tol=0, abs_tol=5e-8),
            "Changed native decimal F1")
        points[name] = (family_points, macro)
        raw_rows += len(decisions)
        if name == CELLS[1]:
            cached_ref = historical[name]["raw_file"]
            checked.append(cached_ref)
            require(record(cached_ref["path"]) == cached_ref, "Changed retained raw file")
            cached = reader.raw(cached_ref["path"], families)
            require((counts, members, decisions) == cached, "Profile scored pair decisions differ from retained")
            raw_rows += len(cached[2])
    new = {r["name"]: r for r in binding["contrasts"]}
    require(len(new) == len(binding["contrasts"]) == 14, "Changed contrast inventory")
    for row in prior["contrasts"]:
        if row["status"] == "native_records_matched":
            require(new[row["name"]] == row, "Prior matched contrast changed")
    contrast = new["P_at_C0_R1"]
    original = [r for r in frozen["comparisons"] if r["name"] == contrast["name"]]
    require(len(original) == 1 and contrast["status"] == "native_records_matched"
        and not contrast["missing_cells"] and not contrast["differing_cells"]
        and contrast["candidate"] == CELLS[1] and contrast["reference"] == CELLS[0]
        and contrast["weights"] == [0, -1, 0, 0, 0, 1, 0, 0]
        and contrast["metrics"] == original[0]["metrics"]
        and contrast["family_differences"] == original[0]["family_differences"], "Changed profile contrast")
    for metric in METRICS:
        values = contrast["metrics"][metric]
        delta = points[CELLS[1]][1][metric] - points[CELLS[0]][1][metric]
        require(math.isclose(float(delta), values["difference"], rel_tol=0, abs_tol=1e-12),
            "Wrong contrast arithmetic")
        differences = [float(points[CELLS[1]][0][f][metric] - points[CELLS[0]][0][f][metric]) for f in families]
        require(values["family_wins"] == sum(d > 1e-10 for d in differences)
            and values["family_ties"] == sum(abs(d) <= 1e-10 for d in differences)
            and values["family_losses"] == sum(d < -1e-10 for d in differences), "Wrong family wins/losses")
        require(all(math.isclose(d, row[metric], rel_tol=0, abs_tol=1e-12)
            for d, row in zip(differences, contrast["family_differences"])), "Wrong family differences")
    for ref in checked:
        require(record(ref["path"]) == ref, "Direct evidence changed during readback")
    return dict(schema="allocated_native_qfo_profile_swiss_rational_readback_v1",
        audit=refs[0], binding=refs[1], prior_binding=refs[2], source=record(__file__), checked_inputs=checked,
        cells=list(CELLS), families_checked=18, native_family_records_checked=36, raw_rows_checked=raw_rows,
        profile_pair_labels_matched=RELATIONS,
        rational_macro_points={c: {k: float(v) for k, v in p[1].items()} for c, p in points.items()},
        contrast=contrast, prior_matched_contrasts_unchanged=True,
        prior_matched_contrasts_checked=sum(r["status"] == "native_records_matched" for r in prior["contrasts"]),
        new_bootstrap_draws=0,
        new_accuracy_or_resource_admission=False, scientific_timings_admitted=False,
        independent_confirmation=False, publication_ready=False, limitations=[
            "Separate CSV/gzip parser and exact Fraction arithmetic, not a new bootstrap or admission.",
            "Intervals checked against frozen output, not independently recalculated or statistically certified.",
            "Only scored SwissTrees pairs; no whole-partition equivalence or other-endpoint uncertainty.",
            "Development exposure, 18-family exchangeability and approximate percentile limits remain.",
            "Initial HMM search stays on; conditional downstream refinement is not a total-HMM control."])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("audit", "binding", "prior-binding", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("audit-sha256", "binding-sha256", "prior-binding-sha256"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), "Output already exists")
    result = verify(args.audit, args.audit_sha256, args.binding, args.binding_sha256,
        args.prior_binding, args.prior_binding_sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: result[k] for k in ("families_checked", "raw_rows_checked", "profile_pair_labels_matched")}))


if __name__ == "__main__":
    main()
