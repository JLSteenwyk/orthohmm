"""Generate bounded manuscript support material from verified event summaries."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "benchmark_tools/results"
REPORT_SHA = "a6dbe37a740ad657f9820cd20f78faf62c75e22e9a8ba28b258a886ee991ac8c"
READER_SHA = "dd399ccbf96943a67cb53061cafa2dd5f71cea3bc78538f4199a7770e2a4bf34"
COHORTS = ("TP-only", "FP-only", "mixed", "no changed scored VGNC pair")
LABELS = ("TP-only", "FP-only", "Mixed", "Unlabeled")
COLORS = ("#12786f", "#b64b37", "#715c9a", "#777777")
FEATURES = ("source_size", "target_size", "source_seed_families", "target_seed_families", "support", "margin",
    "species_overlap_fraction", "forward_hits", "reverse_hits", "forward_average", "reverse_average",
    "forward_maximum", "reverse_maximum", "forward_coverage", "reverse_coverage", "forward_normalized_support", "reverse_normalized_support")
STEM = "accepted_event_support"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def record(path):
    path = Path(path).resolve(strict=True)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1048576), b""):
            digest.update(block)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest.hexdigest())


def check(ref):
    require(record(ref["path"]) == ref, "Changed reporting evidence")


def snapshot(report, reader, report_ref):
    require(report.get("schema") == "native_qfo_candidate_accepted_support_v1"
        and report.get("status") == "descriptive_accepted_event_support_exported"
        and reader.get("schema") == "native_qfo_candidate_accepted_support_readback_v1"
        and reader.get("status") == "accepted_event_features_independently_verified"
        and reader.get("report") == report_ref and reader.get("totals") == report.get("totals")
        and reader.get("localized_summary") == report.get("localized_summary")
        and reader.get("rational_summaries_verified") is True and reader.get("exporter_or_union_kernel_imported") is False,
        "Require exact independently verified accepted-event report")
    require(reader.get("numeric_absolute_tolerance") == reader.get("numeric_relative_tolerance") == 1e-12,
            "Changed independent tolerance")
    for evidence in (report, reader):
        require(all(evidence.get(k) is False for k in ("defaults_changed", "uncertainty_admitted",
            "calibrated_confidence", "causal_mechanism_established", "publication_ready")), "Changed scientific scope")
    require(report.get("failed_r1_timing_remains_ineligible") is True, "Failed timing relabeled")
    groups = report["cohort_summaries"]
    require(len(groups) == 12 and reader.get("cohort_counts") == [
        {k: row[k] for k in ("iteration", "cohort", "events")} for row in groups], "Changed cohort inventory")
    by_key = {}
    for row in groups:
        iteration, cohort, n = row["iteration"], row["cohort"], row["events"]
        require((iteration is None or type(iteration) is int and iteration in (0, 1))
            and cohort in COHORTS and type(n) is int and n >= 0 and (iteration, cohort) not in by_key,
            "Invalid cohort or round identity")
        require(set(row["features"]) == set(FEATURES), "Partial feature inventory")
        for name, values in row["features"].items():
            require(set(values) == {"finite", "missing", "positive_infinity", "minimum", "median", "mean", "maximum"},
                    "Changed feature summary shape")
            require(all(type(values[k]) is int and values[k] >= 0 for k in ("finite", "missing", "positive_infinity"))
                and values["missing"] == 0 and values["finite"] + values["positive_infinity"] == n
                and (name == "margin" or values["positive_infinity"] == 0), "Invalid feature counts")
            numeric = [values[k] for k in ("minimum", "median", "mean", "maximum")]
            if values["finite"] == 0:
                require(all(v is None for v in numeric), "Empty finite summary was imputed")
            else:
                require(all(type(v) in (int, float) and math.isfinite(v) and v >= 0 for v in numeric)
                    and numeric[0] <= numeric[1] <= numeric[3] and numeric[0] <= numeric[2] <= numeric[3],
                    "Invalid finite summary range")
        by_key[iteration, cohort] = row
    require(set(by_key) == {(r, c) for r in (None, 0, 1) for c in COHORTS}, "Missing cohort/round cell")
    totals = report["totals"]
    require(set(totals) == {"accepted_events", "changed_pairs", "implicated_events", "direct_tp_pairs", "direct_fp_pairs",
                          "transitive_tp_pairs", "transitive_fp_pairs"}
        and all(type(v) is int and v >= 0 for v in totals.values()) and totals["accepted_events"] > 0,
        "Invalid diagnostic totals")
    require(sum(by_key[None, c]["events"] for c in COHORTS) == totals["accepted_events"]
        and sum(by_key[None, c]["events"] for c in COHORTS[:3]) == totals["implicated_events"]
        and sum(totals[k] for k in ("direct_tp_pairs", "direct_fp_pairs", "transitive_tp_pairs", "transitive_fp_pairs"))
            == totals["changed_pairs"], "Event and dependent-pair units differ")
    for cohort in COHORTS:
        require(by_key[None, cohort]["events"] == sum(by_key[r, cohort]["events"] for r in (0, 1)), "Round event counts differ")
        for name in FEATURES:
            for key in ("finite", "missing", "positive_infinity"):
                require(by_key[None, cohort]["features"][name][key] == sum(by_key[r, cohort]["features"][name][key] for r in (0, 1)),
                        "Round feature counts differ")
    paths, seen = report["localized_summary"], set()
    path_totals = {k: 0 for k in ("direct_tp_pairs", "direct_fp_pairs", "transitive_tp_pairs", "transitive_fp_pairs")}
    for row in paths:
        key = (row["iteration"], row["candidate_state"], row["connection_path"])
        require(type(key[0]) is int and key[0] in (0, 1) and key[1] in ("TP", "FP")
            and key[2] in ("direct_cross_endpoint", "transitive_union") and key not in seen
            and type(row["pairs"]) is int and row["pairs"] > 0, "Invalid path summary")
        seen.add(key)
        name = ("direct_" if key[2] == "direct_cross_endpoint" else "transitive_") + key[1].lower() + "_pairs"
        path_totals[name] += row["pairs"]
    require(all(v == totals[k] for k, v in path_totals.items()), "Path counts differ")
    return by_key


def compact(value):
    if value is None:
        return "not available"
    return f"{value:g}" if isinstance(value, int) else f"{value:.6f}"


def render(groups, output):
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 9, "svg.hashsalt": "accepted-event-support-v1",
                         "svg.fonttype": "none", "pdf.fonttype": 42}):
        figure, axes = plt.subplots(2, 2, figsize=(12.8, 8.8))
        axis = axes[0, 0]
        axis.axis("off")
        axis.set_title("A  Accepted-event units", loc="left", fontsize=12, pad=16)
        axis.text(.53, .94, "All rounds", ha="right", transform=axis.transAxes)
        axis.text(.75, .94, "Round 0", ha="right", transform=axis.transAxes)
        axis.text(.97, .94, "Round 1", ha="right", transform=axis.transAxes)
        for position, (cohort, label, color) in enumerate(zip(COHORTS, LABELS, COLORS)):
            y = .76 - position * .20
            axis.text(.01, y, label, color=color, fontsize=10, transform=axis.transAxes)
            for x, iteration in ((.53, None), (.75, 0), (.97, 1)):
                axis.text(x, y, f"{groups[iteration, cohort]['events']:,}", ha="right", fontsize=11, transform=axis.transAxes)
        for axis, metric, title, xlabel in ((axes[0, 1], "support", "B  Serialized event support", "Support (observed range and median)"),
            (axes[1, 0], "source_size", "C  Source-group size", "Source members (observed range and median)"),
            (axes[1, 1], "target_size", "D  Target-group size", "Target members (observed range and median)")):
            maximum = 1.0
            for y, (cohort, color) in enumerate(zip(COHORTS, COLORS)):
                values = groups[None, cohort]["features"][metric]
                if values["finite"]:
                    axis.plot([values["minimum"], values["maximum"]], [y, y], color=color, linewidth=1.6)
                    axis.scatter(values["median"], y, s=44, color=color, zorder=3, edgecolor="white", linewidth=.5)
                    maximum = max(maximum, values["maximum"])
                else:
                    axis.text(.02, y, "No finite observations", va="center", transform=axis.get_yaxis_transform())
            axis.set_title(title, loc="left", fontsize=12, pad=16)
            axis.set(xlim=(0, maximum * 1.05), ylim=(3.5, -.5), xlabel=xlabel)
            axis.set_yticks(range(4), LABELS)
            axis.grid(axis="x", color="#e9e9e9", linewidth=.7)
            axis.spines[["top", "right"]].set_visible(False)
        figure.suptitle("Accepted candidate events: descriptive VGNC associations", x=.075, ha="left", y=.98, fontsize=15)
        figure.text(.075, .935, "Each accepted event once. Lines = observed minimum to maximum; dots = median. These are not confidence intervals.")
        figure.text(.075, .106, "Cohorts refer only to changed scored pairs: TP-only, FP-only, mixed, or unlabeled. Unlabeled does not mean correct or true negative.")
        figure.text(.075, .076, "Transitive pairs have no assigned direct event. Accepted-event group support does not establish a direct protein-pair HMM hit.")
        figure.text(.075, .046, "Development-exposed; membership, group size and reference scope confound differences. No threshold, confidence or causal inference.")
        figure.subplots_adjust(left=.10, right=.97, top=.83, bottom=.21, wspace=.40, hspace=.65)
        for suffix in ("pdf", "png", "svg"):
            figure.savefig(output / f"{STEM}.{suffix}", dpi=200, metadata={"Creator": "OrthoHMM benchmark workflow"})
        plt.close(figure)


def manuscript(report, groups):
    totals = report["totals"]
    lines = ["# Accepted Candidate Events: VGNC Support Diagnostic", "", "## Methods", "",
        "This supplementary diagnostic joins the original accepted candidate-event trace to already-verified changed VGNC pair paths. "
        "It reuses the complete partition and original-protein identity evidence, without reexecuting inference, aliases, grouping or scoring. "
        "Every accepted event contributes once to feature summaries, irrespective of the number of dependent scored pairs it introduces. "
        "TP-only events have recovered asserted TPs but no added scored FPs; FP-only events have the converse; mixed events affect both. "
        "Events without a changed scored VGNC pair are unlabeled, not true negatives or biologically correct events.", "",
        "The exporter validates every accepted event, its round, nonoverlapping memberships, declared sizes, unique semantic identity and numeric features. "
        "All direct pair joins check named endpoints and the first-connection round. Transitive paths retain no direct event or assigned direct-event feature. "
        "A separate standard-library reader imports neither exporter nor union kernel and reconstructs joins, tables and rational summary statistics. "
        "Identities and integer counts agree exactly; numerical summaries use absolute and relative tolerance 1e-12. "
        "All source, input and output byte bindings are checked before and after each original diagnostic command.", "",
        "Feature summaries contain finite minimum, median, mean and maximum, missing counts and positive-infinity counts. "
        "All required features are present. The serialized positive-infinity margin denotes no positive alternative; it is counted separately and excluded "
        "from finite summaries, not replaced with a large number. No cutoff fitting, significance tests, bootstrap intervals or default tuning is performed.", "",
        "## Results", "", f"The trace contains {totals['accepted_events']:,} accepted events. The ledger retains all {totals['changed_pairs']:,} "
        f"changed pair paths; {totals['implicated_events']:,} distinct events directly implicate changed scored pairs.", "",
        "| Cohort | Accepted events | Round 0 | Round 1 | Median support | Finite median margin | Positive-infinity margins |",
        "|---|---:|---:|---:|---:|---:|---:|"]
    for cohort, label in zip(COHORTS, LABELS):
        row = groups[None, cohort]
        features = row["features"]
        lines.append(f"| {label} | {row['events']:,} | {groups[0, cohort]['events']:,} | {groups[1, cohort]['events']:,} | "
            f"{compact(features['support']['median'])} | {compact(features['margin']['median'])} | {features['margin']['positive_infinity']:,} |")
    lines += ["", f"Direct paths account for {totals['direct_tp_pairs']:,} recovered TPs and {totals['direct_fp_pairs']:,} added scored FPs. "
        f"The remaining {totals['transitive_tp_pairs']:,} recovered TPs and {totals['transitive_fp_pairs']:,} added FPs have transitive paths. "
        "These are dependent pair counts, not independent event units.", "",
        f"![Accepted-event counts, serialized support and group sizes.]({STEM}.png)", "",
        "Figure S1. Counts use each accepted event once. Lines show observed minimum to maximum and dots show medians, not confidence intervals. "
        "Feature ranges overlap; group-level support does not imply calibrated assignment confidence. The source and target sizes are retained serialized memberships.", "",
        "## Interpretation And Limitations", "",
        "The descriptive cohort support ranges overlap, and mixed events affect both recovered TPs and added FPs. "
        "Higher median support in one cohort is not a validated decision boundary or a demonstrated accuracy gain. "
        "Membership, group size, reference scope and conditioning on accepted events confound differences. "
        "This analysis does not inspect rejected alternatives, recompute eligibility/support criteria, establish a direct pairwise HMM hit, "
        "isolate a biological cause, provide independent VGNC uncertainty or supply independent validation. "
        "It neither changes defaults nor demonstrates superiority over OrthoFinder or any other tool.", "",
        "These analyses use development-exposed QfO observations. Failed R-on timing remains ineligible. "
        "Diagnostic postprocessing observations were collected on a shared Threadripper and are not native inference timings or isolated-performance estimates. "
        "The full publication goal remains incomplete. This companion does not modify or supersede the retained main manuscript, rc5 delivery, "
        "scientific admissions or earlier failed outcomes.", "", "## Claim-To-Evidence Checklist", "",
        "| Claim | Scope | Evidence |", "|---|---|---|"]
    for text, scope, target in (("All changed pair identities and paths retained", "Descriptive changed-pair ledger", "annotated_pairs.tsv"),
        ("Event units counted once, including mixed and unlabeled", "Accepted events only", "event_cohorts.tsv"),
        ("Finite summaries separate positive-infinity margin", "Serialized features, not calibrated confidence", "feature_summaries.tsv"),
        ("Direct and transitive changes distinguished", "No invented direct event for transitive pairs", "pair_paths.tsv")):
        lines.append(f"| {text} | {scope} | [{target}]({target}) |")
    lines += ["", "Machine-readable [claim bindings](claims.json), [original diagnostic report](report.json), "
        "[independent verification](readback.json), [unique implicated events](implicated_events.tsv), "
        "[execution receipt](execution.json) and [prospective protocol](protocol.md) accompany this supplement. "
        "Presence and hash checking do not imply full inference reproduction, transitive dependency closure, permissions or archival readiness.", ""]
    return "\n".join(lines)


def write_table(path, header, rows):
    with path.open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header); writer.writerows(rows)


def generate(report_path, reader_path, output):
    output = Path(output).resolve()
    require(not output.exists() and not output.is_symlink(), "Fresh output namespace required")
    report_ref, reader_ref = record(report_path), record(reader_path)
    require(report_ref["sha256"] == REPORT_SHA and reader_ref["sha256"] == READER_SHA, "Selected diagnostic anchor differs")
    report, reader = json.loads(Path(report_path).read_text()), json.loads(Path(reader_path).read_text())
    groups = snapshot(report, reader, report_ref)
    evidence = [report_ref, reader_ref, report["annotated_pairs"], report["implicated_events"], report["source"],
        report["reader_source"], report["protocol"], record(RESULTS / "native_qfo_candidate_support_execution_20261006_v1.json"), record(__file__)]
    require(reader["source"] == report["reader_source"], "Changed independent reader source")
    for ref in evidence:
        check(ref)
    output.mkdir(parents=True, exist_ok=False)
    copies = {"report.json": report_ref, "readback.json": reader_ref, "annotated_pairs.tsv": report["annotated_pairs"],
              "implicated_events.tsv": report["implicated_events"], "execution.json": evidence[7], "protocol.md": report["protocol"]}
    for name, ref in copies.items():
        with (output / name).open("xb") as handle:
            handle.write(Path(ref["path"]).read_bytes())
    write_table(output / "event_cohorts.tsv", ["iteration", "cohort", "events"],
        [["all" if r is None else r, c, groups[r, c]["events"]] for r in (None, 0, 1) for c in COHORTS])
    write_table(output / "feature_summaries.tsv", ["iteration", "cohort", "events", "feature", "finite", "missing", "positive_infinity",
        "minimum", "median", "mean", "maximum"], [["all" if r is None else r, c, groups[r, c]["events"], f,
        *[groups[r, c]["features"][f][k] for k in ("finite", "missing", "positive_infinity", "minimum", "median", "mean", "maximum")]]
        for r in (None, 0, 1) for c in COHORTS for f in FEATURES])
    write_table(output / "pair_paths.tsv", ["iteration", "candidate_state", "connection_path", "pairs"],
        [[row[k] for k in ("iteration", "candidate_state", "connection_path", "pairs")] for row in report["localized_summary"]])
    claims = dict(schema="native_qfo_support_supplement_claims_v1", report=report_ref, reader=reader_ref,
        unit="accepted event once; pair counts dependent", claims=[
            dict(id="complete_paths", status="independently_verified_descriptive", evidence="annotated_pairs.tsv", pairs=report["totals"]["changed_pairs"]),
            dict(id="event_units", status="independently_verified_descriptive", evidence="event_cohorts.tsv", events=report["totals"]["accepted_events"]),
            dict(id="feature_summaries", status="independently_verified_descriptive", evidence="feature_summaries.tsv", nonfinite_policy="positive infinity margin separate"),
            dict(id="path_semantics", status="independently_verified_descriptive", evidence="pair_paths.tsv", transitive_direct_event_assigned=False)],
        not_established=["biological correctness of unlabeled events", "direct protein-pair HMM hits", "calibrated confidence", "causal error mechanisms",
                         "valid VGNC confidence intervals", "independent confirmation", "optimal thresholds", "OrthoFinder superiority", "publication readiness"])
    with (output / "claims.json").open("x") as handle:
        json.dump(claims, handle, indent=2, sort_keys=True, allow_nan=False); handle.write("\n")
    with (output / "supplement.md").open("x") as handle:
        handle.write(manuscript(report, groups))
    render(groups, output)
    for ref in evidence:
        check(ref)
    result = dict(schema="native_qfo_support_supplement_v1", status="descriptive_supplement_generated", source=record(__file__), evidence=evidence,
        copied_outputs={name: record(output / name) for name in copies}, outputs=[record(p) for p in sorted(output.iterdir())],
        figure_unit="accepted event once", figure_lines="observed minimum to maximum, NOT confidence intervals", figure_dots="median",
        plotted_features=["support", "source_size", "target_size"], cohort_rows=12, feature_rows=12 * len(FEATURES),
        python_version=sys.version, matplotlib_version=matplotlib.__version__, original_trace_or_partition_replayed=False,
        new_accuracy_or_resource_admission=False, new_bootstrap_draws=0, defaults_changed=False, publication_ready=False, visual_review_complete=False)
    with (output / "generation.json").open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False); handle.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=RESULTS / "native_qfo_candidate_support_20261006_v1/report.json")
    parser.add_argument("--readback", type=Path, default=RESULTS / "native_qfo_candidate_support_readback_20261006_v1.json")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = generate(args.report, args.readback, args.output)
    print(json.dumps({k: result[k] for k in ("status", "cohort_rows", "feature_rows", "plotted_features")}))
