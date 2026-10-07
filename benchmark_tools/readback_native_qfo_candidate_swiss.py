"""Independent raw-pair and rational readback of the native candidate contrast."""

import argparse
from collections import Counter
import csv
from fractions import Fraction
import gzip
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.readback_native_qfo_swiss_transitions import HEADER, LABELS, record, require

CELLS = ("p0_c0_r0", "p0_c1_r0")
METRICS = ("F1", "PPV", "TPR")


def raw(path, families):
    counts = {family: Counter({label: 0 for label in LABELS}) for family in families}
    members = {family: set() for family in families}
    decisions = {}
    with gzip.open(path, "rt", newline="") as stream:
        require(stream.readline().rstrip("\r\n") == HEADER, "Changed raw header")
        for row in csv.reader(stream, delimiter="\t"):
            require(len(row) == 4, "Invalid raw row width")
            family, a, b, label = row
            require(family in counts and a and b and a != b and label in LABELS, "Invalid raw relation")
            key = (family, min(a, b), max(a, b))
            require(key not in decisions, "Duplicate raw relation")
            decisions[key] = label
            counts[family][label] += 1
            members[family].update((a, b))
    require(all(sum(c.values()) for c in counts.values()), "Empty raw family")
    return counts, members, decisions


def statistics(counts, families):
    require(families and len(set(families)) == len(families), "Invalid family inventory")
    points = {}
    for family in families:
        row = counts[family]
        require(set(row) == set(LABELS) and all(type(v) is int and v >= 0 for v in row.values())
                and sum(row.values()) > 0, "Invalid confusion counts")
        tp, fp, fn = (Fraction(row[label], 2) + 1 for label in ("TP", "FP", "FN"))
        p, r = tp / (tp + fp), tp / (tp + fn)
        points[family] = dict(F1=2 * p * r / (p + r), PPV=p, TPR=r)
    p, r = (sum(row[m] for row in points.values()) / len(families) for m in ("PPV", "TPR"))
    return points, dict(F1=2 * p * r / (p + r), PPV=p, TPR=r)


def check_cell(cell, counts, members, families):
    require([r["family"] for r in cell["families"]] == families, "Changed family order")
    points, aggregate = statistics(counts, families)
    for row in cell["families"]:
        family = row["family"]
        require(row["counts_without_prior"] == dict(counts[family])
                and row["represented_genes"] == sorted(members[family]), "Changed family counts/members")
        require(set(row["statistics_with_prior"]) == set(METRICS)
                and all(math.isclose(float(points[family][m]), row["statistics_with_prior"][m],
                                     rel_tol=0, abs_tol=1e-12) for m in METRICS), "Changed family statistics")
    require(set(cell["aggregate"]) == set(METRICS) and all(
        math.isclose(float(aggregate[m]), cell["aggregate"][m], rel_tol=0, abs_tol=1e-12)
        for m in METRICS), "Changed macro statistic")
    return points, aggregate


def check_count_copy(original, supplied, parsed):
    require(record(original["path"]) == original
            and record(supplied["path"]) == supplied
            and all(original[k] == supplied[k] for k in ("bytes", "sha256"))
            and json.loads(Path(original["path"]).read_text()) == parsed,
            "Changed bootstrap count copy")


def verify(audit_path, audit_sha, binding_path, binding_sha):
    audit_ref, binding_ref = record(audit_path), record(binding_path)
    require(audit_ref["sha256"] == audit_sha and binding_ref["sha256"] == binding_sha,
            "Changed supplied report")
    audit = json.loads(Path(audit_path).read_text())
    binding = json.loads(Path(binding_path).read_text())
    require(audit["schema"] == "native_qfo_swiss_family_count_audit_v1"
            and audit["status"] == "supplied_native_swiss_family_counts_verified"
            and audit["source"] == record(Path(__file__).with_name("audit_native_qfo_swiss_counts.py"))
            and [r["cell"] for r in audit["cells"]] == [CELLS[1]]
            and all(audit[k] is False for k in ("historical_intervals_attached",
                "new_accuracy_or_resource_admission", "independent_confirmation", "publication_ready"))
            and type(audit["new_bootstrap_draws"]) is int and audit["new_bootstrap_draws"] == 0,
            "Changed candidate count scope/source")
    require(binding["schema"] == "native_qfo_retained_swiss_uncertainty_binding_v1"
            and binding["source"] == record(Path(__file__).with_name("bind_native_qfo_swiss_uncertainty.py"))
            and binding["replicates_reused"] == 100000 and binding["seed_reused"] == 20260922
            and binding["multiplicity_endpoints"] == 42 and binding["alpha"] == .05
            and type(binding["new_bootstrap_draws"]) is int and binding["new_bootstrap_draws"] == 0
            and all(binding[k] is False for k in ("new_accuracy_or_resource_admission",
                "independent_confirmation", "publication_ready")), "Changed retained uncertainty scope")
    require(binding["retained_counts"] == audit["retained_counts"]
            and binding["bound_cells"][CELLS[1]]["count_audit"] == audit_ref,
            "Candidate audit binding differs")
    checked = [audit_ref, binding_ref, audit["source"], binding["source"], *audit["helpers"],
        *audit["checked_inputs"], *binding["evidence"], binding["retained_counts"], binding["bootstrap"],
        record(Path(__file__).with_name("readback_native_qfo_swiss_transitions.py"))]
    for ref in checked:
        require(record(ref["path"]) == ref, "Changed direct evidence")
    retained = json.loads(Path(audit["retained_counts"]["path"]).read_text())
    frozen = json.loads(Path(binding["bootstrap"]["path"]).read_text())
    check_count_copy(frozen["counts"], audit["retained_counts"], retained)
    checked.append(frozen["counts"])
    require(frozen["replicates"] == 100000
            and frozen["seed"] == 20260922 and frozen["multiplicity_endpoints"] == 42,
            "Changed bootstrap provenance")
    families = audit["families"]
    require(families == retained["families"] == binding["families"] and len(families) == 18
            and len(set(families)) == 18 and audit["reference_relation_count"] == 10765,
            "Changed reference universe")
    baseline_ref = binding["bound_cells"][CELLS[0]]["count_audit"]
    require(baseline_ref in checked, "Unchecked baseline audit")
    baseline_audit = json.loads(Path(baseline_ref["path"]).read_text())
    selected = [r for r in baseline_audit["cells"] if r["cell"] == CELLS[0]]
    require(len(selected) == 1, "Missing baseline cell")
    historical = {r["cell"]: r for r in retained["cells"]}
    cells = (selected[0], audit["cells"][0])
    points, truth, raw_rows = {}, None, 0
    for index, cell in zip((6, 8), cells):
        name = cell["cell"]
        require(type(cell["index"]) is int and cell["index"] == index
                and type(cell["native_job_id"]) is int and cell["native_job_id"] > 0
                and binding["bound_cells"][name]["status"] == "native_records_matched"
                and cell["retained_family_records_identical"] is True
                and cell["retained_aggregate_identical"] is True,
                "Unmatched or wrong native cell")
        require(cell["raw_file"] in checked, "Unchecked native raw file")
        counts, members, decisions = raw(cell["raw_file"]["path"], families)
        family_points, macro = check_cell(cell, counts, members, families)
        require(cell["families"] == historical[name]["families"]
                and cell["aggregate"] == historical[name]["aggregate"], "Changed retained family records")
        current_truth = {k: label in ("TP", "FN") for k, label in decisions.items()}
        require(len(decisions) == 10765 and (truth is None or truth == current_truth), "Changed raw truth/coverage")
        truth = current_truth
        require(math.isclose(float(macro["F1"]), cell["native_endpoint_f1"], rel_tol=0, abs_tol=5e-8),
                "Changed native decimal F1")
        points[name] = (family_points, macro)
        raw_rows += len(decisions)
        if name == CELLS[1]:
            cached_ref = historical[name]["raw_file"]
            checked.append(cached_ref)
            require(record(cached_ref["path"]) == cached_ref, "Changed retained raw file")
            cached_counts, cached_members, cached_decisions = raw(cached_ref["path"], families)
            require((counts, members, decisions) == (cached_counts, cached_members, cached_decisions),
                    "Candidate scored pair decisions differ from retained raw evidence")
            raw_rows += len(cached_decisions)
    selected = [r for r in binding["contrasts"] if r["name"] == "C_at_P0_R0"]
    original = [r for r in frozen["comparisons"] if r["name"] == "C_at_P0_R0"]
    require(len(selected) == len(original) == 1, "Missing candidate contrast")
    contrast, old = selected[0], original[0]
    require(contrast["status"] == "native_records_matched" and not contrast["missing_cells"]
            and not contrast["differing_cells"] and contrast["candidate"] == CELLS[1]
            and contrast["reference"] == CELLS[0] and contrast["weights"] == [-1, 0, 1, 0, 0, 0, 0, 0]
            and contrast["metrics"] == old["metrics"]
            and contrast["family_differences"] == old["family_differences"], "Changed retained candidate contrast")
    for metric in METRICS:
        values = contrast["metrics"][metric]
        delta = points[CELLS[1]][1][metric] - points[CELLS[0]][1][metric]
        require(math.isclose(float(delta), values["difference"], rel_tol=0, abs_tol=1e-12), "Wrong contrast arithmetic")
        differences = [float(points[CELLS[1]][0][f][metric] - points[CELLS[0]][0][f][metric]) for f in families]
        require(values["family_wins"] == sum(d > 1e-10 for d in differences)
                and values["family_ties"] == sum(abs(d) <= 1e-10 for d in differences)
                and values["family_losses"] == sum(d < -1e-10 for d in differences), "Wrong family wins/losses")
    for ref in checked:
        require(record(ref["path"]) == ref, "Direct evidence changed during readback")
    return dict(schema="native_qfo_candidate_swiss_rational_readback_v1", audit=audit_ref, binding=binding_ref,
        source=record(__file__), checked_inputs=checked, cells=list(CELLS), families_checked=18,
        native_family_records_checked=36, raw_rows_checked=raw_rows, candidate_pair_labels_matched=10765,
        rational_macro_points={c: {k: float(v) for k, v in p[1].items()} for c, p in points.items()},
        contrast=contrast, new_bootstrap_draws=0, new_accuracy_or_resource_admission=False,
        scientific_timings_admitted=False, independent_confirmation=False, publication_ready=False,
        limitations=["Independent parser/rational arithmetic readback, not a new bootstrap or transitive scientific admission.",
            "Only the scored SwissTrees universe: identical decisions do not prove identical whole partitions or other endpoints.",
            "Retained 42-endpoint-adjusted intervals; 18 development-exposed families with exchangeability/percentile limits.",
            "No isolated timing, recovered-timing repair, tuning or independent biological confirmation."])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("audit", "binding", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("audit-sha256", "binding-sha256"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), "Output already exists")
    result = verify(args.audit, args.audit_sha256, args.binding, args.binding_sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: result[k] for k in ("families_checked", "raw_rows_checked", "candidate_pair_labels_matched")}))


if __name__ == "__main__":
    main()
