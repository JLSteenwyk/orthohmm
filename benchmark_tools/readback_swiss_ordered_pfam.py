"""Independent all-pair order checks and rational retained-count/table readback."""

import argparse
import csv
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess

BINS = ("all", "all_members_usable_same_signature", "all_members_usable_multiple_signatures", "some_members_unusable")
CELLS = ("p0_c0_r0", "p0_c0_r1", "p0_c1_r0")
METRICS = ("F1", "PPV", "TPR")
PINS = dict(annotations="d5e269158c2fb603acc0805a1c140a88758342d7b8e32b3094750ade22fbc06c",
            domain_report="d3235c712e711f005b3182a3621c0079f8a85c3c7abd5db016d5e8ea95884cd6",
            domain_reader="9e77b67ad97e6824e9442691be7a05ba1552d62ffe696ced6034a996ede1945d",
            alignments="ceb4ca3675836e1815fcc1f14c29c8a50947633a8a5e8dd6c67477bdd853065c",
            protocol="c58b217142bf00765cb8733819c4c93d06299eb38eb2f318d384b3fcd9395272",
            counts="55ed72ce771c756424877cfa412da329ff104efb7dca4b35113b0d80703a59e5",
            counts_reader="6f3eef63c0e96e1d85805fcfb70ed9b66fbb1db76a08ca2a105483422a0f6969")
FALSE_FLAGS = ("new_uncertainty", "independent_confirmation", "new_accuracy_or_resource_admission",
               "publication_ready", "scientific_timings_admitted", "raw_scorer_repeated", "alignment_repeated")


def demand(condition, message):
    if not condition:
        raise ValueError(message)


def record(path):
    path = Path(path).resolve(strict=True)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def scope(report):
    demand(all(report[k] is False for k in FALSE_FLAGS) and type(report["new_bootstrap_draws"]) is int
           and report["new_bootstrap_draws"] == 0, "Inflated scientific scope")


def committed_source(ref, repo, commit, name):
    demand(record(ref["path"]) == ref and Path(ref["path"]).read_bytes() ==
           subprocess.check_output(["git", "show", commit + ":benchmark_tools/" + name], cwd=repo), "Unfrozen source")


def protein_descriptor(annotation, length):
    demand(type(length) is int and length > 0 and type(annotation["length"]) is int
           and annotation["length"] == length, "Annotation/sequence length mismatch")
    raw = annotation["ordered_pfam_instances"]
    demand(isinstance(raw, list), "Invalid annotation list")
    domains, seen = {}, set()
    for item in raw:
        demand(set(item) == {"start", "end", "domain"} and type(item["start"]) is int
               and type(item["end"]) is int and 0 <= item["start"] <= item["end"] <= length
               and isinstance(item["domain"], str) and item["domain"].startswith("pfam_"), "Invalid retained instance")
        key = (item["start"], item["end"], item["domain"])
        demand(key not in seen, "Duplicate retained instance")
        seen.add(key)
        domains[item["domain"]] = domains.get(item["domain"], 0) + 1
    demand(annotation["pfam_types"] == sorted(domains) and type(annotation["pfam_type_count"]) is int
           and annotation["pfam_type_count"] == len(domains) and type(annotation["pfam_instance_count"]) is int
           and annotation["pfam_instance_count"] == len(raw)
           and annotation["has_repeated_pfam_type"] is any(v > 1 for v in domains.values()), "Inconsistent annotation summary")
    # Deliberately all pairs, not the primary implementation's adjacent check.
    overlaps = any(max(a["start"], b["start"]) <= min(a["end"], b["end"])
                   for a, b in itertools.combinations(raw, 2))
    ambiguous = overlaps or any(r["start"] == r["end"] for r in raw)
    state = "zero_pfam" if not raw else "order_ambiguous" if ambiguous else "usable"
    ordered = sorted(raw, key=lambda r: (r["start"], r["end"], r["domain"]))
    return dict(length=length, status=state, ordered_instances=ordered,
                ordered_signature=[r["domain"] for r in ordered] if state == "usable" else None,
                type_multiset=[[d, domains[d]] for d in sorted(domains)])


def family_descriptor(members, genes):
    demand(members == sorted(set(members)) and bool(members), "Invalid member list")
    usable = [g for g in members if genes[g]["status"] == "usable"]
    signatures = set(tuple(genes[g]["ordered_signature"]) for g in usable)
    matched_pairs = [(a, b) for a, b in itertools.combinations(usable, 2)
                     if genes[a]["type_multiset"] == genes[b]["type_multiset"]]
    category = BINS[1] if len(signatures) == 1 else BINS[2]
    if len(usable) != len(members):
        category = BINS[3]
    return dict(members=members, member_count=len(members), usable_count=len(usable),
                ambiguous_count=len([g for g in members if genes[g]["status"] == "order_ambiguous"]),
                zero_pfam_count=len([g for g in members if genes[g]["status"] == "zero_pfam"]),
                distinct_usable_signatures=len(signatures), same_multiset_comparable_pairs=len(matched_pairs),
                same_multiset_order_discordant_pairs=len([1 for a, b in matched_pairs
                    if genes[a]["ordered_signature"] != genes[b]["ordered_signature"]]), stratum=category)


def feature_readback(report, inventory, alignments, checked):
    from Bio import SeqIO
    demand(report["schema"] == "swiss_ordered_pfam_features_v1" and report["status"] == "features_constructed_unverified"
           and report["prediction_statistics_evaluated"] is False
           and report["original_annotation_extraction_repeated"] is False, "Incorrect feature scope")
    memberships = alignments["memberships"]
    expected = [g for f in sorted(memberships) for g in memberships[f]]
    demand(report["memberships"] == memberships and len(memberships) == 18
           and len(expected) == len(set(expected)) == 563 and set(expected) == set(inventory["genes"]), "Changed feature population")
    demand([r["family"] for r in alignments["runs"]] == sorted(memberships), "Incomplete alignments")
    lengths = {}
    for run in alignments["runs"]:
        ref = run["alignment"]
        demand(record(ref["path"]) == ref, "Changed retained alignment")
        checked.append(ref)
        records = list(SeqIO.parse(ref["path"], "fasta"))
        demand(sorted(r.id for r in records) == memberships[run["family"]]
               and all(len(r.seq) == run["columns"] for r in records), "Changed alignment membership/columns")
        for sequence in records:
            demand(sequence.id not in lengths, "Duplicate aligned protein")
            lengths[sequence.id] = len(str(sequence.seq).replace("-", ""))
    genes = {g: protein_descriptor(inventory["genes"][g], lengths[g]) for g in sorted(expected)}
    families = {f: family_descriptor(memberships[f], genes) for f in sorted(memberships)}
    bins = {k: sorted(f for f in families if k == "all" or families[f]["stratum"] == k) for k in BINS}
    demand(report["genes"] == genes and report["families"] == families and report["bins"] == bins,
           "Incorrect protein order, family diagnostic or fixed bin")
    return dict(genes_checked=len(genes), families_checked=len(families), bins=bins,
                protein_states={s: sum(g["status"] == s for g in genes.values())
                                for s in ("usable", "order_ambiguous", "zero_pfam")},
                same_multiset_comparable_pairs=sum(f["same_multiset_comparable_pairs"] for f in families.values()),
                same_multiset_order_discordant_pairs=sum(f["same_multiset_order_discordant_pairs"] for f in families.values()))


def point(counts):
    demand(set(counts) == {"TP", "FP", "FN", "TN"} and all(type(v) is int and v >= 0 for v in counts.values())
           and sum(counts.values()) > 0, "Invalid integer counts")
    tp, fp, fn = (Fraction(counts[k] + 2, 2) for k in ("TP", "FP", "FN"))
    p, r = tp / (tp + fp), tp / (tp + fn)
    return dict(F1=2 * p * r / (p + r), PPV=p, TPR=r)


def macro(points):
    if not points:
        return dict.fromkeys(METRICS)
    p, r = (sum(row[m] for row in points) / len(points) for m in ("PPV", "TPR"))
    return dict(F1=2 * p * r / (p + r), PPV=p, TPR=r)


def numeric(actual, expected):
    demand(actual is None if expected is None else type(actual) in (int, float) and math.isfinite(actual)
           and math.isclose(actual, float(expected), abs_tol=1e-12, rel_tol=0), "Incorrect numeric/empty statistic")


def score_readback(report, counts, features):
    demand(report["schema"] == "native_swiss_ordered_pfam_projection_v1" and report["status"] == "projection_constructed_unverified"
           and report["family_rows"] == counts["family_rows"] and report["memberships"] == counts["memberships"] == features["memberships"]
           and report["cells"] == counts["cells"] and report["bins"] == features["bins"], "Changed retained counts, cells or bins")
    demand(report["cells"][1]["timing_eligible"] is False and report["cells"][1]["timing_admitted"] is False, "Failed timing promoted")
    families = sorted(counts["memberships"])
    bins = features["bins"]
    demand(set(bins) == set(BINS) and bins["all"] == families
           and all(v == sorted(set(v)) for v in bins.values())
           and sorted(f for k in BINS[1:] for f in bins[k]) == families, "Invalid bin partition")
    demand([(r["cell"], r["family"]) for r in counts["family_rows"]] == [(c, f) for c in CELLS for f in families], "Incomplete count rows")
    index = {}
    for row in counts["family_rows"]:
        value = point(row["counts_without_prior"])
        index[row["cell"], row["family"]] = value
        for metric in METRICS:
            numeric(row[metric], value[metric])
    rows, differences = [], []
    for name in BINS:
        members = bins[name]
        values = {c: macro([index[c, f] for f in members]) for c in CELLS}
        for cell in CELLS:
            rows.append(dict(stratum=name, cell=cell, families=len(members), family_members=members,
                status="descriptive" if members else "empty",
                prediction_semantics="native_pair" if cell == CELLS[1] else "group_clique", **values[cell]))
        for label, candidate in (("R_at_P0_C0", CELLS[1]), ("C_at_P0_R0", CELLS[2])):
            differences.append(dict(stratum=name, contrast=label, candidate=candidate, reference=CELLS[0],
                families=len(members), family_members=members, status="descriptive" if members else "empty",
                **{m + "_pp": 100 * (values[candidate][m] - values[CELLS[0]][m]) if members else None for m in METRICS}))
    for name, computed, metrics in (("rows", rows, METRICS), ("differences", differences, tuple(m + "_pp" for m in METRICS))):
        demand(len(report[name]) == len(computed), "Incomplete projected rows")
        for actual, expected in zip(report[name], computed):
            demand(set(actual) == set(expected), "Changed projected schema")
            for field, value in expected.items():
                if field in metrics:
                    numeric(actual[field], value)
                else:
                    demand(actual[field] == value, "Changed labels, semantics or memberships")
    return rows, differences


def tables_readback(outputs, rows, differences):
    for key, expected, fields, metrics in (
        ("scores", rows, ("stratum", "cell", "families", "status", *METRICS, "prediction_semantics"), METRICS),
        ("differences", differences, ("stratum", "contrast", "families", "status", *(m + "_pp" for m in METRICS)), tuple(m + "_pp" for m in METRICS))):
        with Path(outputs[key]["path"]).open(newline="", encoding="ascii") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            demand(reader.fieldnames == list(fields), "Incorrect TSV header")
            actual_rows = list(reader)
        demand(len(actual_rows) == len(expected), "Incorrect TSV row count")
        for actual, computed in zip(actual_rows, expected):
            demand(set(actual) == set(fields), "Incorrect TSV fields")
            for field in fields:
                if field in metrics:
                    numeric(None if actual[field] == "" else float(actual[field]), computed[field])
                else:
                    demand(actual[field] == str(computed[field]), "Incorrect TSV label")
    text = Path(outputs["table"]["path"]).read_text(encoding="ascii")
    expected_lines = []
    for records, key, delta in ((rows, "cell", False), (differences, "contrast", True)):
        for row in records:
            values = []
            for m in METRICS:
                value = row[m + "_pp"] if delta else row[m]
                values.append("NA" if value is None else f"{float(value):+.3f}" if delta else f"{100 * float(value):.3f}")
            expected_lines.append("| " + " | ".join([row["stratum"], str(row["families"]), row[key], *values]) + " |")
    actual_lines = [line for line in text.splitlines() if line.startswith("| ")
                    and line.split(" | ")[0][2:] in BINS]
    demand(actual_lines == expected_lines, "Incorrect human numeric rows")
    demand("no new intervals" in text and "Initial HMM search on, downstream profiles off" in text
           and "Failed R1 timing stays ineligible" in text and "biological absence" in text, "Missing table scope")


def verify(stage, path, repo, commit):
    checked = [record(path)]
    report = json.loads(Path(path).read_text())
    scope(report)
    input_order = ("annotations", "domain_report", "domain_reader", "alignments", "protocol") if stage == "features" else (
        "features", "feature_reader", "counts", "counts_reader")
    demand(set(report["inputs"]) == set(input_order), "Unexpected input set")
    documents = {}
    for key, ref in report["inputs"].items():
        demand(record(ref["path"]) == ref and (key not in PINS or ref["sha256"] == PINS[key]), "Changed input binding")
        checked.append(ref)
        # Feature construction/readback bind the historical domain report, never its outcomes.
        if key not in {"domain_report", "protocol"}:
            documents[key] = json.loads(Path(ref["path"]).read_text())
    committed_source(report["source"], repo, report["source_commit"], "prepare_swiss_ordered_pfam.py")
    own_source = record(__file__)
    committed_source(own_source, repo, commit, "readback_swiss_ordered_pfam.py")
    checked += [report["source"], own_source]
    if stage == "features":
        admission = documents["domain_reader"]
        demand(admission["proteins_checked"] == 563 and admission["report"] == report["inputs"]["domain_report"]
               and report["inputs"]["annotations"] in admission["checked_inputs"]
               and documents["annotations"]["prediction_statistics_evaluated"] is False
               and documents["alignments"]["failed_families"] == []
               and documents["alignments"]["prediction_statistics_evaluated"] is False, "Incorrect inherited verification")
        details = feature_readback(report, documents["annotations"], documents["alignments"], checked)
    else:
        feature_reader, count_reader = documents["feature_reader"], documents["counts_reader"]
        demand(feature_reader["report"] == report["inputs"]["features"] and feature_reader["status"] == "ordered_features_verified"
               and feature_reader["genes_checked"] == 563 and feature_reader["families_checked"] == 18
               and count_reader["report"] == report["inputs"]["counts"] and count_reader["family_rows_checked"] == 54,
               "Missing full inherited feature/count readback")
        scope(feature_reader)
        rows, differences = score_readback(report, documents["counts"], documents["features"])
        demand(set(report["outputs"]) == {"scores", "differences", "table"}, "Incorrect output set")
        checked += list(report["outputs"].values())
        for ref in report["outputs"].values():
            demand(record(ref["path"]) == ref, "Changed table binding")
        tables_readback(report["outputs"], rows, differences)
        details = dict(family_rows_checked=len(report["family_rows"]), score_rows_checked=len(rows), differences_checked=len(differences),
                       human_numeric_rows_checked=len(rows) + len(differences), bins=report["bins"])
    demand(report["checked_inputs"] == [*(report["inputs"][k] for k in input_order),
           *([r["alignment"] for r in documents["alignments"]["runs"]] if stage == "features" else []), report["source"]],
           "Incorrect primary direct-input inventory")
    for ref in checked:
        demand(record(ref["path"]) == ref, "Changed input during readback")
    return dict(schema="swiss_ordered_pfam_readback_v1", status="ordered_features_verified" if stage == "features" else "ordered_projection_verified",
                report=checked[0], source=own_source, source_commit=commit, checked_inputs=checked,
                shared_parsers=["json", "Bio.SeqIO" if stage == "features" else "csv"],
                new_bootstrap_draws=0, **dict.fromkeys(FALSE_FLAGS, False), **details)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("features", "scores"))
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = verify(args.stage, args.report, args.repo.resolve(), args.source_commit)
    with args.output.open("x", encoding="ascii") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: v for k, v in result.items() if k.endswith("_checked") or k == "status"}))
