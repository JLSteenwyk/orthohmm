"""Construct ordered-annotation descriptors before separate native-count projection."""

import argparse
from collections import Counter
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess

FEATURE_PINS = {
    "annotations": ("benchmark_tools/results/swiss_domain_annotation_inventory_20260917.json", "d5e269158c2fb603acc0805a1c140a88758342d7b8e32b3094750ade22fbc06c"),
    "domain_report": ("benchmark_tools/results/native_qfo_swiss_domain_strata_20261006_v1/report.json", "d3235c712e711f005b3182a3621c0079f8a85c3c7abd5db016d5e8ea95884cd6"),
    "domain_reader": ("benchmark_tools/results/native_qfo_swiss_domain_strata_readback_20261006_v1.json", "9e77b67ad97e6824e9442691be7a05ba1552d62ffe696ced6034a996ede1945d"),
    "alignments": ("benchmarks/results/swiss_model_divergence_20261007_v1/report.json", "ceb4ca3675836e1815fcc1f14c29c8a50947633a8a5e8dd6c67477bdd853065c"),
    "protocol": ("benchmark_tools/results/SWISS_ORDERED_PFAM_PROTOCOL_20261007.md", "c58b217142bf00765cb8733819c4c93d06299eb38eb2f318d384b3fcd9395272"),
}
COUNT_PINS = {
    "counts": ("benchmark_tools/results/native_qfo_three_cell_strata_20261007_v1/report.json", "55ed72ce771c756424877cfa412da329ff104efb7dca4b35113b0d80703a59e5"),
    "counts_reader": ("benchmark_tools/results/native_qfo_three_cell_strata_readback_20261007_v2.json", "6f3eef63c0e96e1d85805fcfb70ed9b66fbb1db76a08ca2a105483422a0f6969"),
}
BINS = ("all", "all_members_usable_same_signature", "all_members_usable_multiple_signatures", "some_members_unusable")
CELLS = ("p0_c0_r0", "p0_c0_r1", "p0_c1_r0")
METRICS = ("F1", "PPV", "TPR")
SCOPE = dict(new_uncertainty=False, independent_confirmation=False, new_accuracy_or_resource_admission=False,
             publication_ready=False, scientific_timings_admitted=False, raw_scorer_repeated=False,
             alignment_repeated=False, new_bootstrap_draws=0)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def record(path):
    path = Path(path).resolve(strict=True)
    return dict(path=str(path), bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def check(ref):
    require(record(ref["path"]) == ref, "Changed direct binding")


def save(path, value):
    with Path(path).open("x", encoding="ascii") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def pinned_inputs(repo, pins, parse):
    refs, docs = {}, {}
    for key, (name, sha) in pins.items():
        refs[key] = record(repo / name)
        require(refs[key]["sha256"] == sha, "Changed protocol input: " + key)
        if key in parse:
            docs[key] = json.loads(Path(refs[key]["path"]).read_text())
    return refs, docs


def source(repo, commit):
    ref = record(__file__)
    name = Path(__file__).resolve().relative_to(repo).as_posix()
    require(Path(__file__).read_bytes() == subprocess.check_output(["git", "show", commit + ":" + name], cwd=repo),
            "Source was not committed before execution")
    return ref


def alignment_lengths(report):
    from Bio import SeqIO
    lengths, refs = {}, []
    require([r["family"] for r in report["runs"]] == sorted(report["memberships"]), "Incomplete alignment families")
    for run in report["runs"]:
        check(run["alignment"])
        refs.append(run["alignment"])
        sequences = list(SeqIO.parse(run["alignment"]["path"], "fasta"))
        require(sorted(r.id for r in sequences) == report["memberships"][run["family"]]
                and all(len(r.seq) == run["columns"] for r in sequences), "Changed alignment population")
        for seq in sequences:
            require(seq.id not in lengths, "Duplicate canonical protein")
            lengths[seq.id] = sum(c != "-" for c in str(seq.seq))
    return lengths, refs


def gene_feature(annotation, length):
    require(type(length) is int and length > 0 and type(annotation["length"]) is int
            and annotation["length"] == length, "Annotation/sequence length mismatch")
    instances = annotation["ordered_pfam_instances"]
    require(isinstance(instances, list), "Invalid instance list")
    tuples = []
    for row in instances:
        require(set(row) == {"start", "end", "domain"} and type(row["start"]) is int
                and type(row["end"]) is int and 0 <= row["start"] <= row["end"] <= length
                and isinstance(row["domain"], str) and row["domain"].startswith("pfam_"), "Invalid Pfam coordinates/type")
        tuples.append((row["start"], row["end"], row["domain"]))
    require(len(tuples) == len(set(tuples)), "Duplicate domain instance")
    tuples.sort()
    frequency = Counter(t[2] for t in tuples)
    require(annotation["pfam_types"] == sorted(frequency)
            and type(annotation["pfam_type_count"]) is int and type(annotation["pfam_instance_count"]) is int
            and annotation["pfam_type_count"] == len(frequency)
            and annotation["pfam_instance_count"] == len(tuples)
            and annotation["has_repeated_pfam_type"] is (len(tuples) > len(frequency)), "Inconsistent retained Pfam counts")
    ambiguous = any(a == b for a, b, _ in tuples) or any(b[0] <= a[1] for a, b in zip(tuples, tuples[1:]))
    status = "zero_pfam" if not tuples else "order_ambiguous" if ambiguous else "usable"
    return dict(length=length, status=status, ordered_signature=[t[2] for t in tuples] if status == "usable" else None,
                ordered_instances=[dict(start=a, end=b, domain=d) for a, b, d in tuples],
                type_multiset=[[d, frequency[d]] for d in sorted(frequency)])


def family_feature(members, genes):
    require(members == sorted(set(members)) and bool(members), "Invalid canonical member list")
    usable = [genes[g] for g in members if genes[g]["status"] == "usable"]
    signatures = {tuple(g["ordered_signature"]) for g in usable}
    comparable = discordant = 0
    for a, b in itertools.combinations(usable, 2):
        if a["type_multiset"] == b["type_multiset"]:
            comparable += 1
            discordant += a["ordered_signature"] != b["ordered_signature"]
    category = BINS[3] if len(usable) < len(members) else BINS[1] if len(signatures) == 1 else BINS[2]
    return dict(members=members, member_count=len(members), usable_count=len(usable),
                ambiguous_count=sum(genes[g]["status"] == "order_ambiguous" for g in members),
                zero_pfam_count=sum(genes[g]["status"] == "zero_pfam" for g in members),
                distinct_usable_signatures=len(signatures), same_multiset_comparable_pairs=comparable,
                same_multiset_order_discordant_pairs=discordant, stratum=category)


def construct(repo, output, commit):
    refs, docs = pinned_inputs(repo, FEATURE_PINS, {"annotations", "domain_reader", "alignments"})
    src = source(repo, commit)
    inventory, admission, alignments = (docs[k] for k in ("annotations", "domain_reader", "alignments"))
    require(inventory["prediction_statistics_evaluated"] is False and admission["proteins_checked"] == 563
            and admission["report"] == refs["domain_report"] and refs["annotations"] in admission["checked_inputs"]
            and alignments["failed_families"] == [] and alignments["prediction_statistics_evaluated"] is False,
            "Missing inherited annotation/alignment verification")
    memberships = alignments["memberships"]
    expected = [g for f in sorted(memberships) for g in memberships[f]]
    require(len(memberships) == 18 and len(expected) == len(set(expected)) == 563
            and set(expected) == set(inventory["genes"]), "Changed canonical annotation universe")
    lengths, align_refs = alignment_lengths(alignments)
    genes = {g: gene_feature(inventory["genes"][g], lengths[g]) for g in sorted(expected)}
    families = {f: family_feature(memberships[f], genes) for f in sorted(memberships)}
    bins = {name: sorted(f for f in families if name == "all" or families[f]["stratum"] == name) for name in BINS}
    checked = [*refs.values(), *align_refs, src]
    for ref in checked:
        check(ref)
    result = dict(schema="swiss_ordered_pfam_features_v1", status="features_constructed_unverified",
                  inputs=refs, source=src, source_commit=commit, checked_inputs=checked,
                  memberships=memberships, genes=genes, families=families, bins=bins,
                  prediction_statistics_evaluated=False, original_annotation_extraction_repeated=False, **SCOPE)
    save(output, result)
    return result


def macro_counts(rows):
    if not rows:
        return dict.fromkeys(METRICS)
    points = []
    for row in rows:
        c = row["counts_without_prior"]
        require(set(c) == {"TP", "FP", "FN", "TN"} and all(type(v) is int and v >= 0 for v in c.values())
                and sum(c.values()) > 0, "Invalid counts")
        tp, fp, fn = (c[k] / 2 + 1 for k in ("TP", "FP", "FN"))
        points.append((tp / (tp + fp), tp / (tp + fn)))
    p, r = (math.fsum(v[i] for v in points) / len(points) for i in (0, 1))
    return dict(F1=2 * p * r / (p + r), PPV=p, TPR=r)


def projected_rows(counts, bins):
    families = sorted(counts["memberships"])
    require(set(bins) == set(BINS) and bins["all"] == families
            and all(v == sorted(set(v)) for v in bins.values())
            and sorted(f for k in BINS[1:] for f in bins[k]) == families, "Invalid fixed-bin partition")
    require([(r["cell"], r["family"]) for r in counts["family_rows"]] ==
            [(c, f) for c in CELLS for f in families], "Incomplete or duplicate family count rows")
    index = {(r["cell"], r["family"]): r for r in counts["family_rows"]}
    rows, differences = [], []
    for name in BINS:
        members = bins[name]
        values = {cell: macro_counts([index[cell, f] for f in members]) for cell in CELLS}
        for cell in CELLS:
            rows.append(dict(stratum=name, cell=cell, families=len(members), family_members=members,
                status="descriptive" if members else "empty",
                prediction_semantics="native_pair" if cell == CELLS[1] else "group_clique", **values[cell]))
        for contrast, candidate in (("R_at_P0_C0", CELLS[1]), ("C_at_P0_R0", CELLS[2])):
            differences.append(dict(stratum=name, contrast=contrast, candidate=candidate, reference=CELLS[0],
                families=len(members), family_members=members, status="descriptive" if members else "empty",
                **{m + "_pp": 100 * (values[candidate][m] - values[CELLS[0]][m]) if members else None for m in METRICS}))
    return rows, differences


def write_tables(output, rows, differences):
    refs = {}
    for name, records, fields in (
        ("scores", rows, ("stratum", "cell", "families", "status", *METRICS, "prediction_semantics")),
        ("differences", differences, ("stratum", "contrast", "families", "status", *(m + "_pp" for m in METRICS)))):
        path = output / (name + ".tsv")
        with path.open("x", newline="", encoding="ascii") as stream:
            writer = csv.DictWriter(stream, delimiter="\t", fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(records)
        refs[name] = record(path)
    lines = ["# Ordered-Pfam Annotation Strata", "", "Descriptive development-exposed points; no new intervals.", "",
        "| Stratum | Families | Cell | F1 (%) | Precision (%) | Recall (%) |", "| --- | ---: | --- | ---: | ---: | ---: |"]
    for row in rows:
        values = ["NA" if row[m] is None else f"{100 * row[m]:.3f}" for m in METRICS]
        lines.append("| " + " | ".join([row["stratum"], str(row["families"]), row["cell"], *values]) + " |")
    lines += ["", "| Stratum | Families | Contrast | Delta F1 (pp) | Delta Precision (pp) | Delta Recall (pp) |",
              "| --- | ---: | --- | ---: | ---: | ---: |"]
    for row in differences:
        values = ["NA" if row[m + "_pp"] is None else f"{row[m + '_pp']:+.3f}" for m in METRICS]
        lines.append("| " + " | ".join([row["stratum"], str(row["families"]), row["contrast"], *values]) + " |")
    lines += ["", "Initial HMM search on, downstream profiles off; conditional C/R effects, not an interaction.",
              "Predicted domain order and unusable annotations do not establish complete architecture or biological absence.",
              "Failed R1 timing stays ineligible; no raw scoring, bootstrap, inference, default or independent-validation claim."]
    path = output / "TABLE.md"
    with path.open("x", encoding="ascii") as stream:
        stream.write("\n".join(lines) + "\n")
    refs["table"] = record(path)
    return refs


def project(repo, features_path, reader_path, output, commit):
    feature_ref, reader_ref = record(features_path), record(reader_path)
    features, reader = (json.loads(p.read_text()) for p in (features_path, reader_path))
    require(features["schema"] == "swiss_ordered_pfam_features_v1"
            and features["status"] == "features_constructed_unverified"
            and features["prediction_statistics_evaluated"] is False
            and reader["status"] == "ordered_features_verified" and reader["report"] == feature_ref
            and reader["genes_checked"] == 563 and reader["families_checked"] == 18, "Unverified features")
    refs, docs = pinned_inputs(repo, COUNT_PINS, set(COUNT_PINS))
    src = source(repo, commit)
    counts, counts_reader = docs["counts"], docs["counts_reader"]
    require(counts["memberships"] == features["memberships"] and counts_reader["report"] == refs["counts"]
            and counts_reader["family_rows_checked"] == 54 and len(counts["family_rows"]) == 54
            and [(r["cell"], r["family"]) for r in counts["family_rows"]] ==
                [(c, f) for c in CELLS for f in sorted(features["memberships"])]
            and counts["cells"][1]["timing_eligible"] is False and counts["cells"][1]["timing_admitted"] is False,
            "Changed native universe or failed timing")
    rows, differences = projected_rows(counts, features["bins"])
    output.mkdir(parents=True, exist_ok=False)
    tables = write_tables(output, rows, differences)
    checked = [feature_ref, reader_ref, *refs.values(), src]
    for ref in checked:
        check(ref)
    result = dict(schema="native_swiss_ordered_pfam_projection_v1", status="projection_constructed_unverified",
                  inputs=dict(features=feature_ref, feature_reader=reader_ref, **refs), source=src,
                  source_commit=commit, checked_inputs=checked, memberships=counts["memberships"],
                  bins=features["bins"], cells=counts["cells"], family_rows=counts["family_rows"],
                  rows=rows, differences=differences, outputs=tables, **SCOPE)
    save(output / "report.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("features", "scores"))
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--features", type=Path)
    parser.add_argument("--feature-reader", type=Path)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    if args.stage == "scores" and (args.features is None or args.feature_reader is None):
        parser.error("Scores require complete features and separate feature readback")
    result = construct(args.repo.resolve(), args.output, args.source_commit) if args.stage == "features" else project(
        args.repo.resolve(), args.features, args.feature_reader, args.output, args.source_commit)
    print(json.dumps(dict(status=result["status"], bins={k: len(v) for k, v in result["bins"].items()})))
