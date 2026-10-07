"""Extend verified SwissTrees fixed bins using admitted native family counts."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

CELLS = ("p0_c0_r0", "p0_c0_r1", "p0_c1_r0")
METRICS = ("F1", "PPV", "TPR")
COUNTS = ("TP", "FP", "FN", "TN")
SUITES = ("sequence", "domain", "duplication")
CONTRASTS = (("R_at_P0_C0", CELLS[1]), ("C_at_P0_R0", CELLS[2]))
PINS = {
    "sequence": ("native_qfo_swiss_sequence_strata_20261006_v1/report.json",
                 "4650f614d2fdeebd65cbcbf4999c61c83cc549560d2738bfbb5502f686f660e3"),
    "sequence_reader": ("native_qfo_swiss_sequence_strata_readback_20261006.json",
                        "8d0356bc0e159d204939271578e23a6a5b6efcf6198b75a18c0f3d38e0e67c0a"),
    "domain": ("native_qfo_swiss_domain_strata_20261006_v1/report.json",
               "d3235c712e711f005b3182a3621c0079f8a85c3c7abd5db016d5e8ea95884cd6"),
    "domain_reader": ("native_qfo_swiss_domain_strata_readback_20261006_v1.json",
                      "9e77b67ad97e6824e9442691be7a05ba1552d62ffe696ced6034a996ede1945d"),
    "duplication": ("native_qfo_swiss_duplication_strata_20261006_v1/report.json",
                    "e0ecbc9d26738d0d3d573ca31772085b28d145ba2deaa9383c8cff75bd0d67b5"),
    "duplication_reader": ("native_qfo_swiss_duplication_strata_readback_20261006_v1.json",
                           "84eb3d9fe02b5984cecaa06355d1cb1d68c54aec09ab588eef8b3ada42f1957b"),
    "features": ("corrected_swiss_sequence_strata_20260918.json",
                 "912363276fa5456f11aeff9f398fce71aec8d377c8c8eda7b144105945e554f1"),
    "candidate": ("native_qfo_candidate_swiss_counts_20261006_v1.json",
                  "7ed2ea1a7a94d55b28ecbadd7bd2be8e99df64b9fab3c266ebd987e8dbcfc479"),
    "candidate_reader": ("native_qfo_candidate_swiss_readback_20261006_v1.json",
                         "07941ba9c55ba7dc14c5ba53df37d28cb3b2f60f599b9d177b10c9afa99c1329"),
    "protocol": ("NATIVE_QFO_THREE_CELL_STRATA_PROTOCOL_20261007.md",
                 "9e14a818df91a4ab12d6db5c060e556cf872c4eb7d45dfc43f8331fcbbc685ca"),
}
SCOPE = dict(new_bootstrap_draws=0, new_uncertainty=False,
             new_accuracy_or_resource_admission=False, independent_confirmation=False,
             publication_ready=False, scientific_timings_admitted=False,
             raw_evidence_reparsed=False, original_tree_traversal_repeated=False)
LIMITATIONS = [
    "Development-exposed descriptive fixed bins; no subgroup intervals, significance or tuning.",
    "Macro family precision/recall then harmonic F1; not pooled pairs or mean-family F1.",
    "Initial HMM search on and downstream profile refinement off in all three cells.",
    "Two conditional contrasts against baseline; candidate/reconciliation interaction not identified.",
    "Bins overlap and all-family rows repeat; they are not independent findings.",
    "Length/composition proxies are not calibrated divergence or literal fragment truth.",
    "No explicit fragment annotation is not evidence of completeness.",
    "Pfam descriptors are not complete validated domain architectures or domain-loss truth.",
    "Reference informative-node duplication fractions are not ancestral duplication histories; default S is not explicit speciation.",
    "Direct artifact/source checks only; prior raw validation inherited, not repeated or transitively readmitted.",
    "Cell 7 failed timing remains ineligible; no timing repair or new inference resources.",
    "Shared-host contention effects are unknown and potentially tool-dependent, not isolated speed evidence.",
]
SCORE_FIELDS = ("suite", "cell", "stratum", "families", "status", *METRICS, "prediction_semantics")
DIFF_FIELDS = ("suite", "contrast", "candidate", "reference", "stratum", "families", "status", *METRICS)
FAMILY_FIELDS = ("cell", "family", *COUNTS, *METRICS)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def record(path):
    path = Path(path).resolve(strict=True)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest.hexdigest())


def same(actual, expected, message):
    require(set(actual) == set(expected) and all(
        (actual[m] is None if expected[m] is None else
         type(actual[m]) in (int, float) and math.isfinite(actual[m])
         and math.isclose(actual[m], expected[m], rel_tol=0, abs_tol=1e-12))
        for m in expected), message)


def statistics(counts):
    require(set(counts) == set(COUNTS) and all(type(v) is int and v >= 0 for v in counts.values())
            and sum(counts.values()) > 0, "Invalid integer family counts")
    tp, fp, fn = (counts[k] / 2 + 1 for k in COUNTS[:3])
    p, r = tp / (tp + fp), tp / (tp + fn)
    return dict(F1=2 * p * r / (p + r), PPV=p, TPR=r)


def aggregate(values):
    if not values:
        return dict.fromkeys(METRICS)
    p, r = (sum(v[m] for v in values) / len(values) for m in ("PPV", "TPR"))
    return dict(F1=2 * p * r / (p + r), PPV=p, TPR=r)


def load_inputs(repo):
    root = Path(repo).resolve() / "benchmark_tools/results"
    documents, refs, checked = {}, {}, []
    for key, (name, sha) in PINS.items():
        ref = record(root / name)
        require(ref["sha256"] == sha, "Changed frozen input: " + key)
        refs[key] = ref
        checked.append(ref)
        if name.endswith(".json"):
            documents[key] = json.loads(Path(ref["path"]).read_text())
    for key in (*SUITES, *(s + "_reader" for s in SUITES), "features", "candidate", "candidate_reader"):
        source = documents[key]["source"]
        require(record(source["path"]) == source, "Changed inherited source: " + key)
        if source not in checked:
            checked.append(source)
    return documents, refs, checked


def project(docs):
    memberships = docs["features"]["family_memberships"]
    families = sorted(memberships)
    genes = [g for f in families for g in memberships[f]]
    require(len(families) == 18 and len(genes) == len(set(genes)) == 563
            and all(memberships[f] == sorted(set(memberships[f])) for f in families),
            "Changed canonical family/member universe")
    require(docs["features"]["prediction_statistics_evaluated"] is False,
            "Feature inventory used predictions")
    require(docs["domain"]["memberships"] == docs["duplication"]["memberships"] == memberships,
            "Prior canonical memberships disagree")
    prior = docs["sequence"]["family_rows"]
    require([(r["cell"], r["family"]) for r in prior] ==
            [(c, f) for c in CELLS[:2] for f in families], "Changed prior family rows")
    for suite in SUITES:
        require(docs[suite]["schema"] == "native_qfo_swiss_" + suite + "_strata_v1"
                and docs[suite]["family_rows"] == prior, "Prior family counts disagree")
        require(all(docs[suite][k] is False for k in (
            "new_uncertainty", "new_accuracy_or_resource_admission", "independent_confirmation", "publication_ready")),
            "Inflated prior scope")
        reader = docs[suite + "_reader"]
        require(reader["family_rows_checked"] == 36 and reader["raw_rows_checked"] == 21530,
                "Incomplete inherited readback")
    audit, reader = docs["candidate"], docs["candidate_reader"]
    require(audit["schema"] == "native_qfo_swiss_family_count_audit_v1"
            and audit["status"] == "supplied_native_swiss_family_counts_verified"
            and audit["families"] == families and audit["reference_relation_count"] == 10765
            and [c["cell"] for c in audit["cells"]] == [CELLS[2]]
            and type(audit["new_bootstrap_draws"]) is int and audit["new_bootstrap_draws"] == 0
            and all(audit[k] is False for k in ("historical_intervals_attached",
                "new_accuracy_or_resource_admission", "independent_confirmation", "publication_ready")),
            "Changed candidate admission scope")
    require(reader["cells"] == [CELLS[0], CELLS[2]] and reader["families_checked"] == 18
            and reader["native_family_records_checked"] == 36 and reader["candidate_pair_labels_matched"] == 10765,
            "Incomplete candidate count readback")
    candidate = audit["cells"][0]
    require(candidate["index"] == 8 and candidate["native_job_id"] == 22444
            and candidate["retained_family_records_identical"] is True
            and candidate["retained_aggregate_identical"] is True
            and [r["family"] for r in candidate["families"]] == families,
            "Wrong candidate native identity")
    family_rows, values = [], {}
    for row in prior:
        value = statistics(row["counts_without_prior"])
        same({m: row[m] for m in METRICS}, value, "Prior family statistics differ")
        family_rows.append(dict(row))
        values[row["cell"], row["family"]] = value
    for row in candidate["families"]:
        f = row["family"]
        require(row["represented_genes"] == memberships[f], "Candidate membership mismatch")
        value = statistics(row["counts_without_prior"])
        same(row["statistics_with_prior"], value, "Candidate family statistics differ")
        family_rows.append(dict(cell=CELLS[2], family=f, counts_without_prior=row["counts_without_prior"], **value))
        values[CELLS[2], f] = value
    count_index = {(r["cell"], r["family"]): r["counts_without_prior"] for r in family_rows}
    for f in families:
        totals = {(sum(count_index[c, f].values()), count_index[c, f]["TP"] + count_index[c, f]["FN"])
                  for c in CELLS}
        require(len(totals) == 1, "Changed family classified/reference-positive universe")
    require(sum(sum(count_index[CELLS[0], f].values()) for f in families) == 10765,
            "Changed reference relation count")
    same(candidate["aggregate"], aggregate([values[CELLS[2], f] for f in families]),
         "Candidate aggregate differs")
    for c in (CELLS[0], CELLS[2]):
        same(reader["rational_macro_points"][c], aggregate([values[c, f] for f in families]),
             "Prior rational macro readback differs")
    bins, rows, differences = {}, [], []
    for suite, size in zip(SUITES, (11, 5, 4)):
        old_rows = docs[suite]["rows"]
        selected = [r for r in old_rows if r["cell"] == CELLS[0]]
        names = [r["stratum"] for r in selected]
        require(len(names) == len(set(names)) == size and "all" in names,
                "Missing/duplicate prior bins")
        bins[suite] = {r["stratum"]: r["family_members"] for r in selected}
        require(bins[suite]["all"] == families and all(
            m == sorted(set(m)) and set(m) <= set(families) for m in bins[suite].values()),
            "Invalid bin membership")
        expected = (dict(all=families, **docs["features"]["primary_strata"], **docs["features"]["secondary_strata"])
                    if suite == "sequence" else docs[suite]["bins"])
        require(bins[suite] == expected, "Frozen bin membership changed")
        require([(r["cell"], r["stratum"]) for r in old_rows] ==
                [(c, n) for c in CELLS[:2] for n in names], "Wrong prior bin row inventory")
        lookup = {(r["cell"], r["stratum"]): r for r in old_rows}
        for c in CELLS:
            for name, members in bins[suite].items():
                value = aggregate([values[c, f] for f in members])
                if c in CELLS[:2]:
                    old = lookup[c, name]
                    require(old["family_members"] == members and old["families"] == len(members)
                            and old["prediction_semantics"] == ("resolved_native_pairs" if c == CELLS[1] else "group_clique"),
                            "Prior bin inventory/semantics differ")
                    same({m: old[m] for m in METRICS}, value, "Prior bin score does not reproduce")
                rows.append(dict(suite=suite, cell=c, stratum=name, families=len(members), family_members=members,
                                 status="descriptive" if members else "empty_bin",
                                 prediction_semantics="resolved_native_pairs" if c == CELLS[1] else "group_clique", **value))
        by_key = {(r["cell"], r["stratum"]): r for r in rows if r["suite"] == suite}
        for contrast, c in CONTRASTS:
            for name, members in bins[suite].items():
                differences.append(dict(suite=suite, contrast=contrast, candidate=c, reference=CELLS[0], stratum=name,
                    families=len(members), family_members=members, status="descriptive" if members else "empty_bin",
                    **{m: by_key[c, name][m] - by_key[CELLS[0], name][m] if members else None for m in METRICS}))
    cell_sources = [dict(c) for c in docs["domain"]["cells"]]
    require([c["cell"] for c in cell_sources] == list(CELLS[:2])
            and cell_sources[1]["timing_admitted"] is False and cell_sources[1]["timing_eligible"] is False,
            "Failed timing relabeled")
    cell_sources.append({k: candidate[k] for k in ("cell", "index", "native_job_id", "raw_file", "admission")})
    return dict(memberships=memberships, cells=cell_sources, bins=bins, family_rows=family_rows,
                rows=rows, differences=differences)


def table(result):
    lookup = {(r["suite"], r["cell"], r["stratum"]): r for r in result["rows"]}
    changes = {(r["suite"], r["contrast"], r["stratum"]): r for r in result["differences"]}
    lines = ["# Native SwissTrees Fixed Strata: Three Cells", "",
             "F1 in percent; changes in percentage points against p0_c0_r0. NA denotes an empty bin.",
             "R: reconciliation at P0/C0. C: candidate expansion at P0/R0. Initial HMM search on in all cells.",
             "These are development-exposed descriptions, not intervals or evidence of a causal interaction."]
    for suite in SUITES:
        lines.extend(["", "## " + suite.title(), "",
            "| Bin | Families | Baseline F1 | R1 F1 | C1 F1 | R F1 change | R precision change | R recall change | C F1 change | C precision change | C recall change |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"])
        for name in sorted(result["bins"][suite]):
            members = result["bins"][suite][name]
            numbers = [lookup[suite, c, name]["F1"] for c in CELLS]
            numbers.extend(changes[suite, contrast, name][m] for contrast, _ in CONTRASTS for m in METRICS)
            texts = ["NA" if v is None else format(100 * v, ".3f" if i < 3 else "+.3f")
                     for i, v in enumerate(numbers)]
            lines.append("| " + " | ".join([name, str(len(members)), *texts]) + " |")
    return "\n".join([*lines, "", *["- " + s for s in LIMITATIONS], ""])


def write_tsv(path, rows, fields):
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({k: "NA" if r[k] is None else r[k] for k in fields} for r in rows)


def export(repo, output):
    output = Path(output)
    require(not output.exists() and not output.is_symlink(), "Output already exists")
    docs, refs, checked = load_inputs(repo)
    for suite in SUITES:
        require(docs[suite + "_reader"]["report"] == refs[suite], "Prior readback bound to wrong report")
    require(docs["candidate_reader"]["audit"] == refs["candidate"]
            and docs["sequence"]["strata"] == refs["features"], "Changed candidate/feature readback binding")
    result = dict(schema="native_qfo_three_cell_strata_v1", source=record(__file__),
                  inputs=refs, checked_inputs=checked, limitations=LIMITATIONS, **SCOPE, **project(docs))
    for ref in checked:
        require(record(ref["path"]) == ref, "Direct evidence changed before output")
    output.mkdir(parents=True)
    write_tsv(output / "scores.tsv", result["rows"], SCORE_FIELDS)
    write_tsv(output / "differences.tsv", result["differences"], DIFF_FIELDS)
    flattened = [dict(cell=r["cell"], family=r["family"], **r["counts_without_prior"],
                      **{m: r[m] for m in METRICS}) for r in result["family_rows"]]
    write_tsv(output / "family_counts.tsv", flattened, FAMILY_FIELDS)
    with (output / "TABLE.md").open("x") as stream:
        stream.write(table(result))
    result["outputs"] = [record(output / n) for n in ("scores.tsv", "differences.tsv", "family_counts.tsv", "TABLE.md")]
    for ref in checked:
        require(record(ref["path"]) == ref, "Direct evidence changed during export")
    with (output / "report.json").open("x") as stream:
        json.dump(result, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = export(args.repo, args.output)
    print(json.dumps({k: len(result[k]) for k in ("family_rows", "rows", "differences")}))
