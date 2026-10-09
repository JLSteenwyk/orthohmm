"""Project admitted profile-on counts into unchanged, previously checked fixed bins."""

import argparse
import copy
import json
from pathlib import Path

from benchmark_tools.export_native_qfo_three_cell_strata import (
    METRICS, COUNTS, statistics, aggregate, same, write_tsv, record, require,
)


CELLS = ("p0_c0_r0", "p0_c0_r1", "p0_c1_r0", "p1_c0_r1")
SUITES = ("sequence", "domain", "duplication", "model_distance")
CONTRASTS = (("R_at_P0_C0", CELLS[1], CELLS[0]),
             ("C_at_P0_R0", CELLS[2], CELLS[0]),
             ("P_at_C0_R1", CELLS[3], CELLS[1]))
PINS = {
    "fixed": ("native_qfo_three_cell_strata_20261007_v1/report.json", "55ed72ce771c756424877cfa412da329ff104efb7dca4b35113b0d80703a59e5"),
    "fixed_reader": ("native_qfo_three_cell_strata_readback_20261007_v2.json", "6f3eef63c0e96e1d85805fcfb70ed9b66fbb1db76a08ca2a105483422a0f6969"),
    "distance": ("swiss_model_divergence_strata_20261007_v1/report.json", "ddfa27757c4c256df9c881079ac352c37d7e17646f75fbd417b624e7980be78f"),
    "distance_reader": ("swiss_model_divergence_strata_readback_20261007_v1.json", "186c0151137fd825513ef36c27df595b1a0f1cd4b2b4bf68de4b9665ad328d0e"),
    "profile": ("native10_allocated_swiss_counts_20261007_v1.json", "d5bafed830bf93e46b810728df14ffa34185d58ffe3eac78a693b0159d3d0a40"),
    "profile_reader": ("native_qfo_profile_swiss_readback_20261007_v1.json", "3f62b487014cb9ba99ea36b5adc15535f62c7e39c9a99f3cb4eb8c537523005d"),
    "snapshot": ("native_qfo_terminal_failures_20261009_v1/report.json", "863ebc334dbdb5a41d1a9e9bba56be82896f0ff26c98366047928c7ef802966f"),
}
PROTOCOL = "NATIVE_QFO_FOUR_CELL_STRATA_PROTOCOL_20261009.md"
PROTOCOL_SHA = "fa512f814404bca3d88436b2484b32bda03a64d4efebe039a7e551f2791de58c"
KERNEL_SHA = "e6d0980fc71636e4cf0ae0ba4decf13caff54ad502cd182d36b6db773fab6586"
SCOPE = dict(new_bootstrap_draws=0, new_uncertainty=False, independent_confirmation=False,
    new_accuracy_or_resource_admission=False, scientific_timings_admitted=False,
    raw_evidence_reparsed=False, original_tree_traversal_repeated=False, publication_ready=False)
LIMITATIONS = [
    "Retrospective development-exposed fixed bins, not independent confirmation or causal effects.",
    "Initial HMM search stays on; P is downstream refinement, and R also changes pair semantics.",
    "Three conditional contrasts only, not a complete factorial or candidate/reconciliation interaction.",
    "Bins overlap and all-family rows repeat; neither is an independent finding.",
    "No subgroup intervals, significance, bootstrap, tuning, cutoffs or default promotion.",
    "Length/fragment flags do not establish fragment truth or completeness.",
    "Pfam descriptors are incomplete architecture annotations; reference duplication fractions are not ancestral histories.",
    "Estimated WAG+G4 distances include model/alignment/sampling/paralogy effects, not known genealogy or time.",
    "Direct retained-record/source checks, not repeated raw scoring or transitive scientific admission.",
    "Recovered cell timing stays failed/ineligible; shared-host timing distortion is unknown and tool-dependent.",
]
SCORE_FIELDS = ("suite", "cell", "stratum", "families", "status", *METRICS, "prediction_semantics")
DIFF_FIELDS = ("suite", "contrast", "candidate", "reference", "stratum", "families", "status", *METRICS)
FAMILY_FIELDS = ("cell", "family", *COUNTS, *METRICS)


def prepare(repo):
    root = Path(repo).resolve() / "benchmark_tools/results"
    docs, refs, checked = {}, {}, []
    for key, (name, sha) in PINS.items():
        ref = record(root / name)
        require(ref["sha256"] == sha, "Changed four-cell input: " + key)
        docs[key], refs[key] = json.loads(Path(ref["path"]).read_text()), ref
        checked.append(ref)
    protocol = record(root / PROTOCOL)
    require(protocol["sha256"] == PROTOCOL_SHA, "Changed four-cell protocol")
    kernel = record(Path(__file__).with_name("export_native_qfo_three_cell_strata.py"))
    require(kernel["sha256"] == KERNEL_SHA, "Changed statistical kernel")
    checked.extend((protocol, kernel))
    for key, doc in docs.items():
        source = doc["source"]
        require(record(source["path"]) == source, "Changed inherited source: " + key)
        if source not in checked:
            checked.append(source)
    require(docs["fixed_reader"]["report"] == refs["fixed"]
            and docs["distance_reader"]["report"] == refs["distance"]
            and docs["profile_reader"]["audit"] == refs["profile"], "Mixed reader bindings")
    family_values(docs)
    return docs, refs, checked, protocol


def family_values(docs):
    schemas = dict(fixed="native_qfo_three_cell_strata_v1",
        fixed_reader="native_qfo_three_cell_strata_rational_readback_v2",
        distance="native_qfo_swiss_model_divergence_strata_v1",
        distance_reader="swiss_model_divergence_strata_rational_readback_v1",
        profile="allocated_native_qfo_swiss_family_count_audit_v1",
        profile_reader="allocated_native_qfo_profile_swiss_rational_readback_v1",
        snapshot="native_qfo_terminal_failure_reporting_v1")
    for key, schema in schemas.items():
        require(docs[key]["schema"] == schema and docs[key]["publication_ready"] is False,
                "Changed document scope: " + key)
        if key != "snapshot":
            require(type(docs[key]["new_bootstrap_draws"]) is int and docs[key]["new_bootstrap_draws"] == 0
                    and docs[key]["independent_confirmation"] is False
                    and docs[key]["new_accuracy_or_resource_admission"] is False, "Inflated inherited scope")
    for key, scores, differences in (("fixed_reader", 60, 40), ("distance_reader", 9, 6)):
        reader = docs[key]
        require((reader["families_checked"], reader["proteins_checked"], reader["family_rows_checked"],
                 reader["score_rows_checked"], reader["differences_checked"]) == (18, 563, 54, scores, differences),
                "Incomplete prior reader")
    require(docs["profile_reader"]["cells"] == [CELLS[1], CELLS[3]]
            and docs["profile_reader"]["native_family_records_checked"] == 36
            and docs["profile_reader"]["profile_pair_labels_matched"] == 10765,
            "Incomplete profile reader")
    fixed, distance, audit = (docs[k] for k in ("fixed", "distance", "profile"))
    members = fixed["memberships"]
    families = sorted(members)
    genes = [g for family in families for g in members[family]]
    require(len(families) == 18 and len(genes) == len(set(genes)) == 563
            and all(members[f] == sorted(set(members[f])) for f in families)
            and distance["memberships"] == members and audit["families"] == families,
            "Changed canonical family/member universe")
    require(fixed["family_rows"] == distance["family_rows"]
            and [(r["cell"], r["family"]) for r in fixed["family_rows"]] ==
                [(c, f) for c in CELLS[:3] for f in families]
            and fixed["cells"] == distance["cells"] and [r["cell"] for r in fixed["cells"]] == list(CELLS[:3]),
            "Changed inherited family/identity inventory")
    require(audit["selected_index"] == 10 and audit["reference_relation_count"] == 10765
            and audit["status"] == "supplied_allocated_native_swiss_family_counts_verified"
            and [r["cell"] for r in audit["cells"]] == [CELLS[3]], "Wrong profile count admission")
    profile = audit["cells"][0]
    require(profile["index"] == 10 and profile["native_job_id"] == 23902
            and [r["family"] for r in profile["families"]] == families
            and profile["raw_file"] in audit["checked_inputs"], "Wrong profile native identity/count binding")
    snapshot = docs["snapshot"]["rows"]
    require([r["index"] for r in snapshot] == list(range(6, 13))
            and {r["cell"] for r in snapshot if r["accuracy_admitted"]} == set(CELLS), "Changed native snapshot")
    sources = copy.deepcopy(fixed["cells"])
    sources.append({k: profile[k] for k in ("cell", "index", "native_job_id", "admission", "raw_file")})
    sources[-1]["scientific_timings_admitted"] = False
    for source, job in zip(sources, (22435, 22437, 22444, 23902)):
        target = next(r for r in snapshot if r["cell"] == source["cell"])
        require(source["native_job_id"] == target["native_job_id"] == job
                and source["admission"] == target["admission"] and target["accuracy_admitted"] is True,
                "Native admission differs")
    recovered = next(r for r in snapshot if r["cell"] == CELLS[1])
    require(sources[1]["timing_eligible"] is sources[1]["timing_admitted"] is False
            and recovered["timing_eligible"] is recovered["timing_admitted"] is False
            and recovered["resources"] is None, "Recovered timing was repaired")
    require(next(r for r in snapshot if r["cell"] == CELLS[3])["scientific_timings_admitted"] is False,
            "Profile scientific timing was relabeled")
    rows = copy.deepcopy(fixed["family_rows"])
    for row in profile["families"]:
        require(row["represented_genes"] == members[row["family"]], "Profile member set differs")
        value = statistics(row["counts_without_prior"])
        same(row["statistics_with_prior"], value, "Profile family statistic differs")
        rows.append(dict(cell=CELLS[3], family=row["family"], counts_without_prior=copy.deepcopy(row["counts_without_prior"]), **value))
    values, raw = {}, {}
    for row in rows:
        value = statistics(row["counts_without_prior"])
        same({m: row[m] for m in METRICS}, value, "Stored family statistic differs")
        values[row["cell"], row["family"]] = value
        raw[row["cell"], row["family"]] = row["counts_without_prior"]
    for family in families:
        truths = {(raw[c, family]["TP"] + raw[c, family]["FN"],
                   raw[c, family]["FP"] + raw[c, family]["TN"]) for c in CELLS}
        require(len(truths) == 1, "Changed positive/negative reference universe")
    require(sum(sum(raw[CELLS[0], f].values()) for f in families) == 10765, "Changed relation count")
    for cell in CELLS:
        point = aggregate([values[cell, f] for f in families])
        endpoint = next(r for r in snapshot if r["cell"] == cell)["scores"]["SwissTrees"]
        require(abs(point["F1"] - endpoint) <= 5e-8, "Macro endpoint differs from admitted native result")
        if cell in (CELLS[1], CELLS[3]):
            same(docs["profile_reader"]["rational_macro_points"][cell], point, "Profile rational macro differs")
    same(profile["aggregate"], aggregate([values[CELLS[3], f] for f in families]), "Profile aggregate differs")
    return members, rows, values, sources


def project(docs):
    members, family_rows, values, sources = family_values(docs)
    families = sorted(members)
    bins = copy.deepcopy(docs["fixed"]["bins"])
    require(set(bins) == set(SUITES[:3]), "Wrong fixed suites")
    bins["model_distance"] = copy.deepcopy(docs["distance"]["bins"])
    for suite, size in zip(SUITES, (11, 5, 4, 3)):
        require(len(bins[suite]) == size and bins[suite]["all"] == families
                and all(v == sorted(set(v)) and set(v) <= set(families) for v in bins[suite].values()),
                "Changed fixed bin memberships")
    rows, differences = [], []
    for suite in SUITES:
        for cell in CELLS:
            for name, group in sorted(bins[suite].items()):
                rows.append(dict(suite=suite, cell=cell, stratum=name, families=len(group), family_members=group,
                    status="descriptive" if group else "empty_bin",
                    prediction_semantics="resolved_native_pairs" if cell.endswith("r1") else "group_clique",
                    **aggregate([values[cell, f] for f in group])))
    index = {(r["suite"], r["cell"], r["stratum"]): r for r in rows}
    for suite in SUITES:
        for contrast, candidate, reference in CONTRASTS:
            for name, group in sorted(bins[suite].items()):
                differences.append(dict(suite=suite, contrast=contrast, candidate=candidate, reference=reference,
                    stratum=name, families=len(group), family_members=group,
                    status="descriptive" if group else "empty_bin", **{m:
                        index[suite, candidate, name][m] - index[suite, reference, name][m] if group else None for m in METRICS}))
    diff_index = {(r["suite"], r["contrast"], r["stratum"]): r for r in differences}
    for key, suite in (("fixed", None), ("distance", "model_distance")):
        old = docs[key]
        require(len(old["rows"]) == (60 if suite is None else 9)
                and len(old["differences"]) == (40 if suite is None else 6), "Changed prior row inventory")
        old_scores, old_diffs = set(), set()
        for row in old["rows"]:
            ident = (suite or row["suite"], row["cell"], row["stratum"])
            require(ident not in old_scores, "Duplicated prior score")
            old_scores.add(ident)
            new = index[ident]
            require(all(new[k] == row[k] for k in ("families", "family_members", "status", "prediction_semantics")),
                    "Changed prior score identity")
            same({m: row[m] for m in METRICS}, {m: new[m] for m in METRICS}, "Prior bin score does not reproduce")
        for row in old["differences"]:
            ident = (suite or row["suite"], row["contrast"], row["stratum"])
            require(ident not in old_diffs, "Duplicated prior difference")
            old_diffs.add(ident)
            new = diff_index[ident]
            require(all(new[k] == row[k] for k in ("candidate", "reference", "families", "family_members", "status")),
                    "Changed prior contrast identity")
            prior = {m: row[m] if suite is None else None if row[m+"_pp"] is None else row[m+"_pp"]/100 for m in METRICS}
            same(prior, {m: new[m] for m in METRICS}, "Prior contrast/unit conversion does not reproduce")
    return dict(memberships=members, cells=sources, bins=bins, family_rows=family_rows, rows=rows,
        differences=differences, prior_score_rows_reproduced=69, prior_difference_rows_reproduced=46,
        inherited_distance_difference_conversion="original percentage-point fields divided by 100 for raw-unit checks")


def table(result):
    lookup = {(r["suite"], r["cell"], r["stratum"]): r for r in result["rows"]}
    changes = {(r["suite"], r["contrast"], r["stratum"]): r for r in result["differences"]}
    lines = ["# Four-Cell Native SwissTrees Fixed Strata", "",
        "Scores in percent; differences in percentage points. NA is an empty bin, not zero.",
        "R at P0/C0, C at P0/R0, P at C0/R1. Initial HMM search stays on; no subgroup intervals."]
    for suite in SUITES:
        lines.extend(["", "## " + suite.replace("_", " ").title(), "",
            "| Bin | Families | P0/C0/R0 F1 | P0/C0/R1 F1 | P0/C1/R0 F1 | P1/C0/R1 F1 | R F1 | R PPV | R TPR | C F1 | C PPV | C TPR | P F1 | P PPV | P TPR |",
            "|" + "---|" * 15])
        for name, group in sorted(result["bins"][suite].items()):
            values = [lookup[suite, c, name]["F1"] for c in CELLS]
            values.extend(changes[suite, contrast, name][m] for contrast, _, _ in CONTRASTS for m in METRICS)
            rendered = ["NA" if value is None else format(100*value, ".3f" if i < 4 else "+.3f") for i, value in enumerate(values)]
            lines.append("| " + " | ".join([name, str(len(group)), *rendered]) + " |")
    return "\n".join([*lines, "", *["- " + note for note in LIMITATIONS], ""])


def export(repo, output):
    output = Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    docs, refs, checked, protocol = prepare(repo)
    result = dict(schema="native_qfo_four_cell_strata_v1", source=record(__file__), inputs=refs,
        protocol=protocol, checked_inputs=checked, limitations=LIMITATIONS, **SCOPE, **project(docs))
    for ref in [*checked, result["source"]]:
        require(record(ref["path"]) == ref, "Direct evidence changed before writing")
    output.mkdir(parents=True)
    write_tsv(output / "scores.tsv", result["rows"], SCORE_FIELDS)
    write_tsv(output / "differences.tsv", result["differences"], DIFF_FIELDS)
    flat = [dict(cell=r["cell"], family=r["family"], **r["counts_without_prior"], **{m:r[m] for m in METRICS}) for r in result["family_rows"]]
    write_tsv(output / "family_counts.tsv", flat, FAMILY_FIELDS)
    (output / "TABLE.md").write_text(table(result))
    result["outputs"] = [record(output / name) for name in ("scores.tsv", "differences.tsv", "family_counts.tsv", "TABLE.md")]
    for ref in [*checked, result["source"]]:
        require(record(ref["path"]) == ref, "Direct evidence changed while writing")
    with (output / "report.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = export(args.repo, args.output)
    print(json.dumps({key:len(result[key]) for key in ("family_rows", "rows", "differences")}))
