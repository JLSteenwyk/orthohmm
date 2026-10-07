"""Versioned path-normalized independent readback; preserves the failed v1 source."""

import argparse
import csv
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

CELLS = ("p0_c0_r0", "p0_c0_r1", "p0_c1_r0")
METRICS = ("F1", "PPV", "TPR")
SUITES = ("sequence", "domain", "duplication")
CONTRASTS = (("R_at_P0_C0", CELLS[1]), ("C_at_P0_R0", CELLS[2]))
PINS = {
    "sequence": ("native_qfo_swiss_sequence_strata_20261006_v1/report.json", "4650f614d2fdeebd65cbcbf4999c61c83cc549560d2738bfbb5502f686f660e3"),
    "sequence_reader": ("native_qfo_swiss_sequence_strata_readback_20261006.json", "8d0356bc0e159d204939271578e23a6a5b6efcf6198b75a18c0f3d38e0e67c0a"),
    "domain": ("native_qfo_swiss_domain_strata_20261006_v1/report.json", "d3235c712e711f005b3182a3621c0079f8a85c3c7abd5db016d5e8ea95884cd6"),
    "domain_reader": ("native_qfo_swiss_domain_strata_readback_20261006_v1.json", "9e77b67ad97e6824e9442691be7a05ba1552d62ffe696ced6034a996ede1945d"),
    "duplication": ("native_qfo_swiss_duplication_strata_20261006_v1/report.json", "e0ecbc9d26738d0d3d573ca31772085b28d145ba2deaa9383c8cff75bd0d67b5"),
    "duplication_reader": ("native_qfo_swiss_duplication_strata_readback_20261006_v1.json", "84eb3d9fe02b5984cecaa06355d1cb1d68c54aec09ab588eef8b3ada42f1957b"),
    "features": ("corrected_swiss_sequence_strata_20260918.json", "912363276fa5456f11aeff9f398fce71aec8d377c8c8eda7b144105945e554f1"),
    "candidate": ("native_qfo_candidate_swiss_counts_20261006_v1.json", "7ed2ea1a7a94d55b28ecbadd7bd2be8e99df64b9fab3c266ebd987e8dbcfc479"),
    "candidate_reader": ("native_qfo_candidate_swiss_readback_20261006_v1.json", "07941ba9c55ba7dc14c5ba53df37d28cb3b2f60f599b9d177b10c9afa99c1329"),
    "protocol": ("NATIVE_QFO_THREE_CELL_STRATA_PROTOCOL_20261007.md", "9e14a818df91a4ab12d6db5c060e556cf872c4eb7d45dfc43f8331fcbbc685ca"),
}
REPAIR_PINS = {
    "failure": ("native_qfo_three_cell_strata_readback_failure_20261007_v1.json",
                "3bcf3f6eafaf68a1e98a41d00ddf3579497ccaff4c8f553012f8a220d8eaef66"),
    "amendment": ("NATIVE_QFO_THREE_CELL_STRATA_READER_AMENDMENT_20261007.md",
                  "768929b95f482c02e6e38a2d8feaf5b1f872c3376215cd0686f814f47cd2970f"),
}
FLAGS = ("new_uncertainty", "new_accuracy_or_resource_admission", "independent_confirmation",
         "publication_ready", "scientific_timings_admitted", "raw_evidence_reparsed", "original_tree_traversal_repeated")
SCORE_FIELDS = ("suite", "cell", "stratum", "families", "status", *METRICS, "prediction_semantics")
DIFF_FIELDS = ("suite", "contrast", "candidate", "reference", "stratum", "families", "status", *METRICS)
FAMILY_FIELDS = ("cell", "family", "TP", "FP", "FN", "TN", *METRICS)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def record(path):
    path = Path(path).resolve(strict=True)
    sha = hashlib.sha256()
    size = 0
    with path.open("rb") as stream:
        while True:
            data = stream.read(65536)
            if not data:
                break
            size += len(data)
            sha.update(data)
    return {"path": str(path), "sha256": sha.hexdigest(), "bytes": size}


def near(value, expected):
    if expected is None:
        require(value is None, "Empty metric relabeled")
    else:
        require(type(value) in (int, float) and math.isfinite(value)
                and abs(value - float(expected)) <= 1e-12, "Metric does not reproduce")


def points(counts):
    require(set(counts) == {"TP", "FP", "FN", "TN"}
            and all(type(n) is int and n >= 0 for n in counts.values())
            and sum(counts.values()) > 0, "Invalid original family counts")
    p = Fraction(counts["TP"] + 2, counts["TP"] + counts["FP"] + 4)
    r = Fraction(counts["TP"] + 2, counts["TP"] + counts["FN"] + 4)
    return {"PPV": p, "TPR": r, "F1": 2 / (1 / p + 1 / r)}


def macro(values):
    if not values:
        return dict.fromkeys(METRICS)
    p = sum((v["PPV"] for v in values), Fraction()) / len(values)
    r = sum((v["TPR"] for v in values), Fraction()) / len(values)
    return {"PPV": p, "TPR": r, "F1": 2 / (1 / p + 1 / r)}


def validate(docs, report):
    require(report["schema"] == "native_qfo_three_cell_strata_v1"
            and all(report[k] is False for k in FLAGS)
            and type(report["new_bootstrap_draws"]) is int and report["new_bootstrap_draws"] == 0,
            "Inflated projection scope")
    members = docs["features"]["family_memberships"]
    families = sorted(members)
    genes = [g for f in families for g in members[f]]
    require(len(families) == 18 and len(genes) == len(set(genes)) == 563
            and all(members[f] == sorted(set(members[f])) for f in families)
            and members == docs["domain"]["memberships"] == docs["duplication"]["memberships"] == report["memberships"],
            "Canonical member inventory changed")
    require(docs["features"]["prediction_statistics_evaluated"] is False, "Prediction-dependent features")
    old = docs["sequence"]["family_rows"]
    require([(r["cell"], r["family"]) for r in old] == [(c, f) for c in CELLS[:2] for f in families],
            "Prior family inventory changed")
    for suite in SUITES:
        require(docs[suite]["family_rows"] == old and docs[suite]["schema"] ==
                "native_qfo_swiss_" + suite + "_strata_v1", "Prior reports disagree")
        require(all(docs[suite][k] is False for k in FLAGS[:4]), "Inflated inherited scope")
        require(docs[suite + "_reader"]["family_rows_checked"] == 36
                and docs[suite + "_reader"]["raw_rows_checked"] == 21530, "Incomplete prior readback")
    audit = docs["candidate"]
    require(audit["families"] == families and audit["reference_relation_count"] == 10765
            and audit["schema"] == "native_qfo_swiss_family_count_audit_v1"
            and audit["status"] == "supplied_native_swiss_family_counts_verified"
            and len(audit["cells"]) == 1 and audit["cells"][0]["cell"] == CELLS[2]
            and audit["new_bootstrap_draws"] == 0 and type(audit["new_bootstrap_draws"]) is int
            and all(audit[k] is False for k in ("historical_intervals_attached", *FLAGS[1:4])),
            "Candidate audit identity/scope changed")
    candidate = audit["cells"][0]
    inherited = docs["candidate_reader"]
    require(inherited["cells"] == [CELLS[0], CELLS[2]] and inherited["families_checked"] == 18
            and inherited["native_family_records_checked"] == 36
            and inherited["candidate_pair_labels_matched"] == 10765, "Candidate readback incomplete")
    require(candidate["index"] == 8 and candidate["native_job_id"] == 22444
            and candidate["retained_family_records_identical"] is True
            and candidate["retained_aggregate_identical"] is True
            and [r["family"] for r in candidate["families"]] == families,
            "Candidate identity changed")
    counts = {(r["cell"], r["family"]): r["counts_without_prior"] for r in old}
    for row in candidate["families"]:
        require(row["represented_genes"] == members[row["family"]], "Candidate protein inventory changed")
        counts[CELLS[2], row["family"]] = row["counts_without_prior"]
    family_values = {k: points(v) for k, v in counts.items()}
    for row in old:
        for m in METRICS:
            near(row[m], family_values[row["cell"], row["family"]][m])
    for row in candidate["families"]:
        for m in METRICS:
            near(row["statistics_with_prior"][m], family_values[CELLS[2], row["family"]][m])
    for f in families:
        require(len({(sum(counts[c, f].values()), counts[c, f]["TP"] + counts[c, f]["FN"])
                     for c in CELLS}) == 1, "Classified/reference-positive universe changed")
    require(sum(sum(counts[CELLS[0], f].values()) for f in families) == 10765, "Wrong reference total")
    for c in (CELLS[0], CELLS[2]):
        value = macro([family_values[c, f] for f in families])
        for m in METRICS:
            near(inherited["rational_macro_points"][c][m], value[m])
            if c == CELLS[2]:
                near(candidate["aggregate"][m], value[m])
    require([(r["cell"], r["family"]) for r in report["family_rows"]] ==
            [(c, f) for c in CELLS for f in families], "Wrong exported family inventory")
    for row in report["family_rows"]:
        key = (row["cell"], row["family"])
        require(set(row) == {"cell", "family", "counts_without_prior", *METRICS}
                and row["counts_without_prior"] == counts[key], "Changed exported counts")
        for m in METRICS:
            near(row[m], family_values[key][m])
    expected_cells = [dict(c) for c in docs["domain"]["cells"]]
    require([c["cell"] for c in expected_cells] == list(CELLS[:2])
            and expected_cells[1]["timing_admitted"] is False
            and expected_cells[1]["timing_eligible"] is False, "Failed timing repaired")
    expected_cells.append({k: candidate[k] for k in ("cell", "index", "native_job_id", "raw_file", "admission")})
    require(report["cells"] == expected_cells, "Cell provenance changed")
    expected_rows, expected_diffs, expected_bins = [], [], {}
    for suite, size in zip(SUITES, (11, 5, 4)):
        rows = docs[suite]["rows"]
        baseline = [r for r in rows if r["cell"] == CELLS[0]]
        names = [r["stratum"] for r in baseline]
        bins = {r["stratum"]: r["family_members"] for r in baseline}
        require(len(names) == len(set(names)) == size and bins["all"] == families
                and all(v == sorted(set(v)) and set(v) <= set(families) for v in bins.values()),
                "Frozen bins incomplete/invalid")
        frozen = (dict(all=families, **docs["features"]["primary_strata"], **docs["features"]["secondary_strata"])
                  if suite == "sequence" else docs[suite]["bins"])
        require(bins == frozen and [(r["cell"], r["stratum"]) for r in rows] ==
                [(c, n) for c in CELLS[:2] for n in names], "Frozen bins changed")
        expected_bins[suite] = bins
        old_lookup = {(r["cell"], r["stratum"]): r for r in rows}
        values = {}
        for c in CELLS:
            for n, fs in bins.items():
                point = macro([family_values[c, f] for f in fs])
                values[c, n] = point
                semantics = "resolved_native_pairs" if c == CELLS[1] else "group_clique"
                if c in CELLS[:2]:
                    prior = old_lookup[c, n]
                    require(prior["family_members"] == fs and prior["families"] == len(fs)
                            and prior["prediction_semantics"] == semantics, "Prior bin unit/semantics changed")
                    for m in METRICS:
                        near(prior[m], point[m])
                expected_rows.append(dict(suite=suite, cell=c, stratum=n, families=len(fs), family_members=fs,
                    status="descriptive" if fs else "empty_bin", prediction_semantics=semantics, **point))
        for contrast, c in CONTRASTS:
            for n, fs in bins.items():
                expected_diffs.append(dict(suite=suite, contrast=contrast, candidate=c, reference=CELLS[0], stratum=n,
                    families=len(fs), family_members=fs, status="descriptive" if fs else "empty_bin",
                    **{m: values[c, n][m] - values[CELLS[0], n][m] if fs else None for m in METRICS}))
    require(report["bins"] == expected_bins, "Output bins altered")
    for key, expected in (("rows", expected_rows), ("differences", expected_diffs)):
        require(len(report[key]) == len(expected), "Missing/duplicate projected row")
        for row, truth in zip(report[key], expected):
            require(set(row) == set(truth) and all(row[k] == v for k, v in truth.items() if k not in METRICS),
                    "Projected row identity changed")
            for m in METRICS:
                near(row[m], truth[m])
    return dict(families_checked=18, proteins_checked=563, family_rows_checked=54,
                score_rows_checked=len(expected_rows), differences_checked=len(expected_diffs),
                inherited_score_rows_reproduced=40)


def tsv(path, expected, fields):
    with Path(path).open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        require(reader.fieldnames == list(fields), "TSV header altered")
        rows = list(reader)
    require(len(rows) == len(expected), "TSV row inventory altered")
    for row, truth in zip(rows, expected):
        require(set(row) == set(fields), "TSV row width altered")
        for k in fields:
            v = truth[k]
            if k in METRICS and v is not None:
                try:
                    near(float(row[k]), v)
                except (TypeError, ValueError):
                    raise ValueError("TSV metric altered") from None
            else:
                require(row[k] == ("NA" if v is None else str(v)), "TSV identity/count altered")


def human_table(path, report):
    lines = Path(path).read_text().splitlines()
    scores = {(r["suite"], r["cell"], r["stratum"]): r for r in report["rows"]}
    differences = {(r["suite"], r["contrast"], r["stratum"]): r for r in report["differences"]}
    suite, observed = None, []
    for line in lines:
        if line.startswith("## "):
            suite = line[3:].lower()
            require(suite in SUITES, "Unknown human-table suite")
        if line.startswith("| ") and not line.startswith("| Bin |"):
            fields = [s.strip() for s in line.strip("|").split("|")]
            require(suite is not None and len(fields) == 11, "Human-table row altered")
            n = fields[0]
            require(n in report["bins"][suite], "Unknown human-table bin")
            require(fields[1] == str(len(report["bins"][suite][n])), "Human-table family count altered")
            numbers = [scores[suite, c, n]["F1"] for c in CELLS]
            numbers.extend(differences[suite, contrast, n][m] for contrast, _ in CONTRASTS for m in METRICS)
            expected = ["NA" if v is None else format(float(v) * 100, ".3f" if i < 3 else "+.3f")
                        for i, v in enumerate(numbers)]
            require(fields[2:] == expected, "Human-table scores altered")
            observed.append((suite, n))
    require(observed == [(s, n) for s in SUITES for n in report["bins"][s]], "Human-table inventory altered")
    require([line[2:] for line in lines if line.startswith("- ")] == report["limitations"],
            "Human-table limitations changed")


def verify(report_path, expected_sha):
    report_ref = record(report_path)
    require(report_ref["sha256"] == expected_sha, "Changed supplied report")
    report = json.loads(Path(report_path).read_text())
    require(set(report["inputs"]) == set(PINS), "Input identity inventory changed")
    docs, checked = {}, [report_ref]
    for key, (name, sha) in PINS.items():
        ref = report["inputs"][key]
        require(ref["sha256"] == sha and Path(ref["path"]).as_posix().endswith("/benchmark_tools/results/" + name)
                and record(ref["path"]) == ref, "Changed direct pinned evidence")
        checked.append(ref)
        if name.endswith(".json"):
            docs[key] = json.loads(Path(ref["path"]).read_text())
    sources = []
    for key in (*SUITES, *(s + "_reader" for s in SUITES), "features", "candidate", "candidate_reader"):
        ref = docs[key]["source"]
        require(record(ref["path"]) == ref, "Changed inherited source")
        if ref not in sources:
            sources.append(ref)
    require(report["checked_inputs"] == [report["inputs"][k] for k in PINS] + sources,
            "Checked inputs relabeled or omitted")
    checked.extend(sources)
    exporter = Path(__file__).with_name("export_native_qfo_three_cell_strata.py")
    require(report["source"] == record(exporter), "Changed exporter source")
    checked.append(report["source"])
    for suite in SUITES:
        require(docs[suite + "_reader"]["report"] == report["inputs"][suite], "Wrong prior readback binding")
    require(docs["candidate_reader"]["audit"] == report["inputs"]["candidate"]
            and docs["sequence"]["strata"] == report["inputs"]["features"], "Wrong candidate/features binding")
    result = validate(docs, report)
    require(isinstance(report["limitations"], list) and len(report["limitations"]) == 12
            and all(type(s) is str and s for s in report["limitations"]), "Limitations omitted")
    require(len(report["outputs"]) == 4 and [Path(r["path"]).name for r in report["outputs"]] ==
            ["scores.tsv", "differences.tsv", "family_counts.tsv", "TABLE.md"], "Output inventory changed")
    for ref in report["outputs"]:
        require(Path(ref["path"]).parent.resolve() == Path(report_path).parent.resolve()
                and record(ref["path"]) == ref, "Changed output bytes/location")
    checked.extend(report["outputs"])
    tsv(report["outputs"][0]["path"], report["rows"], SCORE_FIELDS)
    tsv(report["outputs"][1]["path"], report["differences"], DIFF_FIELDS)
    flattened = [{"cell": r["cell"], "family": r["family"], **r["counts_without_prior"],
                  **{m: r[m] for m in METRICS}} for r in report["family_rows"]]
    tsv(report["outputs"][2]["path"], flattened, FAMILY_FIELDS)
    human_table(report["outputs"][3]["path"], report)
    for key, (name, sha) in REPAIR_PINS.items():
        ref = record(Path(__file__).parent / "results" / name)
        require(ref["sha256"] == sha, "Changed readback repair evidence")
        checked.append(ref)
        if key == "failure":
            failure = json.loads(Path(ref["path"]).read_text())
            require(failure["report"] == report_ref and failure["execution_exit_code"] == 1
                    and failure["output_written"] is False
                    and failure["report_or_export_overwritten"] is False
                    and failure["source"]["sha256"] == "f59bd8cc924e5cbd0a6ef253a2ed47c243e92d1aaf2155dd3bc1ac286390583f"
                    and record(failure["source"]["path"]) == failure["source"],
                    "Original failed readback history changed")
            checked.append(failure["source"])
    for ref in checked:
        require(record(ref["path"]) == ref, "Evidence changed during readback")
    return dict(schema="native_qfo_three_cell_strata_rational_readback_v2", source=record(__file__),
                report=report_ref, checked_inputs=checked, **result, new_bootstrap_draws=0,
                **dict.fromkeys(FLAGS, False), limitations=report["limitations"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), "Output already exists")
    result = verify(args.report, args.report_sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: result[k] for k in ("family_rows_checked", "score_rows_checked", "differences_checked")}))
