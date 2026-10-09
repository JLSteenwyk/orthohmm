"""Independently check the four-cell projection with exact rational arithmetic."""

import argparse
import csv
from fractions import Fraction
import hashlib
import io
import json
import math
from pathlib import Path


CELLS = ("p0_c0_r0", "p0_c0_r1", "p0_c1_r0", "p1_c0_r1")
SUITES = ("sequence", "domain", "duplication", "model_distance")
METRICS = ("F1", "PPV", "TPR")
COUNTS = ("TP", "FP", "FN", "TN")
CONTRASTS = (("R_at_P0_C0", CELLS[1], CELLS[0]),
             ("C_at_P0_R0", CELLS[2], CELLS[0]), ("P_at_C0_R1", CELLS[3], CELLS[1]))
PINS = {
    "fixed": "55ed72ce771c756424877cfa412da329ff104efb7dca4b35113b0d80703a59e5",
    "fixed_reader": "6f3eef63c0e96e1d85805fcfb70ed9b66fbb1db76a08ca2a105483422a0f6969",
    "distance": "ddfa27757c4c256df9c881079ac352c37d7e17646f75fbd417b624e7980be78f",
    "distance_reader": "186c0151137fd825513ef36c27df595b1a0f1cd4b2b4bf68de4b9665ad328d0e",
    "profile": "d5bafed830bf93e46b810728df14ffa34185d58ffe3eac78a693b0159d3d0a40",
    "profile_reader": "3f62b487014cb9ba99ea36b5adc15535f62c7e39c9a99f3cb4eb8c537523005d",
    "snapshot": "863ebc334dbdb5a41d1a9e9bba56be82896f0ff26c98366047928c7ef802966f",
}
PROTOCOL_SHA = "fa512f814404bca3d88436b2484b32bda03a64d4efebe039a7e551f2791de58c"
SCOPE = ("new_uncertainty", "independent_confirmation", "new_accuracy_or_resource_admission",
         "scientific_timings_admitted", "raw_evidence_reparsed", "original_tree_traversal_repeated", "publication_ready")


def need(condition, message):
    if not condition:
        raise ValueError(message)


def record(path):
    path = Path(path).resolve(strict=True)
    data = path.read_bytes()
    return dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())


def checked(ref):
    need(record(ref["path"]) == ref, "Changed projection source/input/output")


def close(value, exact):
    if exact is None:
        need(value is None, "Empty statistic was imputed")
    else:
        need(type(value) in (int, float) and math.isfinite(value)
             and abs(value-float(exact)) <= 1e-12, "Rational statistic does not match")


def family_statistics(raw):
    need(set(raw) == set(COUNTS) and all(type(v) is int and v >= 0 for v in raw.values())
         and sum(raw.values()) > 0, "Invalid integer family counts")
    p = Fraction(raw["TP"]+2, raw["TP"]+raw["FP"]+4)
    r = Fraction(raw["TP"]+2, raw["TP"]+raw["FN"]+4)
    return dict(F1=2*p*r/(p+r), PPV=p, TPR=r)


def macro(values):
    if not values:
        return dict.fromkeys(METRICS)
    p = sum((v["PPV"] for v in values), Fraction())/len(values)
    r = sum((v["TPR"] for v in values), Fraction())/len(values)
    return dict(F1=2*p*r/(p+r), PPV=p, TPR=r)


def verify(result, docs):
    need(result["schema"] == "native_qfo_four_cell_strata_v1"
         and type(result["new_bootstrap_draws"]) is int and result["new_bootstrap_draws"] == 0
         and all(result[k] is False for k in SCOPE)
         and result["prior_score_rows_reproduced"] == 69 and result["prior_difference_rows_reproduced"] == 46,
         "Inflated or changed projection scope")
    fixed, distance, audit = (docs[k] for k in ("fixed", "distance", "profile"))
    members = fixed["memberships"]
    families = sorted(members)
    genes = [g for f in families for g in members[f]]
    need(len(families) == 18 and len(genes) == len(set(genes)) == 563
         and all(members[f] == sorted(set(members[f])) for f in families)
         and members == result["memberships"] == distance["memberships"], "Canonical members differ")
    bins = dict(fixed["bins"], model_distance=distance["bins"])
    need(result["bins"] == bins and set(bins) == set(SUITES)
         and [len(bins[s]) for s in SUITES] == [11,5,4,3]
         and all(bins[s]["all"] == families and all(v == sorted(set(v)) and set(v) <= set(families)
                  for v in bins[s].values()) for s in SUITES), "Frozen bin membership differs")
    need(fixed["family_rows"] == distance["family_rows"] and fixed["cells"] == distance["cells"],
         "Inherited counts or identities disagree")
    need(audit["selected_index"] == 10 and audit["families"] == families and audit["reference_relation_count"] == 10765
         and [r["cell"] for r in audit["cells"]] == [CELLS[3]], "Changed profile audit")
    profile = audit["cells"][0]
    sources = [*fixed["cells"], {k:profile[k] for k in ("cell","index","native_job_id","admission","raw_file")}]
    sources[-1]["scientific_timings_admitted"] = False
    need(result["cells"] == sources and sources[1]["timing_admitted"] is sources[1]["timing_eligible"] is False,
         "Changed native identity or failed timing")
    snapshot = {r["cell"]:r for r in docs["snapshot"]["rows"]}
    need({c for c,r in snapshot.items() if r["accuracy_admitted"]} == set(CELLS), "Missing cell was admitted")
    for source, job in zip(sources, (22435,22437,22444,23902)):
        target = snapshot[source["cell"]]
        need(source["native_job_id"] == target["native_job_id"] == job
             and source["admission"] == target["admission"], "Native admission does not match")
    need(snapshot[CELLS[1]]["resources"] is None
         and snapshot[CELLS[1]]["timing_admitted"] is snapshot[CELLS[1]]["timing_eligible"] is False
         and snapshot[CELLS[3]]["scientific_timings_admitted"] is False, "Scientific timing was repaired")
    originals = {(r["cell"],r["family"]):r["counts_without_prior"] for r in fixed["family_rows"]}
    for row in profile["families"]:
        need(row["represented_genes"] == members[row["family"]], "Profile members differ")
        originals[CELLS[3],row["family"]] = row["counts_without_prior"]
        for m,v in family_statistics(row["counts_without_prior"]).items():
            close(row["statistics_with_prior"][m],v)
    need([(r["cell"],r["family"]) for r in result["family_rows"]] ==
         [(c,f) for c in CELLS for f in families] and len(originals) == 72, "Family rows incomplete/duplicated")
    values = {}
    for row in result["family_rows"]:
        key = row["cell"],row["family"]
        need(row["counts_without_prior"] == originals[key], "Projected raw count changed")
        values[key] = family_statistics(originals[key])
        for m in METRICS:
            close(row[m],values[key][m])
    for f in families:
        truths = {(originals[c,f]["TP"]+originals[c,f]["FN"], originals[c,f]["FP"]+originals[c,f]["TN"]) for c in CELLS}
        need(len(truths) == 1, "Reference truth universe differs")
    need(sum(sum(originals[CELLS[0],f].values()) for f in families) == 10765, "Relation count differs")
    for c in CELLS:
        point = macro([values[c,f] for f in families])
        need(abs(float(point["F1"])-snapshot[c]["scores"]["SwissTrees"]) <= 5e-8, "Native endpoint differs")
        if c in (CELLS[1],CELLS[3]):
            for m in METRICS:
                close(docs["profile_reader"]["rational_macro_points"][c][m],point[m])
    for m in METRICS:
        close(profile["aggregate"][m],macro([values[CELLS[3],f] for f in families])[m])
    expected_scores = [(s,c,n) for s in SUITES for c in CELLS for n in sorted(bins[s])]
    need([(r["suite"],r["cell"],r["stratum"]) for r in result["rows"]] == expected_scores,
         "Score row inventory/order differs")
    scores = {}
    for row in result["rows"]:
        suite, cell, name = row["suite"],row["cell"],row["stratum"]
        group = bins[suite][name]
        need(row["family_members"] == group and row["families"] == len(group)
             and row["status"] == ("descriptive" if group else "empty_bin")
             and row["prediction_semantics"] == ("resolved_native_pairs" if cell.endswith("r1") else "group_clique"),
             "Score identity/semantics differs")
        scores[suite,cell,name] = macro([values[cell,f] for f in group])
        for m in METRICS:
            close(row[m],scores[suite,cell,name][m])
    expected_diffs = [(s,contrast,n) for s in SUITES for contrast,_,_ in CONTRASTS for n in sorted(bins[s])]
    need([(r["suite"],r["contrast"],r["stratum"]) for r in result["differences"]] == expected_diffs,
         "Difference inventory/order differs")
    differences = {}
    for row in result["differences"]:
        suite,name = row["suite"],row["stratum"]
        candidate,reference = next((c,r) for label,c,r in CONTRASTS if label == row["contrast"])
        group = bins[suite][name]
        need(row["candidate"] == candidate and row["reference"] == reference
             and row["family_members"] == group and row["families"] == len(group)
             and row["status"] == ("descriptive" if group else "empty_bin"), "Contrast identity differs")
        values_exact = {m:scores[suite,candidate,name][m]-scores[suite,reference,name][m] if group else None for m in METRICS}
        differences[suite,row["contrast"],name] = values_exact
        for m in METRICS:
            close(row[m],values_exact[m])
    for key,suite in (("fixed",None),("distance","model_distance")):
        for row in docs[key]["rows"]:
            actual_suite = suite or row["suite"]
            for m in METRICS:
                close(row[m],scores[actual_suite,row["cell"],row["stratum"]][m])
        for row in docs[key]["differences"]:
            actual_suite = suite or row["suite"]
            for m in METRICS:
                value = row[m] if suite is None else None if row[m+"_pp"] is None else row[m+"_pp"]/100
                close(value,differences[actual_suite,row["contrast"],row["stratum"]][m])
    return values,scores,differences


def check_tables(directory, result, values, scores, differences):
    layouts = (("scores.tsv",result["rows"],("suite","cell","stratum","families","status",*METRICS,"prediction_semantics")),
               ("differences.tsv",result["differences"],("suite","contrast","candidate","reference","stratum","families","status",*METRICS)))
    for name,records,fields in layouts:
        parser = csv.DictReader(io.StringIO((directory/name).read_text()),delimiter="\t")
        need(parser.fieldnames == list(fields), "TSV columns differ")
        rows = list(parser)
        need(len(rows) == len(records), "TSV row count differs")
        for text,row in zip(rows,records):
            for field in fields:
                if field in METRICS:
                    if row[field] is None:
                        need(text[field] == "NA", "Empty TSV value was imputed")
                    else:
                        close(float(text[field]),row[field])
                else:
                    need(text[field] == str(row[field]), "TSV identity differs")
    parser = csv.DictReader(io.StringIO((directory/"family_counts.tsv").read_text()),delimiter="\t")
    need(parser.fieldnames == ["cell","family",*COUNTS,*METRICS], "Family TSV columns differ")
    rows = list(parser)
    need(len(rows) == 72, "Family TSV inventory differs")
    for text,row in zip(rows,result["family_rows"]):
        need(text["cell"] == row["cell"] and text["family"] == row["family"]
             and all(text[k] == str(row["counts_without_prior"][k]) for k in COUNTS), "Family TSV counts differ")
        for m in METRICS:
            close(float(text[m]),values[row["cell"],row["family"]][m])
    markdown = (directory/"TABLE.md").read_text()
    expected_lines = []
    for suite in SUITES:
        for name,group in sorted(result["bins"][suite].items()):
            numbers = [scores[suite,c,name]["F1"] for c in CELLS]
            numbers.extend(differences[suite,label,name][m] for label,_,_ in CONTRASTS for m in METRICS)
            text = ["NA" if value is None else format(100*float(value),".3f" if i<4 else "+.3f") for i,value in enumerate(numbers)]
            expected = "| "+" | ".join([name,str(len(group)),*text])+" |"
            expected_lines.append(expected)
    actual_lines = [line for line in markdown.splitlines() if line.startswith("| ") and not line.startswith("| Bin ")]
    need(actual_lines == expected_lines and len(actual_lines) == 23,
         "Human table inventory/ordered rational values differ")
    return dict(score_tsv_rows=92,difference_tsv_rows=69,family_tsv_rows=72,human_rows_checked=len(expected_lines))


def readback(path, sha, output):
    path, output = Path(path).resolve(), Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    report = record(path)
    need(report["sha256"] == sha, "Changed selected projection report")
    result = json.loads(path.read_text())
    need(set(result["inputs"]) == set(PINS) and result["protocol"]["sha256"] == PROTOCOL_SHA, "Changed input/protocol pins")
    refs = [report,result["source"],result["protocol"],*result["checked_inputs"],*result["outputs"]]
    docs = {}
    for key,digest in PINS.items():
        ref = result["inputs"][key]
        need(ref["sha256"] == digest and ref in result["checked_inputs"], "Changed scientific input pin")
        checked(ref)
        docs[key] = json.loads(Path(ref["path"]).read_text())
        need(docs[key]["source"] in result["checked_inputs"], "Inherited source not checked")
    need(docs["fixed_reader"]["report"] == result["inputs"]["fixed"]
         and docs["distance_reader"]["report"] == result["inputs"]["distance"]
         and docs["profile_reader"]["audit"] == result["inputs"]["profile"], "Wrong original reader binding")
    need([Path(r["path"]).name for r in result["outputs"]] == ["scores.tsv","differences.tsv","family_counts.tsv","TABLE.md"]
         and all(Path(r["path"]).parent == path.parent for r in result["outputs"]), "Unbound projection output")
    for ref in refs:
        checked(ref)
    values,scores,differences = verify(result,docs)
    tables = check_tables(path.parent,result,values,scores,differences)
    for ref in refs:
        checked(ref)
    receipt = dict(schema="native_qfo_four_cell_strata_rational_readback_v1", report=report,
        source=record(__file__), checked_inputs=refs, families_checked=18, proteins_checked=563,
        family_rows_checked=72, score_rows_checked=92, differences_checked=69,
        inherited_score_rows_reproduced=69, inherited_difference_rows_reproduced=46,
        **tables, new_bootstrap_draws=0, independent_confirmation=False, publication_ready=False,
        new_uncertainty=False,new_accuracy_or_resource_admission=False,scientific_timings_admitted=False,
        raw_evidence_reparsed=False,original_tree_traversal_repeated=False,
        limitations=["Exact retained-record arithmetic, not repeated raw scoring or transitive admission.",
                     "Fixed overlapping bins are descriptive; no subgroup interval, causal effect or new default."])
    with output.open("x") as stream:
        json.dump(receipt,stream,indent=2,sort_keys=True,allow_nan=False)
        stream.write("\n")
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report",required=True,type=Path)
    parser.add_argument("--report-sha256",required=True)
    parser.add_argument("--output",required=True,type=Path)
    args = parser.parse_args()
    result = readback(args.report,args.report_sha256,args.output)
    print(json.dumps({k:result[k] for k in ("family_rows_checked","score_rows_checked","differences_checked")}))
