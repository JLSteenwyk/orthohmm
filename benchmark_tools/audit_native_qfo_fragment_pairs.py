"""Join fixed historical fragment annotations to all native SwissTrees pairs."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys

import Bio
from Bio import SwissProt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import audit_unisave_fragment_source as annotation_reader
from benchmark_tools import trace_native_qfo_swiss_transitions as native_reader
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

VIEWS = ("historical", "baseline_only")
BINS = ("annotation_positive", "all_matched_unflagged", "missing_without_positive")
LABELS = ("TP", "FP", "FN", "TN")
FIELDS = ("family", "protein_a", "protein_b", "before", "after", "historical", "baseline_only")
PINS = {
    "native_qfo_swiss_pair_transitions_20261006_v2.json":
        "fe003f2cbc4285ea56cd80a92b71703b5c9e135ad1c8f244504dcfeab0911b46",
    "swiss_historical_fragment_admission_22117.json":
        "a480f32666bb96de3395213c7cad024c170318ca110c389f15e4756f8052adf8",
    "corrected_swiss_sequence_strata_20260918.json":
        "912363276fa5456f11aeff9f398fce71aec8d377c8c8eda7b144105945e554f1",
    "native_qfo_swiss_pair_transition_sql_readback_20261006.json":
        "09249d97c6ed3da6c840586552b4f5a72480cc46bbb944bc60124f5382164ffb",
}


def state(annotation, baseline_only=False):
    if annotation is None:
        return None
    native_reader.require(type(annotation["fragment_flag"]) is bool
                          and isinstance(annotation["incomplete_sequence_features"], list)
                          and annotation["selection_class"] in ("baseline_release", "later_sequence_version"),
                          "Invalid admitted annotation state")
    if baseline_only and annotation["selection_class"] != "baseline_release":
        return None
    return annotation["fragment_flag"] or bool(annotation["incomplete_sequence_features"])


def pair_bin(a, b):
    return (BINS[0] if a is True or b is True else BINS[2] if a is None or b is None else BINS[1])


def summarize(view, name, table):
    before = {a: sum(table[a + "->" + b] for b in LABELS) for a in LABELS}
    after = {b: sum(table[a + "->" + b] for a in LABELS) for b in LABELS}
    tp, fp = table["TP->FN"], table["FP->TN"]
    return dict(view=view, bin=name, relations=sum(table.values()), transitions=dict(table),
                before_counts=before, after_counts=after, removed_true_positives=tp, removed_false_positives=fp,
                tp_removal_fraction=tp / before["TP"] if before["TP"] else None,
                fp_removal_fraction=fp / before["FP"] if before["FP"] else None)


def project(before, after, annotations):
    native_reader.require(before.keys() == after.keys(), "Changed native pair universe")
    proteins = {g for (_, a, b) in before for g in (a, b)}
    native_reader.require(proteins == set(annotations), "Changed annotation/native accession universe")
    states = {view: {g: state(a, view == "baseline_only") for g, a in annotations.items()} for view in VIEWS}
    tables = {(v, name): Counter({a + "->" + b: 0 for a in LABELS for b in LABELS}) for v in VIEWS for name in BINS}
    ledger = []
    for key in sorted(before):
        a, b = before[key], after[key]
        native_reader.require(a in LABELS and b in LABELS and (a in ("TP", "FN")) == (b in ("TP", "FN")),
                              "Invalid decisions or changed reference truth")
        names = [pair_bin(states[v][key[1]], states[v][key[2]]) for v in VIEWS]
        for view, name in zip(VIEWS, names):
            tables[view, name][a + "->" + b] += 1
        ledger.append((*key, a, b, *names))
    rows = [summarize(view, name, tables[view, name]) for view in VIEWS for name in BINS]
    return rows, ledger


def audit(repo, output):
    root = Path(repo).resolve() / "benchmark_tools/results"
    output = Path(output)
    native_reader.require(not output.exists() and not output.is_symlink(), "Output already exists")
    checked, refs, documents = [], {}, {}
    for name, sha in PINS.items():
        ref = record(root / name)
        native_reader.require(ref["sha256"] == sha, "Changed frozen input: " + name)
        refs[name], documents[name] = ref, json.loads(Path(ref["path"]).read_text())
        checked.append(ref)
    native, admission, features, readback = [documents[name] for name in PINS]
    native_reader.require(native["schema"] == "native_qfo_swiss_pair_transitions_v1"
                          and native["source"] == record(native_reader.__file__)
                          and readback["report"] == refs[list(PINS)[0]] and readback["raw_rows_checked"] == 21530
                          and all(native[k] is False for k in ("new_scoring_or_admission", "uncertainty_admitted",
                              "scientific_timings_admitted", "independent_confirmation", "publication_ready")),
                          "Wrong native transition source/scope")
    native_reader.require(admission["status"] == "historical_annotation_panel_checked_with_explicit_missingness"
                          and admission["annotation_panel_admitted"] is True and admission["matched"] == 563
                          and admission["missing"] == 0 and admission["prediction_statistics_evaluated"] is False
                          and admission["publication_ready"] is False
                          and admission["families"] == features["family_memberships"], "Wrong annotation admission/membership")
    native_reader.require(all(r in admission["records"] for r in features["fasta_inputs"]),
                          "Historical annotation input identities differ")
    checked.extend(native["checked_inputs"])
    checked.extend([record(annotation_reader.__file__), record(native_reader.__file__), record(SwissProt.__file__),
                    record(root / "NATIVE_QFO_FRAGMENT_PAIR_PROTOCOL_20261006.md")])
    native_reader.require(record(annotation_reader.__file__)["sha256"] ==
                          "0642beeee11199853518b328b15bdb12d297d3a1324bcd3327815ccbc7355680", "Changed annotation reader")
    for ref in checked:
        check(ref)
    annotations = admission["annotations"]
    entries = [r for r in admission["records"] if Path(r["path"]).name == "entry.txt"]
    native_reader.require(len(entries) == len(annotations) == 563 and
                          {Path(r["path"]).parent.name for r in entries} == set(annotations), "Wrong selected entry inventory")
    for ref in entries:
        check(ref)
        gene = Path(ref["path"]).parent.name
        a = annotations[gene]
        observed = annotation_reader.inspect_entry(Path(ref["path"]).read_text(), gene, a["entry_version"],
            a["sequence_version"], a["sequence_sha256"], a["taxid"])
        native_reader.require(all(observed[k] == a[k] for k in observed if k != "limitation"),
                              "Selected raw annotation differs")
        check(ref)
        checked.append(ref)
    families = sorted(admission["families"])
    native_reader.require(len(families) == 18 and [r["cell"] for r in native["cells"]] ==
                          ["p0_c0_r0", "p0_c0_r1"] and native["cells"][1]["resources"] is None
                          and native["cells"][1]["timing_eligible"] is False
                          and native["cells"][1]["timing_admitted"] is False, "Wrong cells or repaired timing")
    labels = []
    for cell in native["cells"]:
        native_reader.require(cell["raw"] in checked and cell["admission"] in checked, "Missing native source binding")
        scientific = json.loads(Path(cell["admission"]["path"]).read_text())
        native_reader.require(scientific["accuracy_admitted"] is True
                              and scientific["conversion"]["input_fastas"] == features["fasta_inputs"],
                              "Native/historical FASTA identities differ")
        label = native_reader.read_labels(Path(cell["raw"]["path"]), families)
        members = {f: set() for f in families}
        for f, a, b in label:
            members[f].update((a, b))
        native_reader.require(all(sorted(members[f]) == admission["families"][f] for f in families),
                              "Native/historical represented genes differ")
        labels.append(label)
    rows, ledger = project(*labels, annotations)
    native_reader.require(len(ledger) == 10765 and all(
        sum(row["transitions"][k] for row in rows if row["view"] == view) == v
        for view in VIEWS for k, v in native["comparison"]["transitions"].items()), "Native transition totals differ")
    for ref in checked:
        check(ref)
    output.mkdir(parents=True)
    with (output / "pairs.tsv").open("x", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(FIELDS)
        writer.writerows(ledger)
    lines = ["# Native Fragment-Annotation Pair Changes", "",
        "Descriptive endpoint-level counts. Unflagged is not proven complete; baseline-only treats later-version annotations as missing.",
        "No pair-IID intervals, significance, native F1 substitution or causal explanation.", "",
        "| View | Endpoint bin | Relations | R0 TP | Removed TP | TP removal % | R0 FP | Removed FP | FP removal % |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for row in rows:
        rates = ["NA" if row[k] is None else f"{100 * row[k]:.3f}" for k in ("tp_removal_fraction", "fp_removal_fraction")]
        lines.append(f"| {row['view']} | {row['bin']} | {row['relations']} | {row['before_counts']['TP']} | "
                     f"{row['removed_true_positives']} | {rates[0]} | {row['before_counts']['FP']} | "
                     f"{row['removed_false_positives']} | {rates[1]} |")
    lines.extend(["", "All 16 transition cells per view/bin, both cell marginals and the full 10,765-pair ledger are retained in report.json/pairs.tsv.",
                  "Annotation selection/history admission is inherited; selected raw entries are newly rechecked. Parsers share Bio.SwissProt.",
                  "Initial HMM search on in both native cells; failed R1 timing remains ineligible. Shared-host postprocessing costs are not inference timing; contention effects are unknown and potentially tool-dependent.", ""])
    (output / "TABLE.md").write_text("\n".join(lines))
    result = dict(schema="native_qfo_fragment_pair_audit_v1", source=record(__file__),
        native=refs[list(PINS)[0]], annotations=refs[list(PINS)[1]], features=refs[list(PINS)[2]],
        original_transition_readback=refs[list(PINS)[3]], checked_inputs=checked, entry_sources=entries,
        source_history_readmission=False, unselected_acquisition_records_rechecked=False,
        source_acquisition_records_inherited=sum(r not in checked for r in admission["records"]),
        annotation_parser="Bio.SwissProt", biopython_version=Bio.__version__, cells=native["cells"],
        rows=rows, paired_relations=len(ledger), raw_rows_checked=2 * len(ledger),
        matched_proteins=563, annotation_positive_proteins=sum(state(a) is True for a in annotations.values()),
        later_version_proteins=sum(a["selection_class"] == "later_sequence_version" for a in annotations.values()),
        new_bootstrap_draws=0, new_scoring_or_admission=False, uncertainty_admitted=False,
        scientific_timings_admitted=False, independent_confirmation=False, publication_ready=False,
        outputs=[record(output / name) for name in ("pairs.tsv", "TABLE.md")],
        limitations=["Retrospective endpoint-level annotations, not experimental truth, causal fragment effects or method tuning.",
                     "Unflagged is not proven complete; a positive family does not make every pair positive.",
                     "Both annotation checks use Bio.SwissProt; not entirely independent parser validation.",
                     "No new sampling law/intervals; dependent reference pairs are not independent binomial trials.",
                     "Original history selection and FASTA consumption admission are reused, not reacquired or wholly re-admitted.",
                     "Shared-host postprocessing only, unknown potentially tool-dependent contention; no isolated speed ranking.",
                     "Failed R1 timing remains ineligible; incomplete native factorial and publication work remain open."])
    for ref in checked:
        check(ref)
    with (output / "report.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.repo, args.output)
    print(json.dumps({k: result[k] for k in ("paired_relations", "matched_proteins", "annotation_positive_proteins")}))
