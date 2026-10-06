"""Independent SQLite grouping of native fragment-annotation pair changes."""

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sqlite3
import sys

import Bio
from Bio import SwissProt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.readback_native_qfo_swiss_sequence_strata import require, record, equal

VIEWS = ("historical", "baseline_only")
BINS = ("annotation_positive", "all_matched_unflagged", "missing_without_positive")
LABELS = ("TP", "FP", "FN", "TN")
FIELDS = ("family", "protein_a", "protein_b", "before", "after", "historical", "baseline_only")
HEADER = "# Dataset<tab>Protein ID 1<tab>Protein ID 2<tab>Correctness (TP:True positive, FP: False positive, TN: True Negative. FN: False Negative)"


def annotation_states(a):
    if a is None:
        return None, None
    require(type(a["fragment_flag"]) is bool and isinstance(a["incomplete_sequence_features"], list)
            and a["selection_class"] in ("baseline_release", "later_sequence_version"), "Invalid annotation state")
    positive = int(a["fragment_flag"] or len(a["incomplete_sequence_features"]) > 0)
    return positive, positive if a["selection_class"] == "baseline_release" else None


def sql_projection(before, after, annotations):
    require(before.keys() == after.keys() and set(annotations) == {g for _, a, b in before for g in (a, b)},
            "Changed native pair/annotation universe")
    with sqlite3.connect(":memory:") as db:
        db.execute("CREATE TABLE decisions(family TEXT,a TEXT,b TEXT,before TEXT,after TEXT,PRIMARY KEY(family,a,b))")
        for key, label in before.items():
            other = after[key]
            require(label in LABELS and other in LABELS and (label in ('TP', 'FN')) == (other in ('TP', 'FN')),
                    "Invalid decisions or changed truth")
            db.execute("INSERT INTO decisions VALUES(?,?,?,?,?)", (*key, label, other))
        db.execute("CREATE TABLE annotation(gene TEXT PRIMARY KEY,historical INTEGER,baseline_only INTEGER)")
        db.executemany("INSERT INTO annotation VALUES(?,?,?)", [(g, *annotation_states(a)) for g, a in annotations.items()])
        cases = [f"CASE WHEN l.{v}=1 OR r.{v}=1 THEN 'annotation_positive' "
                 f"WHEN l.{v} IS NULL OR r.{v} IS NULL THEN 'missing_without_positive' "
                 f"ELSE 'all_matched_unflagged' END AS {v}" for v in VIEWS]
        db.execute("CREATE VIEW paired AS SELECT d.*, " + ", ".join(cases) +
                   " FROM decisions d JOIN annotation l ON d.a=l.gene JOIN annotation r ON d.b=r.gene")
        ledger = list(db.execute("SELECT family,a,b,before,after,historical,baseline_only FROM paired ORDER BY family,a,b"))
        require(len(ledger) == len(before), "Lost annotation join rows")
        rows = []
        for view in VIEWS:
            grouped = {(s, a, b): n for s, a, b, n in db.execute(
                f"SELECT {view},before,after,COUNT(*) FROM paired GROUP BY {view},before,after")}
            for name in BINS:
                transitions = {a + "->" + b: grouped.get((name, a, b), 0) for a in LABELS for b in LABELS}
                left = {a: db.execute(f"SELECT COUNT(*) FROM paired WHERE {view}=? AND before=?", (name, a)).fetchone()[0]
                        for a in LABELS}
                right = {b: db.execute(f"SELECT COUNT(*) FROM paired WHERE {view}=? AND after=?", (name, b)).fetchone()[0]
                         for b in LABELS}
                tp, fp = transitions["TP->FN"], transitions["FP->TN"]
                rows.append(dict(view=view, bin=name, relations=sum(transitions.values()), transitions=transitions,
                    before_counts=left, after_counts=right, removed_true_positives=tp, removed_false_positives=fp,
                    tp_removal_fraction=tp / left["TP"] if left["TP"] else None,
                    fp_removal_fraction=fp / left["FP"] if left["FP"] else None))
    return rows, ledger


def parse_labels(path, memberships):
    labels, members = {}, {f: set() for f in memberships}
    with gzip.open(path, "rt", newline="") as stream:
        require(stream.readline().rstrip("\r\n") == HEADER, "Invalid raw header")
        for row in csv.reader(stream, delimiter="\t"):
            require(len(row) == 4, "Invalid native row width")
            f, a, b, label = row
            require(f in memberships and a and b and a != b and label in LABELS, "Invalid native decision")
            key = (f, min(a, b), max(a, b))
            require(key not in labels, "Duplicate canonical native pair")
            labels[key] = label
            members[f].update((a, b))
    require(all(sorted(members[f]) == memberships[f] for f in memberships), "Wrong represented native genes")
    return labels


def verify(path, digest):
    ref = record(path)
    require(ref["sha256"] == digest, "Changed fragment audit")
    report = json.loads(Path(path).read_text())
    require(report["schema"] == "native_qfo_fragment_pair_audit_v1"
            and report["annotation_parser"] == "Bio.SwissProt" and report["biopython_version"] == Bio.__version__
            and report["source"] == record(Path(__file__).with_name("audit_native_qfo_fragment_pairs.py"))
            and report["native"]["sha256"] == "fe003f2cbc4285ea56cd80a92b71703b5c9e135ad1c8f244504dcfeab0911b46"
            and report["annotations"]["sha256"] == "a480f32666bb96de3395213c7cad024c170318ca110c389f15e4756f8052adf8"
            and report["features"]["sha256"] == "912363276fa5456f11aeff9f398fce71aec8d377c8c8eda7b144105945e554f1"
            and type(report["new_bootstrap_draws"]) is int and report["new_bootstrap_draws"] == 0
            and all(report[k] is False for k in ("new_scoring_or_admission", "uncertainty_admitted",
                "scientific_timings_admitted", "independent_confirmation", "publication_ready",
                "source_history_readmission", "unselected_acquisition_records_rechecked")), "Wrong fragment audit source/scope")
    checked = [ref, report["source"], *report["checked_inputs"], *report["outputs"]]
    for item in checked:
        require(record(item["path"]) == item, "Changed direct evidence")
    native = json.loads(Path(report["native"]["path"]).read_text())
    admission = json.loads(Path(report["annotations"]["path"]).read_text())
    features = json.loads(Path(report["features"]["path"]).read_text())
    original_reader = json.loads(Path(report["original_transition_readback"]["path"]).read_text())
    require(report["original_transition_readback"]["sha256"] == "09249d97c6ed3da6c840586552b4f5a72480cc46bbb944bc60124f5382164ffb"
            and original_reader["report"] == report["native"] and original_reader["raw_rows_checked"] == 21530,
            "Changed original transition readback")
    require(report["cells"] == native["cells"] and len(report["cells"]) == 2
            and [r["cell"] for r in report["cells"]] == ["p0_c0_r0", "p0_c0_r1"]
            and report["cells"][1]["resources"] is None and report["cells"][1]["timing_eligible"] is False
            and report["cells"][1]["timing_admitted"] is False, "Wrong native cells or repaired timing")
    members, annotations = admission["families"], admission["annotations"]
    require(members == features["family_memberships"] and len(members) == 18 and len(annotations) == 563
            and admission["annotation_panel_admitted"] is True and admission["matched"] == 563 and admission["missing"] == 0
            and all(r in admission["records"] for r in features["fasta_inputs"]), "Wrong historical/native feature binding")
    entries = [r for r in admission["records"] if Path(r["path"]).name == "entry.txt"]
    require(entries == report["entry_sources"] and len(entries) == 563 and all(r in checked for r in entries),
            "Wrong selected entry source inventory")
    for item in entries:
        gene = Path(item["path"]).parent.name
        a = annotations[gene]
        with open(item["path"]) as stream:
            entry = SwissProt.read(stream)
        flags = [f.strip()[7:] for f in entry.description.split(";") if f.strip().startswith("Flags: ")]
        incomplete = [dict(type=f.type, location=str(f.location), qualifiers=f.qualifiers)
                      for f in entry.features if f.type in ("NON_TER", "NON_CONS")]
        require(entry.accessions[0] == gene and entry.taxonomy_id == [a["taxid"]]
                and entry.sequence_update[1] == a["sequence_version"] and entry.annotation_update[1] == a["entry_version"]
                and hashlib.sha256(entry.sequence.encode("ascii")).hexdigest() == a["sequence_sha256"]
                and len(entry.sequence) == a["length"] and entry.annotation_update[0] == a["annotation_date"]
                and entry.data_class == a["data_class"] and flags == a["flags"]
                and any(f in ("Fragment", "Fragments") for f in flags) == a["fragment_flag"]
                and incomplete == a["incomplete_sequence_features"], "Selected raw annotation disagrees")
    require(report["annotation_positive_proteins"] == sum(annotation_states(a)[0] == 1 for a in annotations.values())
            and report["later_version_proteins"] == sum(a["selection_class"] == "later_sequence_version" for a in annotations.values())
            and report["matched_proteins"] == 563
            and report["source_acquisition_records_inherited"] == sum(r not in report["checked_inputs"] for r in admission["records"]),
            "Wrong annotation scope totals")
    labels = []
    for cell in native["cells"]:
        require(cell["raw"] in checked and cell["admission"] in checked, "Missing native binding")
        scientific = json.loads(Path(cell["admission"]["path"]).read_text())
        require(scientific["accuracy_admitted"] is True and
                scientific["conversion"]["input_fastas"] == features["fasta_inputs"], "Different native input identities")
        labels.append(parse_labels(cell["raw"]["path"], members))
    rows, ledger = sql_projection(*labels, annotations)
    require(len(ledger) == report["paired_relations"] == 10765 and report["raw_rows_checked"] == 21530,
            "Wrong paired relation count")
    require(len(report["rows"]) == len(rows), "Wrong fragment summary row count")
    for actual, expected in zip(report["rows"], rows):
        require(set(actual) == set(expected) and all(actual[k] == v for k, v in expected.items()
                if k not in ("tp_removal_fraction", "fp_removal_fraction"))
                and all(equal(actual[k], expected[k]) for k in ("tp_removal_fraction", "fp_removal_fraction")),
                "Wrong SQL transition/marginal/rate readback")
    for view in VIEWS:
        require(all(sum(r["transitions"][k] for r in rows if r["view"] == view) == n
                    for k, n in native["comparison"]["transitions"].items()), "Original native totals differ")
    output = next(r for r in report["outputs"] if Path(r["path"]).name == "pairs.tsv")
    with open(output["path"], newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        require(next(reader) == list(FIELDS) and [tuple(r) for r in reader] == ledger, "Wrong full paired ledger")
    for item in checked:
        require(record(item["path"]) == item, "Evidence changed during readback")
    return dict(schema="native_qfo_fragment_pair_sql_readback_v1", source=record(__file__), report=ref,
        checked_inputs=checked, annotation_parser="Bio.SwissProt", biopython_version=Bio.__version__,
        shared_annotation_parser=True, entries_checked=563, raw_rows_checked=21530, paired_ledger_rows_checked=10765,
        transitions_checked=96, marginal_count_cells_checked=48, summary_rows_checked=6, new_scoring_or_admission=False,
        uncertainty_admitted=False, scientific_timings_admitted=False, independent_confirmation=False, publication_ready=False,
        limitations=["Independent raw CSV/SQL grouping and annotation decision logic; Bio.SwissProt parser is shared.",
                     "Direct diagnostic checks, not reacquired history admission or independent biological truth.",
                     "No pair-IID inference, causal explanation or timing repair."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists() and not args.output.is_symlink(), "Output already exists")
    result = verify(args.report, args.sha256)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: result[k] for k in ("entries_checked", "raw_rows_checked", "paired_ledger_rows_checked")}))
