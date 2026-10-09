"""Check the new metadata boundary, then reuse independent rational arithmetic."""

import argparse
import json
from pathlib import Path

from benchmark_tools import readback_native_qfo_four_cell_strata as rational


AMENDMENT_SHA = "0edde329061ad42a9a3fe14d0e2a55373ffb285fa3e17a6715a72b4cac1b0fd4"
EXPORT_COMPONENT_SHA = "024092bb653947ca2e6d1498498fe51273cd3718795abb0806367977073eac10"
READER_COMPONENT_SHA = "9fb89edea6ea7f38f2cac25bb237c0cc9d475455d04e361ddab3eb9de60ef842"


def verify(result, docs):
    rational.need(result["schema"] == "native_qfo_four_cell_strata_v2", "Wrong prospective report schema")
    scores = {(r["suite"], r["cell"], r["stratum"]):r for r in result["rows"]}
    differences = {(r["suite"], r["contrast"], r["stratum"]):r for r in result["differences"]}
    changes = []
    for key in ("fixed", "distance"):
        doc = docs[key]
        bins = doc["bins"]
        expected = ({(s, c, n) for s in bins for c in rational.CELLS[:3] for n in bins[s]}
                    if key == "fixed" else
                    {("model_distance", c, n) for c in rational.CELLS[:3] for n in bins})
        seen = set()
        for row in doc["rows"]:
            suite = row["suite"] if key == "fixed" else "model_distance"
            ident = suite, row["cell"], row["stratum"]
            rational.need(ident in expected and ident not in seen and ident in scores, "Wrong legacy row identity")
            seen.add(ident)
            normalized = "resolved_native_pairs" if row["cell"] == rational.CELLS[1] else "group_clique"
            retained = "native_pair" if key == "distance" and row["cell"] == rational.CELLS[1] else normalized
            rational.need(row["prediction_semantics"] == retained
                          and scores[ident]["prediction_semantics"] == normalized, "Legacy vocabulary does not match")
            rational.need(all(row[k] == scores[ident][k] for k in ("families", "family_members", "status")),
                          "Legacy score identity changed")
            if retained != normalized:
                changes.append(dict(input=key, suite=suite, cell=row["cell"], stratum=row["stratum"],
                                    retained=retained, normalized=normalized))
        rational.need(seen == expected, "Incomplete legacy scores")
        expected_diffs = {(s, label, n) for s, _, n in expected for label, _, _ in rational.CONTRASTS[:2]}
        seen_diffs = set()
        for row in doc["differences"]:
            suite = row["suite"] if key == "fixed" else "model_distance"
            ident = suite, row["contrast"], row["stratum"]
            rational.need(ident in expected_diffs and ident not in seen_diffs and ident in differences,
                          "Wrong legacy contrast identity")
            seen_diffs.add(ident)
            rational.need(all(row[k] == differences[ident][k] for k in
                              ("candidate", "reference", "families", "family_members", "status")),
                          "Legacy contrast metadata changed")
        rational.need(seen_diffs == expected_diffs, "Incomplete legacy contrasts")
    rational.need(len(changes) == 3 and result["compatibility"] ==
                  dict(scope="in_memory_metadata_only", rows=changes, row_count=3,
                       original_input_bytes_changed=False, numerical_values_changed=False),
                  "Untruthful compatibility mapping")
    # Only the reused arithmetic helper's schema gate needs this internal view.
    arithmetic_view = dict(result, schema="native_qfo_four_cell_strata_v1")
    return rational.verify(arithmetic_view, docs)


def readback(path, sha, output):
    path, output = Path(path).resolve(), Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    report = rational.record(path)
    rational.need(report["sha256"] == sha, "Changed selected report")
    result = json.loads(path.read_text())
    rational.need(set(result["inputs"]) == set(rational.PINS)
                  and result["protocol"]["sha256"] == rational.PROTOCOL_SHA
                  and result["amendment"]["sha256"] == AMENDMENT_SHA, "Changed input/protocol/amendment pins")
    rational.need(len(result["components"]) == 1 and result["components"][0]["sha256"] == EXPORT_COMPONENT_SHA
                  and result["components"][0] in result["checked_inputs"]
                  and result["amendment"] in result["checked_inputs"], "Unbound compatibility components")
    component = rational.record(rational.__file__)
    rational.need(component["sha256"] == READER_COMPONENT_SHA, "Changed reused rational reader")
    source = rational.record(__file__)
    refs = [report, result["source"], result["protocol"], *result["checked_inputs"],
            *result["outputs"], component, source]
    docs = {}
    for key, digest in rational.PINS.items():
        ref = result["inputs"][key]
        rational.need(ref["sha256"] == digest and ref in result["checked_inputs"], "Changed scientific input pin")
        rational.checked(ref)
        docs[key] = json.loads(Path(ref["path"]).read_text())
        rational.need(docs[key]["source"] in result["checked_inputs"], "Inherited source not checked")
    rational.need(docs["fixed_reader"]["report"] == result["inputs"]["fixed"]
                  and docs["distance_reader"]["report"] == result["inputs"]["distance"]
                  and docs["profile_reader"]["audit"] == result["inputs"]["profile"], "Wrong original reader binding")
    rational.need([Path(r["path"]).name for r in result["outputs"]] ==
                  ["scores.tsv", "differences.tsv", "family_counts.tsv", "TABLE.md"]
                  and all(Path(r["path"]).parent == path.parent for r in result["outputs"]), "Unbound output")
    for ref in refs:
        rational.checked(ref)
    records = verify(result, docs)
    tables = rational.check_tables(path.parent, result, *records)
    for ref in refs:
        rational.checked(ref)
    receipt = dict(schema="native_qfo_four_cell_strata_rational_readback_v2", report=report,
                   source=source, components=[component], checked_inputs=refs,
                   report_schema_checked="native_qfo_four_cell_strata_v2",
                   internal_arithmetic_contract="native_qfo_four_cell_strata_v1",
                   compatibility_rows_checked=3, families_checked=18, proteins_checked=563,
                   family_rows_checked=72, score_rows_checked=92, differences_checked=69,
                   inherited_score_rows_reproduced=69, inherited_difference_rows_reproduced=46,
                   **tables, new_bootstrap_draws=0, independent_confirmation=False, publication_ready=False,
                   new_uncertainty=False, new_accuracy_or_resource_admission=False, scientific_timings_admitted=False,
                   raw_evidence_reparsed=False, original_tree_traversal_repeated=False,
                   limitations=["Explicit metadata adapter plus unchanged independent rational/table kernel.",
                                "Retained-record validation, not repeated raw scoring or transitive admission.",
                                "Overlapping bins are descriptive, not causal, independent or new subgroup intervals."])
    with output.open("x") as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = readback(args.report, args.report_sha256, args.output)
    print(json.dumps({k:result[k] for k in ("family_rows_checked", "score_rows_checked", "differences_checked")}))
