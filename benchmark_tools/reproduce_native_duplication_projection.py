"""Replay the native duplication projection using relocated direct inputs only."""

import argparse
import csv
import gzip
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import readback_native_qfo_swiss_duplication_strata as reader
from benchmark_tools import readback_native_qfo_swiss_sequence_strata as raw
from benchmark_tools import relocate_swiss_raw_sources as relocation

REPORT_SHA = "e0ecbc9d26738d0d3d573ca31772085b28d145ba2deaa9383c8cff75bd0d67b5"
SOURCE_PINS = {
    "readback_native_qfo_swiss_duplication_strata.py":
        "1a362d6a406e0bd29e865741d5363b58fd17ba69b923cd8a47ed010b380a1800",
    "readback_native_qfo_swiss_sequence_strata.py":
        "722be2e7d40982f70ecb03aca13b1658fc004a392deacde4cfdb81306a6da396",
    "relocate_swiss_raw_sources.py":
        "5bf4f1e4d2c11fb62b1ba34a95d28d89e8919b15c9c773ad5ca205199fc42d95",
    "prepare_ob_candidate_neighborhood.py":
        "03843ed17fa9c44ea1ce5249782cdf009d8c7db9ceb85c62aa85bc2b4156ec29",
}


def relocated_lookup(originals, restored):
    raw.require(len(originals) == len(restored), "Incomplete relocation")
    pairs = {}
    for original, observed in zip(originals, restored):
        raw.require(all(original[k] == observed[k] for k in ("bytes", "sha256")), "Changed relocated identity")
        previous = pairs.setdefault(original["path"], (original, observed))
        raw.require(previous == (original, observed), "Conflicting repeated relocation")

    def locate(ref):
        raw.require(ref["path"] in pairs and pairs[ref["path"]][0] == ref, "Unbound logical input; no original-path fallback")
        return Path(pairs[ref["path"]][1]["path"])

    return locate


def check_scope(report):
    raw.require(report["schema"] == "native_qfo_swiss_duplication_strata_v1"
                and type(report["new_bootstrap_draws"]) is int and report["new_bootstrap_draws"] == 0
                and all(report[k] is False for k in ("original_tree_traversal_repeated", "new_uncertainty",
                    "new_accuracy_or_resource_admission", "independent_confirmation", "publication_ready")),
                "Wrong original projection scope")


def verify_cells(report, native, binding, load, locate):
    raw.require(report["family_rows"] == native["family_rows"]
                and [c["cell"] for c in report["cells"]] == list(raw.CELLS)
                and set(binding["bound_cells"]) == set(raw.CELLS), "Changed native cell/count inventory")
    counts, truth_anchor, rows = {}, None, 0
    for cell in report["cells"]:
        name = cell["cell"]
        audit = load(binding["bound_cells"][name]["count_audit"])
        selected = [r for r in audit["cells"] if r["cell"] == name]
        raw.require(len(selected) == 1 and all(cell[k] == selected[0][k] for k in cell), "Changed original count binding")
        if name == raw.CELLS[1]:
            raw.require(cell["timing_eligible"] is False and cell["timing_admitted"] is False
                        and selected[0]["resources"] is None, "Failed R1 timing relabeled")
        counts[name], truth = raw.raw_counts(locate(cell["raw_file"]), report["memberships"])
        raw.require(truth_anchor is None or truth == truth_anchor, "Changed native raw truth")
        truth_anchor = truth
        rows += len(truth)
        raw.require(all(raw.equal(selected[0]["aggregate"][k], v)
                        for k, v in raw.stats(counts[name], list(report["memberships"])).items()), "Changed native aggregate")
    raw.require(rows == report["raw_relations_checked"], "Changed native raw row coverage")
    return counts, rows


def verify_tsv(report, path):
    with Path(path).open(newline="") as stream:
        table = csv.DictReader(stream, delimiter="\t")
        fields = ["cell", "stratum", "families", "status", *raw.METRICS, "prediction_semantics"]
        raw.require(table.fieldnames == fields, "Wrong original TSV columns")
        rows = list(table)
    raw.require(rows == [{k: "NA" if row[k] is None else str(row[k]) for k in fields} for row in report["rows"]],
                "Changed original TSV values")


def reproduce(projection, bindings, bindings_sha, scores, table, output):
    output = Path(output)
    raw.require(not output.exists() and not output.is_symlink(), "Output already exists")
    projection_ref = raw.record(projection)
    raw.require(projection_ref["sha256"] == REPORT_SHA, "Changed canonical projection")
    report = json.loads(Path(projection).read_text())
    check_scope(report)
    binding = relocation.frozen_json(bindings, bindings_sha)
    raw.require(all(binding["source"][k] == projection_ref[k] for k in ("bytes", "sha256")), "Wrong relocation source")
    restored, relocation_info = relocation.restore_inputs(report["checked_inputs"], binding["source"], bindings, bindings_sha)
    locate = relocated_lookup(report["checked_inputs"], restored)
    checked = [projection_ref, raw.record(__file__), relocation_info["binding"], *restored]
    for name, sha in SOURCE_PINS.items():
        ref = raw.record(Path(__file__).with_name(name))
        raw.require(ref["sha256"] == sha, "Changed replay implementation: " + name)
        checked.append(ref)
    primary_source = raw.record(Path(__file__).with_name("export_native_qfo_swiss_duplication_strata.py"))
    raw.require(all(primary_source[k] == report["source"][k] for k in ("bytes", "sha256")), "Changed original primary source")
    checked.append(primary_source)

    def load(ref):
        return json.loads(locate(ref).read_text())

    native, feature, mapping = (load(report[k]) for k in ("native", "features", "mapping"))
    feature_binding = load(native["binding"])
    counts, rows = verify_cells(report, native, feature_binding, load, locate)
    with gzip.open(locate(report["identifiers"]), "rt") as stream:
        identifiers = json.load(stream)["mapping"]
    reader.check_memberships(report["memberships"], mapping, identifiers)
    del identifiers
    reader.verify_projection(report, counts, feature)
    for path, name in ((scores, "scores.tsv"), (table, "TABLE.md")):
        original = [r for r in report["outputs"] if Path(r["path"]).name == name]
        ref = raw.record(path)
        raw.require(len(original) == 1 and all(ref[k] == original[0][k] for k in ("bytes", "sha256")), "Changed canonical output")
        checked.append(ref)
    verify_tsv(report, scores)
    recalculated = [dict(cell=cell, stratum=name, **raw.stats(counts[cell], members))
                    for cell in raw.CELLS for name, members in reader.independent_bins(report["memberships"], feature).items()]
    for ref in checked:
        raw.require(raw.record(ref["path"]) == ref, "Relocated replay evidence changed")
    result = dict(schema="relocated_native_duplication_projection_replay_v1", source=raw.record(__file__),
                  original_projection=binding["source"], relocated_projection=projection_ref,
                  relocation=relocation_info, checked_inputs=checked, raw_rows_checked=rows,
                  families_checked=len(report["memberships"]), genes_checked=sum(map(len, report["memberships"].values())),
                  family_rows_checked=len(report["family_rows"]), projection_rows_checked=len(report["rows"]),
                  differences_checked=len(report["differences"]), recalculated_rows=recalculated,
                  original_absolute_inputs_accessed=False, primary_exporter_rerun=False, new_bootstrap_draws=0,
                  inference_reproduced=False, original_tree_traversal_repeated=False,
                  new_accuracy_or_resource_admission=False, independent_confirmation=False,
                  redistribution_authorized=False, publication_ready=False,
                  limitations=["Replay of original independent integer-bin/raw-statistic functions on relocated direct inputs.",
                      "Statistics checked within original1e-12 tolerance; canonical TSV bytes/values checked, not primary JSON regenerated.",
                      "Original logical paths retained as provenance; no original-path lookup or transitive dependency admission.",
                      "No inference, original tree extraction, native scoring, timing repair, new intervals or biological replication.",
                      "Same-host execution, not cross-platform validation, data-rights clearance or final publication archive."])
    with output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("projection", "bindings", "scores", "table", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--bindings-sha256", required=True)
    args = parser.parse_args()
    result = reproduce(args.projection, args.bindings, args.bindings_sha256, args.scores, args.table, args.output)
    print(json.dumps({k: result[k] for k in ("raw_rows_checked", "family_rows_checked", "projection_rows_checked", "differences_checked")}))
