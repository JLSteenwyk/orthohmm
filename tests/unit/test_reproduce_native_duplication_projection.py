import copy
import csv
import gzip
import json
from pathlib import Path
import shutil

import pytest

from benchmark_tools import reproduce_native_duplication_projection as replay


def write_json(path, value):
    path.write_text(json.dumps(value, sort_keys=True) + "\n")
    return replay.raw.record(path)


def fixture(tmp_path, monkeypatch):
    original, exported = tmp_path / "original", tmp_path / "exported"
    original.mkdir()
    exported.mkdir()
    members = dict(F=["A", "B"], G=["C", "D"])
    count_vectors, cells, audits, records, family_rows = {}, [], {}, [], []
    for i, cell in enumerate(replay.raw.CELLS):
        path = original / (cell + ".gz")
        labels = ("TP", "FP") if i == 0 else ("FN", "TN")
        with gzip.open(path, "wt") as stream:
            stream.write(replay.raw.HEADER + "\n")
            stream.write("F\tA\tB\t" + labels[0] + "\nG\tC\tD\t" + labels[1] + "\n")
        raw_ref = replay.raw.record(path)
        records.append(raw_ref)
        counts, _ = replay.raw.raw_counts(path, members)
        count_vectors[cell] = counts
        family_rows.extend(dict(cell=cell, family=f, counts_without_prior=dict(counts[f]),
                                **replay.raw.stats(counts, [f])) for f in members)
        declaration = dict(cell=cell, raw_file=raw_ref, admission=dict(path="/not/accessed/admission", bytes=0, sha256="0" * 64),
                           native_job_id=100 + i)
        if i:
            declaration.update(timing_eligible=False, timing_admitted=False)
        cells.append(declaration)
        audit = dict(cells=[dict(copy.deepcopy(declaration), resources=None,
                                 aggregate=replay.raw.stats(counts, list(members)))])
        audits[cell] = write_json(original / (cell + "_audit.json"), audit)
        records.append(audits[cell])
    binding_ref = write_json(original / "binding.json", dict(bound_cells={c: dict(count_audit=r) for c, r in audits.items()}))
    native_ref = write_json(original / "native.json", dict(binding=binding_ref, family_rows=family_rows))
    records.extend((binding_ref, native_ref))
    feature = dict(families={}, median_fraction="1/2", primary_strata=dict(
        lower_duplication_fraction=["F"], upper_duplication_fraction=["G"], missing_duplication_fraction=[]))
    for i, family in enumerate(members):
        feature["families"][family] = dict(explicit_duplication_nodes=i, explicit_speciation_nodes=0,
            default_speciation_nodes=1 - i, informative_nodes=1, child_overlap_nodes=0, mapped_members=2,
            duplication_fraction=str(i))
    feature_ref = write_json(original / "feature.json", feature)
    mapping_ref = write_json(original / "mapping.json", dict(families=dict(
        F=dict(mapped_labels=dict(ENSA=1, ENSB=2), exact_match=True, mapped_members=2),
        G=dict(mapped_labels=dict(ENSC=3, ENSD=4, ALIAS=4), exact_match=True, mapped_members=2))))
    identifiers_path = original / "identifiers.gz"
    with gzip.open(identifiers_path, "wt") as stream:
        json.dump(dict(mapping=dict(A=1, B=2, C=3, D=4)), stream)
    records.extend((feature_ref, mapping_ref, replay.raw.record(identifiers_path)))
    bins = replay.reader.independent_bins(members, feature)
    rows = [dict(cell=cell, stratum=name, families=len(values), family_members=values,
                 status="descriptive" if values else "empty_bin",
                 prediction_semantics="group_clique" if cell == replay.raw.CELLS[0] else "resolved_native_pairs",
                 **replay.raw.stats(count_vectors[cell], values))
            for cell in replay.raw.CELLS for name, values in bins.items()]
    differences = [dict(stratum=name, families=len(values), family_members=values,
                        status="descriptive" if values else "empty_bin",
                        **{k: rows[i + 4][k] - rows[i][k] if values else None for k in replay.raw.METRICS})
                   for i, (name, values) in enumerate(bins.items())]
    scores = original / "scores.tsv"
    fields = ["cell", "stratum", "families", "status", *replay.raw.METRICS, "prediction_semantics"]
    with scores.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: "NA" if v is None else v for k, v in row.items()} for row in rows)
    table = original / "TABLE.md"
    table.write_text("Fixture canonical table\n")
    report = dict(schema="native_qfo_swiss_duplication_strata_v1", source=replay.raw.record(
        Path(replay.__file__).with_name("export_native_qfo_swiss_duplication_strata.py")),
        checked_inputs=records, cells=cells, family_rows=family_rows, memberships=members, bins=bins, rows=rows,
        differences=differences, native=native_ref, features=feature_ref, mapping=mapping_ref,
        identifiers=records[-1], outputs=[replay.raw.record(scores), replay.raw.record(table)], raw_relations_checked=4,
        new_bootstrap_draws=0, original_tree_traversal_repeated=False, new_uncertainty=False,
        new_accuracy_or_resource_admission=False, independent_confirmation=False, publication_ready=False)
    projection = original / "report.json"
    projection_ref = write_json(projection, report)
    monkeypatch.setattr(replay, "REPORT_SHA", projection_ref["sha256"])
    binding = replay.relocation.stage(projection, projection_ref["sha256"], "checked_inputs", tmp_path / "relocated")
    for path in (projection, scores, table):
        shutil.copyfile(path, exported / path.name)
    shutil.rmtree(original)
    return report, exported, binding, tmp_path / "result.json"


def test_full_replay_with_original_data_unavailable(tmp_path, monkeypatch):
    report, exported, binding, output = fixture(tmp_path, monkeypatch)
    result = replay.reproduce(exported / "report.json", binding["path"], binding["sha256"],
                              exported / "scores.tsv", exported / "TABLE.md", output)
    assert not (tmp_path / "original").exists()
    assert result["raw_rows_checked"] == 4 and result["family_rows_checked"] == 4
    assert result["projection_rows_checked"] == 8 and result["differences_checked"] == 4
    assert result["original_absolute_inputs_accessed"] is False
    assert result["primary_exporter_rerun"] is False and result["publication_ready"] is False
    assert result["recalculated_rows"][3]["F1"] is None


@pytest.mark.parametrize("target", ["report.json", "scores.tsv", "TABLE.md", "raw", "binding"])
def test_replay_refuses_changed_payloads(tmp_path, monkeypatch, target):
    report, exported, binding, output = fixture(tmp_path, monkeypatch)
    if target == "raw":
        path = Path(binding["path"]).parent / "inputs" / report["cells"][0]["raw_file"]["sha256"]
    elif target == "binding":
        path = Path(binding["path"])
    else:
        path = exported / target
    with path.open("ab") as stream:
        stream.write(b"tampered")
    with pytest.raises((ValueError, json.JSONDecodeError)):
        replay.reproduce(exported / "report.json", binding["path"], binding["sha256"],
                         exported / "scores.tsv", exported / "TABLE.md", output)
    assert not output.exists()


@pytest.mark.parametrize("key,value", [("new_bootstrap_draws", False), ("new_bootstrap_draws", 1),
    ("original_tree_traversal_repeated", True), ("new_uncertainty", True),
    ("new_accuracy_or_resource_admission", True), ("independent_confirmation", True), ("publication_ready", True)])
def test_scope_refusals(key, value):
    report = dict(schema="native_qfo_swiss_duplication_strata_v1", new_bootstrap_draws=0,
                  original_tree_traversal_repeated=False, new_uncertainty=False,
                  new_accuracy_or_resource_admission=False, independent_confirmation=False, publication_ready=False)
    report[key] = value
    with pytest.raises(ValueError):
        replay.check_scope(report)


def test_lookup_has_no_original_path_fallback(tmp_path):
    original = dict(path="/unavailable/original", bytes=1, sha256="0" * 64)
    observed = dict(original, path=str(tmp_path / "payload"))
    locate = replay.relocated_lookup([original, original], [observed, observed])
    assert locate(original) == tmp_path / "payload"
    with pytest.raises(ValueError):
        locate(dict(original, path="/another/original"))
    with pytest.raises(ValueError):
        locate(dict(original, sha256="1" * 64))


@pytest.mark.parametrize("change", ["length", "digest", "conflict"])
def test_lookup_refuses_invalid_bindings(change):
    original = dict(path="/logical", bytes=1, sha256="0" * 64)
    observed = dict(original, path="/relocated")
    if change == "length":
        originals, restored = [original], []
    elif change == "digest":
        originals, restored = [original], [dict(observed, sha256="1" * 64)]
    else:
        originals, restored = [original, original], [observed, dict(observed, path="/other")]
    with pytest.raises(ValueError):
        replay.relocated_lookup(originals, restored)


def test_existing_output_refused_before_input_access(tmp_path):
    with pytest.raises(ValueError, match="Output already exists"):
        replay.reproduce(tmp_path / "missing", tmp_path / "missing", "0" * 64,
                         tmp_path / "missing", tmp_path / "missing", tmp_path)


def test_dangling_output_symlink_refused(tmp_path):
    output = tmp_path / "output"
    output.symlink_to(tmp_path / "missing")
    with pytest.raises(ValueError, match="Output already exists"):
        replay.reproduce(tmp_path / "missing", tmp_path / "missing", "0" * 64,
                         tmp_path / "missing", tmp_path / "missing", output)


def test_changed_executable_helper_refused(tmp_path, monkeypatch):
    _, exported, binding, output = fixture(tmp_path, monkeypatch)
    monkeypatch.setattr(replay, "SOURCE_PINS", {"readback_native_qfo_swiss_duplication_strata.py": "0" * 64})
    with pytest.raises(ValueError, match="Changed replay implementation"):
        replay.reproduce(exported / "report.json", binding["path"], binding["sha256"],
                         exported / "scores.tsv", exported / "TABLE.md", output)
    assert not output.exists()


@pytest.mark.parametrize("change", ["job", "admission", "duplicate", "resources", "aggregate", "timing", "cells"])
def test_native_binding_or_timing_mutations_refused(tmp_path, monkeypatch, change):
    report, _, binding, _ = fixture(tmp_path, monkeypatch)
    document = replay.relocation.frozen_json(binding["path"], binding["sha256"])
    restored, _ = replay.relocation.restore_inputs(report["checked_inputs"], document["source"],
                                                  binding["path"], binding["sha256"])
    locate = replay.relocated_lookup(report["checked_inputs"], restored)
    native = json.loads(locate(report["native"]).read_text())
    original_binding = json.loads(locate(native["binding"]).read_text())
    audit_ref = original_binding["bound_cells"][replay.raw.CELLS[1]]["count_audit"]
    audit = json.loads(locate(audit_ref).read_text())
    selected = audit["cells"][0]
    if change == "job":
        selected["native_job_id"] += 1
    elif change == "admission":
        selected["admission"]["sha256"] = "1" * 64
    elif change == "duplicate":
        audit["cells"].append(copy.deepcopy(selected))
    elif change == "resources":
        selected["resources"] = {}
    elif change == "aggregate":
        selected["aggregate"]["F1"] += .1
    elif change == "timing":
        report["cells"][1]["timing_eligible"] = True
        selected["timing_eligible"] = True
    else:
        report["cells"] = report["cells"][:1]

    def load(ref):
        return audit if ref == audit_ref else json.loads(locate(ref).read_text())

    with pytest.raises(ValueError):
        replay.verify_cells(report, native, original_binding, load, locate)
