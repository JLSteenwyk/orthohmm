"""Independent raw/statistic readback and corrupted native-bin receipt refusals."""

from collections import Counter
import csv
import gzip
import json
from pathlib import Path

import pytest

from benchmark_tools import readback_native_qfo_swiss_sequence_strata as reader


def test_actual_retained_projection_and_independent_receipt():
    results = Path(reader.__file__).parent / "results"
    receipt = json.loads((results / "native_qfo_swiss_sequence_strata_readback_20261006.json").read_text())
    assert receipt["source"] == reader.record(reader.__file__)
    verified = reader.verify(Path(receipt["report"]["path"]), receipt["report"]["sha256"])
    assert verified == receipt
    assert receipt["raw_rows_checked"] == 21530 and receipt["family_rows_checked"] == 36
    assert receipt["projection_rows_checked"] == 22 and receipt["differences_checked"] == 11


def test_independent_statistic_is_macro_not_pooled_or_mean_f1():
    counts = dict(a=dict(TP=10, FP=100, FN=0, TN=0), b=dict(TP=100, FP=0, FN=100, TN=0))
    actual = reader.stats(counts, ["a", "b"])
    p = (12 / 114 + 102 / 104) / 2
    r = (12 / 14 + 102 / 204) / 2
    assert actual["F1"] == pytest.approx(2 * p * r / (p + r), abs=1e-12)
    mean_f1 = sum(reader.stats(counts, [f])["F1"] for f in counts) / 2
    pooled_f1 = reader.stats(dict(pooled={k: sum(c[k] for c in counts.values()) for k in counts["a"]}), ["pooled"])["F1"]
    assert abs(actual["F1"] - mean_f1) > .01 and abs(actual["F1"] - pooled_f1) > .01
    assert reader.stats(counts, []) == dict.fromkeys(reader.METRICS)
    assert not reader.equal(True, 1.0) and not reader.equal(float("nan"), .5)


@pytest.mark.parametrize("fault", (None, "header", "duplicate", "label", "membership", "truth", "width"))
def test_independent_raw_enumeration(tmp_path, fault):
    raw = tmp_path / "raw.gz"
    header = reader.HEADER if fault != "header" else "wrong"
    line = "f\ta\tb\t" + ("wrong" if fault == "label" else "TP") + "\n"
    if fault == "width":
        line = "f\ta\tb\n"
    elif fault == "truth":
        line = "f\ta\ta\tTP\n"
    with gzip.open(raw, "wt") as stream:
        stream.write(header + "\n" + line + (line if fault == "duplicate" else ""))
    members = dict(f=["a", "c" if fault == "membership" else "b"])
    if fault is None:
        counts, truth = reader.raw_counts(raw, members)
        assert counts == dict(f=dict(TP=1, FP=0, FN=0, TN=0)) and truth == {("f", "a", "b"): True}
    else:
        with pytest.raises(ValueError):
            reader.raw_counts(raw, members)


@pytest.mark.parametrize("fault", (None, "scope", "source", "stratum_members", "stratum_score", "empty", "semantics",
                                 "difference", "family", "counts", "rows", "relation_count", "table", "columns"))
def test_full_independent_readback_synthetic_reports(tmp_path, monkeypatch, fault):
    def dump(name, value):
        path = tmp_path / name
        path.write_text(json.dumps(value))
        return reader.record(path)

    members = {f"f{i:02d}": [f"g{i:02d}a", f"g{i:02d}b"] for i in range(18)}
    bins = dict(all=sorted(members), lower_entropy=sorted(members)[:9], higher_entropy=sorted(members)[9:], missing_entropy=[],
        concentrated=[], not_concentrated=sorted(members), missing=[], short_relative=sorted(members)[:7],
        not_short_relative=sorted(members)[7:], explicit_fragment=[], no_explicit_fragment=sorted(members))
    features = dict(family_memberships=members, primary_strata={k: bins[k] for k in
        ("lower_entropy", "higher_entropy", "missing_entropy")}, secondary_strata={k: v for k, v in bins.items()
        if k not in ("all", "lower_entropy", "higher_entropy", "missing_entropy")})
    feature_ref = dump("feature.json", features)
    monkeypatch.setattr(reader, "STRATA_SHA", feature_ref["sha256"])
    audits, evidence, counts_by_cell = {}, [], {}
    for index, cell in enumerate(reader.CELLS):
        raw = tmp_path / f"raw{index}.gz"
        counts = {}
        with gzip.open(raw, "wt") as stream:
            stream.write(reader.HEADER + "\n")
            for i, (family, genes) in enumerate(members.items()):
                label = "TP" if (i + index) % 2 else "FN"
                stream.write(f"{family}\t{genes[0]}\t{genes[1]}\t{label}\n")
                counts[family] = Counter({k: int(k == label) for k in ("TP", "FP", "FN", "TN")})
        raw_ref = reader.record(raw)
        audit = dump(f"audit{index}.json", dict(cells=[dict(cell=cell, raw_file=raw_ref,
                                                          aggregate=reader.stats(counts, list(members)))]))
        evidence.extend([raw_ref, audit])
        audits[cell] = dict(count_audit=audit)
        counts_by_cell[cell] = counts
    binding_ref = dump("binding.json", dict(bound_cells=audits))
    rows = [dict(cell=cell, stratum=name, families=len(selected), family_members=selected,
        status="descriptive" if selected else "empty_bin",
        prediction_semantics="group_clique" if cell == reader.CELLS[0] else "resolved_native_pairs",
        **reader.stats(counts_by_cell[cell], selected)) for cell in reader.CELLS for name, selected in bins.items()]
    families = [dict(cell=cell, family=family, counts_without_prior=counts_by_cell[cell][family],
                     **reader.stats(counts_by_cell[cell], [family])) for cell in reader.CELLS for family in members]
    differences = []
    for name, selected in bins.items():
        left, right = (reader.stats(counts_by_cell[cell], selected) for cell in reader.CELLS)
        differences.append(dict(stratum=name, families=len(selected), family_members=selected,
            status="descriptive" if selected else "empty_bin",
            **{k: None if left[k] is None else right[k] - left[k] for k in reader.METRICS}))
    table = tmp_path / "scores.tsv"
    fields = ["cell", "stratum", "families", "status", *reader.METRICS, "prediction_semantics"]
    if fault == "columns":
        fields.remove("PPV")
    with table.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({k: "NA" if v is None else v for k, v in row.items()} for row in rows)
    if fault == "table":
        table.write_text(table.read_text().replace("p0_c0_r0", "wrong", 1))
    markdown = tmp_path / "TABLE.md"
    markdown.write_text("Fixture table, not actual publication output.\n")
    report = dict(schema="native_qfo_swiss_sequence_strata_v1",
        source=reader.record(Path(reader.__file__).with_name("export_native_qfo_swiss_sequence_strata.py")),
        strata=feature_ref, binding=binding_ref, checked_inputs=evidence,
        outputs=[reader.record(table), reader.record(markdown)], rows=rows, differences=differences,
        family_rows=families, raw_relations_checked=36, new_uncertainty=False,
        new_accuracy_or_resource_admission=False, independent_confirmation=False, publication_ready=False)
    if fault == "scope":
        report["new_uncertainty"] = True
    elif fault == "source":
        report["source"] = reader.record(__file__)
    elif fault == "stratum_members":
        rows[0]["family_members"] = []
    elif fault == "stratum_score":
        rows[0]["F1"] += .01
    elif fault == "empty":
        next(r for r in rows if not r["families"])["F1"] = 0.
    elif fault == "semantics":
        rows[0]["prediction_semantics"] = "resolved_native_pairs"
    elif fault == "difference":
        differences[0]["F1"] += .01
    elif fault == "family":
        families[0]["F1"] += .01
    elif fault == "counts":
        families[0]["counts_without_prior"]["TP"] += 1
    elif fault == "rows":
        rows.pop()
    elif fault == "relation_count":
        report["raw_relations_checked"] += 1
    ref = dump("report.json", report)
    if fault is None:
        result = reader.verify(Path(ref["path"]), ref["sha256"])
        assert result["raw_rows_checked"] == 36 and result["projection_rows_checked"] == 22
        assert result["family_rows_checked"] == 36 and result["differences_checked"] == 11
        assert result["new_uncertainty"] is result["publication_ready"] is False
    else:
        with pytest.raises(ValueError):
            reader.verify(Path(ref["path"]), ref["sha256"])
