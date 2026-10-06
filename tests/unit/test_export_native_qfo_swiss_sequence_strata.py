"""Native frozen-bin projections and source/identity refusals, not new admission."""

import copy
import gzip
import json
import math
from pathlib import Path

import pytest

from benchmark_tools import export_native_qfo_swiss_sequence_strata as exporter
from benchmark_tools.audit_qfo_swiss_counts import HEADER


@pytest.fixture
def sample():
    memberships = {f"f{i:02d}": [f"g{i:02d}a", f"g{i:02d}b"] for i in range(18)}
    genes = {gene: dict(canonical_residues=100, length=100, canonical_entropy_bits=3.0 + i / 50,
                       explicit_fragment_description=False)
             for i, members in enumerate(memberships.values()) for gene in members}
    strata = dict(family_memberships=memberships, genes=genes,
                  **exporter.features.define_strata(memberships, genes))
    cells = []
    for cell in exporter.CELLS:
        rows = []
        for i, (family, members) in enumerate(memberships.items()):
            counts = dict(TP=i + 1, FP=17 - i, FN=i, TN=4)
            if cell.endswith("r1"):
                counts = dict(counts, FP=0, FN=i + 1, TP=i)
            rows.append(dict(family=family, represented_genes=members, counts_without_prior=counts,
                             statistics_with_prior=exporter.family_statistics(counts)))
        aggregate = exporter.aggregate([r["statistics_with_prior"] for r in rows])
        cells.append(dict(cell=cell, families=rows, aggregate=aggregate, native_endpoint_f1=aggregate["F1"]))
    return cells, strata


def test_independent_half_count_and_macro_oracle(sample):
    cells, strata = sample
    rows, differences, family_rows = exporter.project(cells, strata)
    assert len(rows) == 22 and len(differences) == 11 and len(family_rows) == 36
    for row in rows:
        source = next(c for c in cells if c["cell"] == row["cell"])
        selected = [r for r in source["families"] if r["family"] in row["family_members"]]
        if not selected:
            assert row["status"] == "empty_bin" and all(row[k] is None for k in exporter.METRICS)
            continue
        precision, recall, f1s = [], [], []
        for family in selected:
            count = family["counts_without_prior"]
            p = (count["TP"] + 2) / (count["TP"] + count["FP"] + 4)
            r = (count["TP"] + 2) / (count["TP"] + count["FN"] + 4)
            precision.append(p)
            recall.append(r)
            f1s.append(2 * p * r / (p + r))
        p, r = math.fsum(precision) / len(selected), math.fsum(recall) / len(selected)
        assert row["PPV"] == pytest.approx(p, abs=1e-12)
        assert row["TPR"] == pytest.approx(r, abs=1e-12)
        assert row["F1"] == pytest.approx(2 * p * r / (p + r), abs=1e-12)
        if row["stratum"] == "all":
            assert abs(row["F1"] - math.fsum(f1s) / len(f1s)) > 1e-5
    by_key = {(r["cell"], r["stratum"]): r for r in rows}
    for row in differences:
        for metric in exporter.METRICS:
            left, right = (by_key[cell, row["stratum"]][metric] for cell in exporter.CELLS)
            assert row[metric] == (None if left is None else right - left)


@pytest.mark.parametrize("fault", ("cell", "coverage", "members", "counts", "family_stat", "aggregate", "endpoint",
                                 "bin", "feature", "extra_family"))
def test_changed_projection_refuses(sample, fault):
    cells, strata = copy.deepcopy(sample)
    if fault == "cell":
        cells.reverse()
    elif fault == "coverage":
        cells[0]["families"].pop()
    elif fault == "members":
        cells[0]["families"][0]["represented_genes"] = ["wrong"]
    elif fault == "counts":
        cells[0]["families"][0]["counts_without_prior"]["TP"] += 1
    elif fault == "family_stat":
        cells[0]["families"][0]["statistics_with_prior"]["F1"] += .01
    elif fault == "aggregate":
        cells[0]["aggregate"]["F1"] += .01
    elif fault == "endpoint":
        cells[0]["native_endpoint_f1"] += .01
    elif fault == "bin":
        strata["primary_strata"]["lower_entropy"].pop()
    elif fault == "feature":
        strata["genes"]["g00a"]["length"] = 10
    elif fault == "extra_family":
        strata["family_memberships"]["extra"] = ["g00a"]
    with pytest.raises(ValueError):
        exporter.project(cells, strata)


@pytest.mark.parametrize("counts", (dict(TP=True, FP=0, FN=0, TN=0), dict(TP=-1, FP=0, FN=0, TN=0),
                                   dict(TP=1, FP=0, FN=0), dict(TP=1., FP=0, FN=0, TN=0)))
def test_invalid_count_types_refuse(counts):
    with pytest.raises(ValueError, match="counts"):
        exporter.family_statistics(counts)


@pytest.mark.parametrize("fault", (None, "scope", "source", "schema", "input_identity", "admission", "timing", "raw"))
def test_complete_export_synthetic_bindings(tmp_path, monkeypatch, sample, fault):
    def dump(name, value):
        path = tmp_path / name
        path.write_text(json.dumps(value))
        return exporter.record(path)

    cells, strata = copy.deepcopy(sample)
    for cell in cells:
        for row in cell["families"]:
            counts = dict(TP=1, FP=0, FN=0, TN=0)
            row.update(counts_without_prior=counts, statistics_with_prior=exporter.family_statistics(counts))
        cell["aggregate"] = exporter.aggregate([r["statistics_with_prior"] for r in cell["families"]])
        cell["native_endpoint_f1"] = cell["aggregate"]["F1"]
    fasta_inputs = [dict(path="retained-only.fasta", bytes=1, sha256="0" * 64)]
    strata.update(status="corrected_swiss_sequence_strata_prepared_unscored", prediction_statistics_evaluated=False,
        publication_ready=False, source=exporter.record(exporter.features.__file__),
        helper=exporter.record(Path(exporter.__file__).with_name("inventory_swiss_sequences.py")),
        protocol=exporter.record(Path(exporter.__file__).parent / "results/CORRECTED_SWISS_SEQUENCE_STRATA_PROTOCOL_20260918.md"),
        inputs=[], fasta_inputs=fasta_inputs)
    feature_ref = dump("features.json", strata)
    monkeypatch.setattr(exporter, "STRATA_SHA", feature_ref["sha256"])
    plan_ref = dump("plan.json", dict(runs=[{}] * 6 + [dict(inputs=fasta_inputs), dict(inputs=fasta_inputs)]))
    bound, snapshot_rows = {}, []
    retained_ref = dict(path="historical-counts.json", bytes=1, sha256="0" * 64)
    for index, cell in enumerate(cells):
        name = cell["cell"]
        raw = tmp_path / f"raw{index}.gz"
        with gzip.open(raw, "wt") as stream:
            stream.write(HEADER + "\n")
            for family, members in strata["family_memberships"].items():
                stream.write(f"{family}\t{members[0]}\t{members[1]}\t{'FN' if fault == 'raw' and index else 'TP'}\n")
        raw_ref = exporter.record(raw)
        admission = dict(accuracy_admitted=fault != "admission", resources=None, scientific_timings_admitted=False,
            eligible_for_timing_comparison=False, conversion=dict(cell=name, plan=plan_ref, input_fastas=fasta_inputs))
        if fault == "input_identity":
            admission["conversion"]["input_fastas"] = []
        if fault == "timing" and index:
            admission["scientific_timings_admitted"] = True
        admission_ref = dump(f"admission{index}.json", admission)
        cell.update(admission=admission_ref, index=index + 6, native_job_id=index + 10, raw_file=raw_ref)
        if index:
            cell.update(resources=None, timing_admitted=False, timing_eligible=False)
        audit = dict(schema="native_qfo_swiss_family_count_audit_v1" if index == 0 else
                           "recovered_native_qfo_swiss_family_count_audit_v1",
            source=exporter.record(exporter.ordinary.__file__ if index == 0 else exporter.recovered.__file__),
            status="supplied_native_swiss_family_counts_verified" if index == 0 else
                   "supplied_recovered_native_swiss_family_counts_verified", families=sorted(strata["family_memberships"]),
            retained_counts=retained_ref, historical_intervals_attached=False, new_accuracy_or_resource_admission=False,
            independent_confirmation=False, publication_ready=False, cells=[cell], checked_inputs=[raw_ref],
            reference_relation_count=18)
        if fault == "schema":
            audit["schema"] = "wrong"
        audit_ref = dump(f"audit{index}.json", audit)
        bound[name] = {k: cell[k] for k in ("admission", "index", "native_job_id", "native_endpoint_f1")}
        bound[name].update(count_audit=audit_ref, count_aggregate=cell["aggregate"])
        snapshot_rows.append(dict(cell=name, admission=admission_ref, scores=dict(SwissTrees=cell["native_endpoint_f1"]),
            status="supplied_native_admission" if index == 0 else "supplied_recovered_scientific_admission"))
    snapshot_ref = dump("snapshot.json", dict(source=exporter.record(exporter.reporter.__file__),
        schema="native_qfo_scientific_reporting_snapshot_v1", plan=plan_ref, rows=snapshot_rows))
    binding = dict(schema="native_qfo_retained_swiss_uncertainty_binding_v1", source=exporter.record(exporter.binder.__file__),
        bound_cells=bound, snapshot=snapshot_ref, families=sorted(strata["family_memberships"]), retained_counts=retained_ref,
        new_accuracy_or_resource_admission=False, independent_confirmation=False, publication_ready=False)
    if fault == "scope":
        binding["publication_ready"] = True
    elif fault == "source":
        binding["source"] = exporter.record(__file__)
    binding_ref = dump("binding.json", binding)
    output = tmp_path / "export"
    if fault is None:
        result = exporter.export(Path(binding_ref["path"]), binding_ref["sha256"], Path(feature_ref["path"]), output)
        assert len(result["rows"]) == 22 and result["raw_relations_checked"] == 36
        assert result["new_uncertainty"] is False and (output / "TABLE.md").is_file()
        assert all(exporter.record(r["path"]) == r for r in result["outputs"])
        with pytest.raises(ValueError, match="already exists"):
            exporter.export(Path(binding_ref["path"]), binding_ref["sha256"], Path(feature_ref["path"]), output)
    else:
        with pytest.raises(ValueError):
            exporter.export(Path(binding_ref["path"]), binding_ref["sha256"], Path(feature_ref["path"]), output)
        assert not output.exists()
