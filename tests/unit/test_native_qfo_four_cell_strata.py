"""Synthetic subgroup counts test the extension before selected projection."""

import copy
from fractions import Fraction
import json
from pathlib import Path

import pytest

from benchmark_tools import export_native_qfo_four_cell_strata as current
from benchmark_tools import readback_native_qfo_four_cell_strata as independent


ROOT = Path(__file__).resolve().parents[2]


def literal(raw):
    p = Fraction(raw["TP"]+2,raw["TP"]+raw["FP"]+4)
    r = Fraction(raw["TP"]+2,raw["TP"]+raw["FN"]+4)
    return dict(F1=float(2*p*r/(p+r)),PPV=float(p),TPR=float(r))


def mean(records):
    if not records:
        return dict.fromkeys(current.METRICS)
    p = sum(r["PPV"] for r in records)/len(records)
    r = sum(r["TPR"] for r in records)/len(records)
    return dict(F1=2*p*r/(p+r),PPV=p,TPR=r)


def fixture():
    # Invented labels, predictions and bin assignments; fixed universe sizes.
    sizes = (27,32,14,53,63,28,17,15,56,13,50,34,35,12,43,12,31,28)
    families = [f"F{i:02d}" for i in range(18)]
    members = {f:[f"{f}_g{j:03d}" for j in range(n)] for f,n in zip(families,sizes)}
    even,odd = families[::2],families[1::2]
    bins = dict(sequence=dict(all=families,concentrated=[],explicit_fragment=[],higher_entropy=even,
        lower_entropy=odd,missing=[],missing_entropy=[],no_explicit_fragment=families,
        not_concentrated=families,not_short_relative=families[7:],short_relative=families[:7]),
        domain=dict(all=families,median_pfam_types_at_least_two=families[:6],median_pfam_types_below_two=families[6:],
                    repeated_type_fraction_at_least_quarter=families[:3],repeated_type_fraction_below_quarter=families[3:]),
        duplication=dict(all=families,lower_duplication_fraction=even,upper_duplication_fraction=odd,missing_duplication_fraction=[]),
        model_distance=dict(all=families,lower_or_equal_median=families[:9],higher_than_median=families[9:]))
    count_rows = []
    for c,cell in enumerate(current.CELLS):
        for i,(family,n) in enumerate(zip(families,sizes)):
            total = n*(n-1)//2
            positive = total//2
            tp = positive-((i+1)*(c+2) % max(2,positive//2))
            fp = (i+3*c) % max(2,(total-positive)//3)
            raw = dict(TP=tp,FP=fp,FN=positive-tp,TN=total-positive-fp)
            count_rows.append(dict(cell=cell,family=family,counts_without_prior=raw,**literal(raw)))
    values = {(r["cell"],r["family"]):r for r in count_rows}
    rawref = lambda cell: dict(path="/synthetic/"+cell,bytes=1,sha256="0"*64)
    cells = [dict(cell=cell,native_job_id=job,admission=rawref("admission_"+cell),raw_file=rawref("raw_"+cell))
             for cell,job in zip(current.CELLS,(22435,22437,22444,23902))]
    cells[1].update(timing_eligible=False,timing_admitted=False)
    cells[3]["index"] = 10
    common = dict(publication_ready=False,independent_confirmation=False,new_accuracy_or_resource_admission=False,new_bootstrap_draws=0)
    fixed = dict(common,schema="native_qfo_three_cell_strata_v1",memberships=members,cells=copy.deepcopy(cells[:3]),
        bins={k:v for k,v in bins.items() if k != "model_distance"},family_rows=copy.deepcopy(count_rows[:54]),rows=[],differences=[])
    distance = dict(common,schema="native_qfo_swiss_model_divergence_strata_v1",memberships=members,
        cells=copy.deepcopy(cells[:3]),bins=bins["model_distance"],family_rows=copy.deepcopy(count_rows[:54]),rows=[],differences=[])
    for suite in current.SUITES:
        target = distance if suite == "model_distance" else fixed
        scores = {}
        for cell in current.CELLS[:3]:
            for name,group in sorted(bins[suite].items()):
                score = mean([values[cell,f] for f in group])
                scores[cell,name] = score
                row = dict(cell=cell,stratum=name,families=len(group),family_members=group,
                    status="descriptive" if group else "empty_bin",
                    prediction_semantics="resolved_native_pairs" if cell.endswith("r1") else "group_clique",**score)
                if target is fixed:
                    row["suite"] = suite
                target["rows"].append(row)
        for label,candidate,reference in current.CONTRASTS[:2]:
            for name,group in sorted(bins[suite].items()):
                row = dict(contrast=label,candidate=candidate,reference=reference,stratum=name,
                    families=len(group),family_members=group,status="descriptive" if group else "empty_bin")
                for m in current.METRICS:
                    change = scores[candidate,name][m]-scores[reference,name][m] if group else None
                    row[m if target is fixed else m+"_pp"] = change if target is fixed or change is None else 100*change
                if target is fixed:
                    row["suite"] = suite
                target["differences"].append(row)
    profile = dict(copy.deepcopy(cells[3]),aggregate=mean([values[current.CELLS[3],f] for f in families]),families=[])
    for row in count_rows[54:]:
        profile["families"].append(dict(family=row["family"],represented_genes=members[row["family"]],
            counts_without_prior=copy.deepcopy(row["counts_without_prior"]),statistics_with_prior={m:row[m] for m in current.METRICS}))
    audit = dict(common,schema="allocated_native_qfo_swiss_family_count_audit_v1",families=families,
        selected_index=10,reference_relation_count=10765,status="supplied_allocated_native_swiss_family_counts_verified",
        checked_inputs=[profile["raw_file"]],cells=[profile])
    readers = {}
    for key,schema,scores,diffs in (("fixed_reader","native_qfo_three_cell_strata_rational_readback_v2",60,40),
            ("distance_reader","swiss_model_divergence_strata_rational_readback_v1",9,6)):
        readers[key] = dict(common,schema=schema,families_checked=18,proteins_checked=563,family_rows_checked=54,
            score_rows_checked=scores,differences_checked=diffs)
    readers["profile_reader"] = dict(common,schema="allocated_native_qfo_profile_swiss_rational_readback_v1",
        cells=[current.CELLS[1],current.CELLS[3]],native_family_records_checked=36,profile_pair_labels_matched=10765,
        rational_macro_points={c:mean([values[c,f] for f in families]) for c in (current.CELLS[1],current.CELLS[3])})
    snapshot = dict(schema="native_qfo_terminal_failure_reporting_v1",publication_ready=False,rows=[])
    for index,cell in zip(range(6,13),("p0_c0_r0","p0_c0_r1","p0_c1_r0","p0_c1_r1","p1_c0_r1","p1_c1_r0","p1_c1_r1")):
        source = next((r for r in cells if r["cell"] == cell),None)
        row = dict(cell=cell,index=index,accuracy_admitted=source is not None,
            native_job_id=source["native_job_id"] if source else None,admission=copy.deepcopy(source["admission"]) if source else None,
            scores=dict(SwissTrees=mean([values[cell,f] for f in families])["F1"]) if source else dict(SwissTrees=None))
        if cell == current.CELLS[1]:
            row.update(resources=None,timing_eligible=False,timing_admitted=False)
        if cell == current.CELLS[3]:
            row["scientific_timings_admitted"] = False
        snapshot["rows"].append(row)
    return dict(fixed=fixed,distance=distance,profile=audit,snapshot=snapshot,**readers)


def projected(docs):
    return dict(schema="native_qfo_four_cell_strata_v1",**current.SCOPE,**current.project(docs))


def test_synthetic_projection_complete_and_literal_rational_readback():
    docs = fixture()
    before = copy.deepcopy(docs)
    result = projected(docs)
    assert docs == before
    assert [len(result[k]) for k in ("family_rows","rows","differences")] == [72,92,69]
    assert result["prior_score_rows_reproduced"] == 69 and result["prior_difference_rows_reproduced"] == 46
    values,scores,differences = independent.verify(result,docs)
    assert len(values) == 72 and len(scores) == 92 and len(differences) == 69
    for row in result["rows"]:
        if row["families"] == 0:
            assert all(row[m] is None for m in current.METRICS)
    raw = {k:sum(r["counts_without_prior"][k] for r in result["family_rows"] if r["cell"] == current.CELLS[0]) for k in current.COUNTS}
    assert abs(scores["sequence",current.CELLS[0],"all"]["F1"]-Fraction(literal(raw)["F1"])) > 1e-5
    for row in result["differences"]:
        if row["contrast"] == "P_at_C0_R1":
            assert row["reference"] == current.CELLS[1] and row["candidate"] == current.CELLS[3]


@pytest.mark.parametrize("change", ["schema","ready","draws","reader","admission","job","timing",
    "members","truth","negative","boolean","prior_stat","prior_unit","duplicate_prior","bin"])
def test_exporter_refuses_changed_counts_scope_or_bindings(change):
    docs = fixture()
    profile = docs["profile"]["cells"][0]
    if change == "schema": docs["distance"]["schema"] = "other"
    elif change == "ready": docs["profile"]["publication_ready"] = True
    elif change == "draws": docs["profile"]["new_bootstrap_draws"] = True
    elif change == "reader": docs["fixed_reader"]["score_rows_checked"] -= 1
    elif change == "admission": profile["admission"]["sha256"] = "1"*64
    elif change == "job": profile["native_job_id"] += 1
    elif change == "timing": docs["snapshot"]["rows"][1]["timing_admitted"] = True
    elif change == "members": profile["families"][0]["represented_genes"] = ["missing"]
    elif change == "truth": profile["families"][0]["counts_without_prior"]["FN"] += 1
    elif change == "negative": profile["families"][0]["counts_without_prior"]["FP"] = -1
    elif change == "boolean": profile["families"][0]["counts_without_prior"]["TP"] = True
    elif change == "prior_stat": docs["fixed"]["rows"][0]["F1"] += .01
    elif change == "prior_unit": docs["distance"]["differences"][0]["F1_pp"] /= 100
    elif change == "duplicate_prior": docs["fixed"]["rows"][1] = copy.deepcopy(docs["fixed"]["rows"][0])
    elif change == "bin": docs["fixed"]["bins"]["sequence"]["lower_entropy"] = []
    with pytest.raises(ValueError): current.project(docs)


@pytest.mark.parametrize("change", ["scope","cells","bin","family","count","score","empty","difference",
    "reference","duplicate","semantics","units","nan"])
def test_independent_reader_refuses_tampered_projection(change):
    docs = fixture()
    result = projected(docs)
    if change == "scope": result["publication_ready"] = True
    elif change == "cells": result["cells"][1]["timing_eligible"] = True
    elif change == "bin": result["bins"]["model_distance"]["lower_or_equal_median"] = []
    elif change == "family": result["family_rows"].pop()
    elif change == "count": result["family_rows"][0]["counts_without_prior"]["TP"] += 1
    elif change == "score": result["rows"][0]["F1"] += .01
    elif change == "empty": next(r for r in result["rows"] if r["families"] == 0)["F1"] = 0
    elif change == "difference": result["differences"][0]["TPR"] += .01
    elif change == "reference": next(r for r in result["differences"] if r["contrast"] == "P_at_C0_R1")["reference"] = current.CELLS[0]
    elif change == "duplicate": result["rows"][1] = copy.deepcopy(result["rows"][0])
    elif change == "semantics": result["rows"][0]["prediction_semantics"] = "native"
    elif change == "units": result["rows"][0]["F1"] *= 100
    elif change == "nan": result["rows"][0]["F1"] = float("nan")
    with pytest.raises(ValueError): independent.verify(result,docs)


def export_fixture(tmp_path,monkeypatch):
    docs = fixture()
    ref = current.record(__file__)
    monkeypatch.setattr(current,"prepare",lambda root:(docs,dict(fixture=ref),[ref],ref))
    output = tmp_path/"synthetic"
    return output,current.export(ROOT,output),docs


def test_synthetic_serialization_and_all_tables(tmp_path,monkeypatch):
    output,result,docs = export_fixture(tmp_path,monkeypatch)
    records = independent.verify(result,docs)
    counts = independent.check_tables(output,result,*records)
    assert counts == dict(score_tsv_rows=92,difference_tsv_rows=69,family_tsv_rows=72,human_rows_checked=23)
    assert json.loads((output/"report.json").read_text())["schema"] == "native_qfo_four_cell_strata_v1"


@pytest.mark.parametrize("name",["scores.tsv","differences.tsv","family_counts.tsv","TABLE.md"])
def test_table_tampering_detected(tmp_path,monkeypatch,name):
    output,result,docs = export_fixture(tmp_path,monkeypatch)
    path = output/name
    path.write_text(path.read_text().replace("F00","unknown") if name == "family_counts.tsv" else path.read_text().replace("all","tampered",1))
    with pytest.raises(ValueError): independent.check_tables(output,result,*independent.verify(result,docs))


def test_prepare_real_metadata_only_without_selected_subgroups():
    docs,refs,checked,protocol = current.prepare(ROOT)
    assert len(refs) == 7 and protocol["sha256"] == current.PROTOCOL_SHA
    members,rows,_,sources = current.family_values(docs)
    assert len(members) == 18 and len(rows) == 72 and len(sources) == 4
    assert len(checked) >= 9


def test_existing_output_guards_precede_loading(tmp_path,monkeypatch):
    monkeypatch.setattr(current,"prepare",lambda root:pytest.fail("must not load"))
    with pytest.raises(FileExistsError): current.export(ROOT,tmp_path)
    with pytest.raises(FileExistsError): independent.readback(tmp_path/"absent.json","0"*64,tmp_path)


def fully_bound_synthetic_export(tmp_path,monkeypatch):
    docs = fixture()
    source = current.record(__file__)
    refs = {}
    for key in ("fixed","distance","profile","snapshot","fixed_reader","distance_reader","profile_reader"):
        docs[key]["source"] = source
        if key == "fixed_reader": docs[key]["report"] = refs["fixed"]
        if key == "distance_reader": docs[key]["report"] = refs["distance"]
        if key == "profile_reader": docs[key]["audit"] = refs["profile"]
        path = tmp_path/("invented_"+key+".json")
        path.write_text(json.dumps(docs[key],sort_keys=True,allow_nan=False))
        refs[key] = current.record(path)
    protocol = current.record(ROOT/"benchmark_tools/results"/current.PROTOCOL)
    checked = [*refs.values(),source,protocol]
    monkeypatch.setattr(current,"prepare",lambda root:(docs,refs,checked,protocol))
    monkeypatch.setattr(independent,"PINS",{key:ref["sha256"] for key,ref in refs.items()})
    output = tmp_path/"invented_projection"
    current.export(ROOT,output)
    return output,current.record(output/"report.json")


def test_complete_readback_entrypoint_on_invented_bound_artifacts(tmp_path,monkeypatch):
    output,report = fully_bound_synthetic_export(tmp_path,monkeypatch)
    receipt_path = tmp_path/"invented_readback.json"
    receipt = independent.readback(output/"report.json",report["sha256"],receipt_path)
    assert receipt["report"] == report and receipt["family_rows_checked"] == 72
    assert receipt["score_rows_checked"] == 92 and receipt["differences_checked"] == 69
    assert receipt["human_rows_checked"] == 23 and receipt["publication_ready"] is False
    assert json.loads(receipt_path.read_text())["schema"] == "native_qfo_four_cell_strata_rational_readback_v1"


@pytest.mark.parametrize("change",["report_sha","bound_input","output"])
def test_readback_provenance_or_output_tampering_refused(tmp_path,monkeypatch,change):
    output,report = fully_bound_synthetic_export(tmp_path,monkeypatch)
    if change == "report_sha": report["sha256"] = "0"*64
    elif change == "bound_input":
        path = tmp_path/"invented_profile.json"
        path.write_text(path.read_text()+" ")
    else:
        path = output/"scores.tsv"
        path.write_text(path.read_text()+" ")
    receipt_path = tmp_path/"refused.json"
    with pytest.raises(ValueError): independent.readback(output/"report.json",report["sha256"],receipt_path)
    assert not receipt_path.exists()
