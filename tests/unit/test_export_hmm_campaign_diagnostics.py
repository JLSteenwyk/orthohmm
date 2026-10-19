from copy import deepcopy
import csv
import json

import pytest

from benchmark_tools import export_hmm_campaign_diagnostics as export
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_export_hillclimb_diagnostics import retained_ob
from tests.unit.test_hmm_transfer_campaign_report import evidence


def fixture():
    weight, local, pilot, producer, following = evidence()
    report = export.campaign.merge_transfer(export.campaign.recovery.collate(weight, local, pilot, producer),
                                            pilot, producer, following)
    for method in report["methods"]:
        for cell in method["cells"]:
            cell.setdefault("status", "independent_content_readback_verified")
            if cell["dataset"] == "ThreeKingdoms":
                cell["details"] = dict(f_score=cell["scores"][export.base.ENDPOINTS[-1]], precision=.7, recall=.7,
                    reference_genes=10, reference_genes_in_prediction=9, reference_gene_coverage=.9,
                    reference_orthogroups=3, true_positive_gene_pairs=7, false_positive_gene_pairs=3,
                    false_negative_gene_pairs=3, predicted_orthogroups=10)
    canonical = report["methods"][0]
    legacy, retained = retained_ob()
    legacy = legacy["methods"][0]
    report["context"] = [{"key": "canonical", "label": "Canonical context", "status": "retained_context_not_reaudited_here",
        "scores": deepcopy(canonical["scores"]),
        "population_details": [{k: deepcopy(v) for k, v in c.items() if k != "scores"} for c in canonical["cells"]]},
        {"key": legacy["key"], "label": legacy["label"], "status": "retained_context_not_reaudited_here",
         "scores": legacy["scores"], "population_details": legacy["cells"]}]
    report["retained_fas_samples"] = {"populations": [{"mode": mode,
        "FAS_reported_eligible_pool": 100 if mode == "none" else 95, "FAS_numeric_sample_pairs": 40}
        for mode in ("none", "henikoff", "henikoff_cluster")]}
    return report, retained


def resolve(value, path):
    for key in path: value = value[key]
    return value


def test_all_endpoint_slots_precisions_populations_and_real_source_pointers_are_preserved():
    report, retained = fixture()
    before = deepcopy(report)
    current, history = export.build(report, "campaign.json", retained, "retained.json")
    assert report == before
    assert len(current) == 8 * (len(report["methods"]) + len(report["context"]))
    for row in [*current, *history]:
        if row.get("score") is None: continue
        assert resolve(report, row["score_source_location"]) == row["score"]
        source = retained if row["detail_source_document"] == "retained.json" else report
        assert isinstance(resolve(source, row["detail_source_location"]), dict)
        if row.get("coverage_source_location"): assert isinstance(resolve(report, row["coverage_source_location"]), dict)
        if row.get("numeric_population_source_location"):
            assert resolve(report, row["numeric_population_source_location"]) == row["numeric_scored_pairs"]
        if row["statistic"] == "conditional_similarity":
            assert row["precision"] is row["recall"] is None
    native = next(r for r in current if r["method"] == "canonical" and r["dataset"] == "OrthoBench")
    assert native["score_source_location"][:2] == ["context", 0]
    assert native["score_source_location"] == ["context", 0, "scores", "OrthoBench: OrthoBench F1"]
    assert native["detail_source_location"][:3] == ["context", 0, "population_details"]
    assert native["evidence_scope"] == "retained_context_not_reaudited"
    of = next(r for r in current if r["method"] == "orthofinder_3_1_5_full" and r["dataset"] == "OrthoBench")
    assert of["precision"] == .9 and of["recall"] == .65 and of["benchmark_reported_gene_coverage"] == .821
    assert "later debug" in of["semantics"] and of["detail_source_document"] == "retained.json"
    weighted = next(r for r in current if r["method"] == "profile_weighting:henikoff" and r["endpoint"] == "FAS")
    assert weighted["assessed_relations"] == 95 and weighted["numeric_scored_pairs"] == 40
    ec = next(r for r in current if r["method"] == "profile_weighting:henikoff" and r["endpoint"] == "EC")
    assert "numeric_scored_pairs" not in ec
    repaired = next(r for r in current if r["method"] == "projection_recovery:projection_off" and r["endpoint"] == "FAS")
    assert repaired["assessed_relations"] == 100 and repaired["numeric_scored_pairs"] == 40
    assert repaired["method_role"] == "advancing_candidate"
    tk = next(r for r in current if r["method"] == "profile_weighting:none" and r["dataset"] == "ThreeKingdoms")
    assert tk["reference_groups"] == 3 and tk["reference_gene_coverage"] == .9
    failed = [r for r in current if r["method"] == "selective_original:failed_panel"]
    assert len(failed) == 8 and all(r["score"] is None for r in failed)
    stopped = [r for r in current if r["method"] == "projection_recovery:projection_hmmer_reselect"]
    assert sum(r["score"] is not None for r in stopped) == 1


def test_missing_precision_remains_na_with_a_resolvable_parent_pointer():
    report, retained = fixture()
    current, _ = export.build(report, "campaign.json", retained, "retained.json")
    row = next(r for r in current if r["method"] == "profile_weighting:none" and r["endpoint"] == "VGNC F")
    assert row["score"] is not None and row["precision"] is row["recall"] is None
    assert resolve(report, row["detail_source_location"]) == report["methods"][0]["cells"][1]["details"]


@pytest.mark.parametrize("kind", ["duplicate", "modes", "pool", "count", "count_type", "functional_precision", "score", "context_score", "retained_prediction"])
def test_changed_counts_ownership_or_functional_definitions_are_rejected(kind):
    report, retained = fixture()
    if kind == "duplicate": report["context"][0]["key"] = report["methods"][0]["key"]
    elif kind == "modes": report["retained_fas_samples"]["populations"].reverse()
    elif kind == "pool": report["retained_fas_samples"]["populations"][1]["FAS_reported_eligible_pool"] += 1
    elif kind == "count": report["retained_fas_samples"]["populations"][1]["FAS_numeric_sample_pairs"] = 96
    elif kind == "count_type": report["retained_fas_samples"]["populations"][1]["FAS_numeric_sample_pairs"] = True
    elif kind == "functional_precision": report["methods"][0]["cells"][1]["details"]["EC"]["precision"] = .9
    elif kind == "score": report["methods"][1]["scores"]["QfO: EC"] += .01
    elif kind == "context_score": report["context"][0]["population_details"][0]["scores"] = {"OrthoBench: OrthoBench F1": 0.}
    else: retained["rows"][0]["orthobench_retained_evidence"]["prediction_provenance"] = {}
    with pytest.raises(ValueError): export.build(report, "campaign.json", retained, "retained.json")


@pytest.mark.parametrize("field,value", [("sample_pairs", -1), ("sample_pairs", 101), ("sample_pairs", float("nan")),
    ("sample_pairs", True), ("numeric_scored_pairs", 101)])
def test_transfer_numeric_samples_cannot_exceed_eligible_counts_or_have_invalid_values(field, value):
    report, retained = fixture()
    repaired = next(m for m in report["methods"] if m["key"] == "projection_recovery:projection_off")
    repaired["cells"][1]["functional_populations"]["FAS" if field == "sample_pairs" else "EC"][field] = value
    with pytest.raises(ValueError): export.build(report, "campaign.json", retained, "retained.json")


@pytest.mark.parametrize("kind", [None, "pending", "revision", "source", "parent", "submitted", "job", "dependency", "cpus", "memory", "wall", "retained"])
def test_admission_requires_immutable_completed_compiler_and_committed_retained_details(tmp_path, monkeypatch, kind):
    retained_path = tmp_path / export.RETAINED
    retained_path.parent.mkdir(parents=True)
    retained_path.write_text("{}\n")
    entry = {"runner_commit": export.COMPILER, "scientific_reader_source": "benchmarks/work/hmm_transfer_campaign_report_v1_20261019",
        "transfer_reader_submission": export.campaign.READBACK_LAUNCH,
        "job": {"submitted": True, "job_id": "24527", "dependency": "afterany:24525", "cpus": 1,
                "memory_allocation": "16G", "wall_limit": "02:00:00"}}
    if kind == "revision": entry["runner_commit"] = "other"
    elif kind == "source": entry["scientific_reader_source"] = "other"
    elif kind == "parent": entry["transfer_reader_submission"] = "other"
    elif kind == "submitted": entry["job"]["submitted"] = False
    elif kind in ("job", "dependency", "cpus", "memory", "wall"):
        key = {"job": "job_id", "dependency": "dependency", "cpus": "cpus", "memory": "memory_allocation", "wall": "wall_limit"}[kind]
        entry["job"][key] = "other"
    calls = []
    def git(command, **kw):
        calls.append(command)
        return (b"changed\n" if kind == "retained" else retained_path.read_bytes()) if command[-1].endswith(export.RETAINED) else json.dumps({"complete_individual_effects_compilation": entry}).encode()
    monkeypatch.setattr(export.subprocess, "check_output", git)
    def read(repo, actual, status):
        assert actual == entry and status == export.STATUS
        if kind == "pending": raise ValueError("Compiler pending")
        return {"resource_table": record(retained_path)}, [], {"reader": {"path": "campaign.json"}}
    monkeypatch.setattr(export.base, "read_admission", read)
    if kind:
        with pytest.raises(ValueError): export.admit(tmp_path)
    else:
        report, retained, pins, provenance = export.admit(tmp_path)
        assert retained == {} and len(pins) == 2
        assert provenance["launch"]["revision"] == export.LAUNCH
        assert provenance["retained"] == record(retained_path)
    assert calls[0] == ["git", "show", export.LAUNCH + ":" + export.base.LEDGER]


@pytest.mark.parametrize("kind", [None, "pending", "changed_pin"])
def test_actual_export_writers_preserve_endpoint_scopes_and_pending_guard(tmp_path, monkeypatch, kind):
    report, retained = fixture()
    source, old = tmp_path / "campaign.json", tmp_path / "retained.json"
    source.write_text(json.dumps(report))
    old.write_text(json.dumps(retained))
    pins = [record(source), record(old)]
    if kind == "changed_pin": source.write_text("changed\n")
    def admit(repo):
        if kind == "pending": raise ValueError("Compiler pending")
        return report, retained, pins, {"reader": pins[0], "retained": pins[1]}
    monkeypatch.setattr(export, "admit", admit)
    monkeypatch.setattr(export.campaign, "reporter_inventory", lambda _: ("exporter", []))
    output = tmp_path / "output"
    if kind:
        with pytest.raises(ValueError): export.run(tmp_path, output)
        assert not output.exists()
    else:
        result = export.run(tmp_path, output)
        assert json.loads((output / "summary.json").read_text()) == result
        assert result["new_scoring_or_sampling"] is result["defaults_promoted"] is result["goal_completed"] is False
        assert result["reporter_revision"] == "exporter"
        for pin in result["tables"]: assert pin == record(pin["path"])
        rows = list(csv.DictReader((output / "endpoint_diagnostics.tsv").open(), delimiter="\t"))
        assert len(rows) == len(result["current"])
        failed = [row for row in rows if row["method"] == "selective_original:failed_panel"]
        assert all(row["score"] == row["precision"] == row["numeric_scored_pairs"] == "NA" for row in failed)
        with pytest.raises(FileExistsError): export.run(tmp_path, output)


def test_existing_export_rejected_before_source_or_admission_reads(tmp_path):
    with pytest.raises(FileExistsError): export.run(tmp_path, tmp_path)
