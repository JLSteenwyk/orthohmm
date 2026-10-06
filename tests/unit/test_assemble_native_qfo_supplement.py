"""Self-contained report assembly fixtures, not native scientific validation."""

import copy
import csv
import json
from pathlib import Path

import pytest

from benchmark_tools import assemble_native_qfo_supplement as module


@pytest.fixture
def reports(tmp_path):
    source = tmp_path / "inert_source.txt"
    source.write_text("Synthetic identity fixture, never executed.\n")
    source_ref = module.record(source)
    data = {key: dict(schema=schema, publication_ready=False, source=source_ref) for key, schema in module.SCHEMAS.items()}
    rows = []
    for index, cell in enumerate(module.CELLS, 6):
        scores = {name: .6 + .1 * (index - 6) for name in module.ENDPOINTS}
        rows.append(dict(index=index, cell=cell, accuracy_admitted=True, scores=scores,
            secondary_mean=sum(scores.values()) / 6, resources=None, timing_eligible=False, timing_admitted=False))
    rows += [dict(index=i, cell="missing" + str(i), accuracy_admitted=False,
        scores={name: None for name in module.ENDPOINTS}) for i in range(8, 13)]
    data["snapshot"].update(rows=rows, new_scoring_or_admission=False, recovered_inference_resources_admitted=False)
    data["intervals"].update(new_bootstrap_draws=0, replicates_reused=100000, seed_reused=20260922,
        multiplicity_endpoints=42, families=["f" + str(i) for i in range(18)], contrasts=[dict(name="R_at_P0_C0",
            status="native_records_matched", metrics={name: dict(difference=.1, bonferroni_percentile_ci=[-.1, .3],
                family_wins=10, family_ties=3, family_losses=5) for name in ("F1", "PPV", "TPR")})])
    for name in ("vgnc", "graph", "functional", "transitions"):
        data[name].update(uncertainty_admitted=False, new_scoring_or_admission=False)
    data["vgnc"].update(methods=[dict(cell=row["cell"], counts=dict(TP=1, FP=2, FN=3),
        metrics=dict(precision=.1, recall=.2, f1=row["scores"]["VGNC"])) for row in rows[:2]],
        transition_counts=[dict(r0="TP", r1="FN", pairs=1)])
    data["transitions"].update(comparison=dict(changed_relations=2, transitions={"TP->FN": 1, "FP->TN": 1, "FN->TP": 0}))
    data["graph"].update(graph_original_admission_established=False, changed_pairs=2,
        summary=[dict(cell=cell, before=before, graph_support="direct_edge", pairs=1)
            for cell in module.CELLS for before in ("TP", "FP")])
    data["search"].update(changed_pairs=2)
    data["reconciliation"].update(changed_pairs_traced=2, node_annotations_previously_inventoried=False)
    data["strata"].update(differences=[dict(stratum="higher_entropy", families=9, F1=.1, PPV=.2, TPR=-.01),
        dict(stratum="empty", families=0, F1=None, PPV=None, TPR=None)])
    data["functional"].update(comparisons=[dict(metric=name, result={
        "right_" + suffix: 30, "left_" + suffix: 20, "shared_" + suffix: 10,
        "shared_pairs_with_different_serialized_scores": 0}) for name, suffix in
        (("GO", "pairs"), ("EC", "pairs"), ("FAS", "sample_pairs"))])
    for name in ("score_figure", "mechanism_figure"):
        outputs = []
        for kind in ("png", "pdf"):
            path = tmp_path / (name + "." + kind)
            path.write_bytes(b"Synthetic non-rendered asset fixture.")
            outputs.append(module.record(path))
        data[name]["outputs"] = outputs
    # Write dependency order, so each synthetic record is genuinely byte-bound.
    refs = {}

    def bind(name):
        path = tmp_path / (name + ".json")
        path.write_text(json.dumps(data[name]))
        refs[name] = module.record(path)

    bind("snapshot")
    for key in ("intervals", "functional", "transitions", "strata", "vgnc", "score_figure"):
        data[key]["snapshot"] = refs["snapshot"]
    bind("intervals")
    data["strata"]["binding"] = data["score_figure"]["swiss_binding"] = refs["intervals"]
    bind("transitions")
    data["reconciliation"]["transition"] = refs["transitions"]
    bind("reconciliation")
    data["search"]["localization"] = refs["reconciliation"]
    bind("search")
    for key in ("functional", "strata", "vgnc", "score_figure"):
        bind(key)
    for primary in ("functional", "transitions", "reconciliation", "search", "strata", "vgnc", "score_figure"):
        key = primary + "_reader"
        field = "composition" if primary == "functional" else "manifest" if primary.endswith("figure") else "report"
        data[key][field] = refs[primary]
        bind(key)
    data["graph"].update(search=refs["search"], search_readback=refs["search_reader"])
    bind("graph")
    data["graph_reader"]["report"] = refs["graph"]
    bind("graph_reader")
    data["mechanism_figure"]["inputs"] = {key: [refs[key], refs[key + "_reader"]]
        for key in ("reconciliation", "search", "strata")}
    bind("mechanism_figure")
    data["mechanism_figure_reader"]["manifest"] = refs["mechanism_figure"]
    bind("mechanism_figure_reader")
    selected = {k: dict(r, path=str(Path(r["path"]).relative_to(tmp_path))) for k, r in refs.items()}
    selection = tmp_path / "selection.json"
    selection.write_text(json.dumps(dict(schema="native_qfo_supplement_selection_v1", publication_ready=False, inputs=selected)))
    return data, refs, module.record(selection)


def test_complete_synthetic_assembly(tmp_path, reports):
    _, _, selection = reports
    output = tmp_path / "supplement"
    result = module.assemble(tmp_path, selection["path"], selection["sha256"], output)
    assert result["publication_ready"] is False and result["scientific_evidence_replayed"] is False
    assert result["table_rows"] == dict(scores=7, swiss_intervals=3, vgnc_counts=2, vgnc_transitions=1,
        swiss_transitions=3, graph_support=4, strata_differences=2, functional_overlap=3)
    text = (output / "supplement.md").read_text()
    for phrase in ("2 of 7", "not a total-HMM ablation", "42 originally planned endpoints", "not official QfO F1",
        "includes zero", "not a TN", "not validated biological histories",
        "unknown and potentially tool-dependent", "failed R1 timing remains ineligible", "Not submission-ready"):
        assert phrase in text
    with (output / "strata_differences.tsv").open() as stream:
        rows = list(csv.reader(stream, delimiter="\t"))
    assert rows[-1] == ["empty", "0", "Unavailable", "Unavailable", "Unavailable"]
    for ref in result["outputs"]:
        module.check(ref)
    with pytest.raises(ValueError, match="Output already exists"):
        module.assemble(tmp_path, selection["path"], selection["sha256"], output)


@pytest.mark.parametrize("fault", ("cohort", "missing_score", "score_nan", "mean", "timing", "scope", "schema",
    "reader", "snapshot", "chain", "figure_chain", "multiplicity", "seed", "draws", "families", "contrast",
    "vgnc", "changed_pairs", "graph_sum", "disconnected", "graph_admission", "node_admission"))
def test_changed_scientific_contract_refuses(reports, fault):
    data, refs, _ = reports
    data = copy.deepcopy(data)
    if fault == "cohort":
        data["snapshot"]["rows"][0]["cell"] = "wrong"
    elif fault == "missing_score":
        data["snapshot"]["rows"][2]["scores"]["GO"] = .1
    elif fault == "score_nan":
        data["snapshot"]["rows"][0]["scores"]["GO"] = float("nan")
    elif fault == "mean":
        data["snapshot"]["rows"][0]["secondary_mean"] += .1
    elif fault == "timing":
        data["snapshot"]["rows"][1]["timing_eligible"] = True
    elif fault in ("scope", "schema"):
        data["graph"]["publication_ready" if fault == "scope" else "schema"] = True
    elif fault == "reader":
        data["graph_reader"]["report"] = refs["strata"]
    elif fault == "snapshot":
        data["vgnc"]["snapshot"] = refs["strata"]
    elif fault == "chain":
        data["search"]["localization"] = refs["vgnc"]
    elif fault == "figure_chain":
        data["mechanism_figure"]["inputs"]["search"].reverse()
    elif fault in ("multiplicity", "seed", "draws", "families", "contrast"):
        key = {"multiplicity": "multiplicity_endpoints", "seed": "seed_reused", "draws": "new_bootstrap_draws",
            "families": "families", "contrast": "contrasts"}[fault]
        data["intervals"][key] = [] if fault in ("families", "contrast") else 1
    elif fault == "vgnc":
        data["vgnc"]["methods"][0]["metrics"]["f1"] += .1
    elif fault == "changed_pairs":
        data["search"]["changed_pairs"] += 1
    elif fault == "graph_sum":
        data["graph"]["summary"][0]["pairs"] += 1
    elif fault == "disconnected":
        data["graph"]["summary"][0]["graph_support"] = "disconnected"
    elif fault == "graph_admission":
        data["graph"]["graph_original_admission_established"] = True
    else:
        data["reconciliation"]["node_annotations_previously_inventoried"] = True
    with pytest.raises(ValueError):
        module.validate(data, refs)


@pytest.mark.parametrize("fault", ("hash", "input", "unsafe", "missing", "symlink", "outside"))
def test_selection_and_output_refusals(tmp_path, reports, fault):
    _, refs, selection_ref = reports
    path = Path(selection_ref["path"])
    selected = json.loads(path.read_text())
    output = tmp_path / "supplement"
    if fault == "hash":
        selection_ref["sha256"] = "0" * 64
    elif fault == "input":
        Path(refs["graph"]["path"]).write_text("changed")
    elif fault == "unsafe":
        selected["inputs"]["graph"]["path"] = "../escape"
    elif fault == "missing":
        selected["inputs"].pop("graph")
    elif fault == "symlink":
        output.symlink_to(tmp_path / "missing")
    else:
        output = tmp_path.parent / "outside"
    if fault in ("unsafe", "missing"):
        path.write_text(json.dumps(selected))
        selection_ref = module.record(path)
    with pytest.raises(ValueError):
        module.assemble(tmp_path, path, selection_ref["sha256"], output)


def test_table_arithmetic_and_no_family_f1_averaging(reports):
    data, refs, _ = reports
    rows, contrast = module.validate(data, refs)
    values = module.tables(data, rows, contrast)
    assert values["scores"][1][-1][1] == "Secondary, not F1"
    assert values["swiss_intervals"][1][0] == ["F1", 10., -10., 30., 10, 3, 5]
    assert values["graph_support"][1][0] == ["p0_c0_r0", "TP", 1, 0, 0]
    assert len(values["swiss_transitions"][1]) == 3


@pytest.mark.parametrize("value", (False, 0., -1))
def test_draw_count_is_integer_zero_not_a_boolean(reports, value):
    data, refs, _ = reports
    data["intervals"]["new_bootstrap_draws"] = value
    with pytest.raises(ValueError, match="retained uncertainty protocol"):
        module.validate(data, refs)
