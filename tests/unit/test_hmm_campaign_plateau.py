from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import report_hmm_campaign_plateau as report
from tests.unit.test_weighted_projection_followup_selection import admitted


def evidence():
    freeze = json.loads((Path(__file__).resolve().parents[2] / report.FREEZE).read_text())
    review = dict(status="bounded_postfreeze_confirmation_review_complete_no_independent_panel_admitted",
        freeze_revision=report.FREEZE_REVISION, logical_probes=6,
        request_policy=deepcopy(freeze["postfreeze_confirmation_review"]["request_bounds"]),
        requests=[dict(url=url, method="HEAD", body_bytes_read=0)
                  for url in freeze["postfreeze_confirmation_review"]["head_urls"]],
        defaults_promoted=False, goal_completed=False, independent_confirmation_established=False,
        limitations_nonblocking=True, new_accuracy_evaluation=False,
        new_reference_assignments_read=False, new_reference_assignments_downloaded=False)
    pilot, producer = admitted((False, False))
    pilot["uncertainty"] = dict(multiplicity_endpoints=18, contrasts=[
        dict(candidate=c, control=b, result={"fixture": "unchanged"})
        for c, b in report.specification()["uncertainty"]["contrasts"]])
    controls = {report.configuration(c["label"])["reference_key"]: c for c in pilot["cells"][:4]}
    full = {freeze["baseline"]["key"], *report.CANDIDATES, "projection_recovery:legacy_off"}
    methods = []
    for key in report.INDIVIDUAL_KEYS:
        failed = key.endswith(":failed_panel")
        scores = {} if failed else {k: .7 for k in (report.base.ENDPOINTS if key in full else report.base.ENDPOINTS[:1])}
        control = controls.get(key)
        if control:
            scores[report.base.ENDPOINTS[0]] = control["scores"][report.base.ENDPOINTS[0]]
        details = deepcopy(control["details"]) if control else {"precision": .7, "recall": .7}
        cells = [] if failed else [dict(dataset="OrthoBench", label=key,
            status="independent_content_readback_verified", details=details,
            scores={report.base.ENDPOINTS[0]: scores[report.base.ENDPOINTS[0]]})]
        if key in full:
            qfo = {k: v for k, v in scores.items() if k.startswith("QfO:")}
            cells += [dict(dataset="QfO", label=key, status="independent_content_readback_verified",
                scores=qfo, details={name: dict(similarity=.7, assessed_relations=100) for name in ("EC", "GO", "FAS")},
                functional_populations={name: {"sample_pairs" if name == "FAS" else "numeric_scored_pairs": 40}
                                       for name in ("EC", "GO", "FAS")}),
                dict(dataset="ThreeKingdoms", label=key, status="independent_content_readback_verified",
                     scores={report.base.ENDPOINTS[-1]: .7}, details=dict(f_score=.7, precision=.7, recall=.7))]
        methods.append(dict(key=key, label=key, scores=scores, cells=cells,
            status="complete_development_vector" if key in full else "failed_incomplete_panel_not_accuracy_evidence"
                   if failed else "stopped_after_orthobench_point_gate",
            missing_endpoints=[k for k in report.base.ENDPOINTS if k not in scores]))
        if key in report.CANDIDATES:
            methods[-1]["full_vector_screen"] = {"eligible": True}
    canonical = methods[0]
    context = [dict(key="canonical_collection:canonical", label="Canonical", status="retained_context_not_reaudited_here",
        scores=deepcopy(canonical["scores"]), population_details=deepcopy(canonical["cells"]))]
    individual = dict(status=report.diagnostics.STATUS, accuracy_result_admitted=True,
        defaults_promoted=False, goal_completed=False, context=context, methods=methods,
        uncertainty={"weighting": "retained original intervals"},
        retained_fas_samples={"populations": [dict(mode=mode, FAS_reported_eligible_pool=100,
            FAS_numeric_sample_pairs=40) for mode in ("none", "henikoff", "henikoff_cluster")]})
    pilot["context"] = deepcopy(context)
    return individual, pilot, producer, freeze, review


def test_preserves_all_individual_effects_and_all_six_interaction_roles_without_imputation():
    args = evidence()
    before = deepcopy(args)
    result = report.assemble(*args)
    assert args == before
    assert result["methods"][:14] == args[0]["methods"]
    assert result["uncertainty"] == args[0]["uncertainty"]
    assert result["interaction"]["uncertainty"] == args[1]["uncertainty"]
    assert result["interaction"]["selection"]["execution_arms"] == []
    for method, cell in zip(result["methods"][14:], args[1]["cells"]):
        assert method["cells"] == [cell] and method["scores"] == cell["scores"]
        assert method["role"] == report.configuration(cell["label"])["role"]
        assert method["missing_endpoints"] == list(report.base.ENDPOINTS[1:])
        assert "full_vector_screen" not in method
    assert not result["defaults_promoted"] and not result["goal_completed"]


@pytest.mark.parametrize("kind", ["status", "admission", "default", "complete", "missing_method", "order",
    "context", "baseline", "missing_endpoint", "candidate_screen", "pending_individual", "control_point",
    "contrast", "multiplicity"])
def test_partial_changed_or_unfinished_scientific_evidence_is_rejected(kind):
    individual, pilot, producer, freeze, review = evidence()
    if kind == "status": individual["status"] = "running"
    elif kind == "admission": individual["accuracy_result_admitted"] = False
    elif kind == "default": individual["defaults_promoted"] = True
    elif kind == "complete": individual["goal_completed"] = True
    elif kind == "missing_method": individual["methods"].pop()
    elif kind == "order": individual["methods"].reverse()
    elif kind == "context": pilot["context"] = []
    elif kind == "baseline": individual["methods"][0]["scores"][report.base.ENDPOINTS[1]] += .01
    elif kind == "missing_endpoint": individual["methods"][1]["scores"].pop(report.base.ENDPOINTS[1])
    elif kind == "candidate_screen": individual["methods"][1]["full_vector_screen"]["eligible"] = False
    elif kind == "pending_individual": individual["methods"][6]["status"] = "full_vector_extension_required"
    elif kind == "control_point": individual["methods"][1]["cells"][0]["details"]["precision"] += .01
    elif kind == "contrast": pilot["uncertainty"]["contrasts"].pop()
    else: pilot["uncertainty"]["multiplicity_endpoints"] = 3
    with pytest.raises(ValueError): report.assemble(individual, pilot, producer, freeze, review)


@pytest.mark.parametrize("passing", [(True, False), (False, True), (True, True)])
def test_any_passing_interaction_requires_real_followup_not_negative_closure(passing):
    individual, _, _, freeze, review = evidence()
    pilot, producer = admitted(passing)
    pilot["context"] = deepcopy(individual["context"])
    with pytest.raises(ValueError, match="full-vector follow-up"):
        report.assemble(individual, pilot, producer, freeze, review)


@pytest.mark.parametrize("kind", ["candidate", "policy", "review_status", "revision", "count", "url", "get",
    "body", "bounds", "blocking", "labels", "download", "accuracy", "independent", "default", "complete"])
def test_confirmation_review_cannot_change_candidates_read_labels_or_claim_validation(kind):
    _, _, _, freeze, review = evidence()
    if kind == "candidate": freeze["retained_candidates"].pop()
    elif kind == "policy": freeze["retained_candidates"][0]["profile_projection_policy"] = "residue_preserving"
    elif kind == "review_status": review["status"] = "running"
    elif kind == "revision": review["freeze_revision"] = "other"
    elif kind == "count": review["requests"].pop()
    elif kind == "url": review["requests"][0]["url"] = "other"
    elif kind == "get": review["requests"][0]["method"] = "GET"
    elif kind == "body": review["requests"][0]["body_bytes_read"] = 1
    elif kind == "bounds": review["request_policy"]["tls_verification"] = False
    elif kind == "blocking": review["limitations_nonblocking"] = False
    else:
        field = dict(labels="new_reference_assignments_read", download="new_reference_assignments_downloaded",
            accuracy="new_accuracy_evaluation", independent="independent_confirmation_established",
            default="defaults_promoted", complete="goal_completed")[kind]
        review[field] = True
    with pytest.raises(ValueError): report.validate_review(freeze, review)


def test_actual_diagnostic_helpers_keep_every_score_and_population_pointer_resolvable():
    result = report.assemble(*evidence())
    current, history = report.diagnostics.build(result, "final.json", {"rows": []}, "retained.json")
    assert len(current) == 8 * (len(result["context"]) + len(result["methods"]))
    report.verify_diagnostic_pointers(result, {"rows": []}, [*current, *history], "final.json", "retained.json")
    assert all(r["precision"] is r["recall"] is None for r in current
               if r["statistic"] == "conditional_similarity" and r["score"] is not None)
    corrupted = deepcopy(next(r for r in current if r["score"] is not None))
    corrupted["score"] += .01
    with pytest.raises(ValueError, match="score pointer"):
        report.verify_diagnostic_pointers(result, {"rows": []}, [corrupted], "final.json", "retained.json")


def test_existing_attempt_and_unfrozen_main_are_not_overwritten(tmp_path):
    with pytest.raises(FileExistsError): report.run(tmp_path, tmp_path)
    source = Path(report.__file__).resolve().parents[1]
    with pytest.raises(ValueError, match="frozen committed"):
        report.run(source, tmp_path / "new")
    assert not (tmp_path / "new").exists()
